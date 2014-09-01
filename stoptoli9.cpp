#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <sstream>
#include <iostream>

#include "TTree.h"
#include "TFile.h"
#include "TError.h"

struct dataparts{
  bool coinov;
  int run, trgId;
  float ctmqtqall, ctrmsts;
  float fido_qiv, fido_qid;
  int fido_nidtubes, fido_nivtubes;

  int ids_didfit;
  float ids_chi2;
  float ids_end_x, ids_end_y, ids_end_z;
  float ids_entr_x, ids_entr_y, ids_entr_z;
  float ids_gclen;
  float ids_ivlen, ids_buflen;

  double deltaT;
  double trgtime;
  float ctX[3];
  float ctEvisID;
  float qrms, qdiff;
};

const unsigned int noff = 95;
const char * turnoff[noff] = {
"ctaplanarity",
"ctFlagMu",
"ctfwhm",
"ctgoodness",
"ctIDMuDeltaT",
"ctIVMuDeltaT",
"ctlightflux",
"ctmqtq",
"ctmqtqflag",
"ctnbadch",
"ctnbadchIV",
"ctnpe",
"ctnpeIV",
"ctnpulse",
"ctnpulseIV",
"ctphi",
"ctq",
"ctqIV",
"ctqtot",
"ctqtotIV",
"ctR",
"ctrho",
"ctsphericity",
"ctt2tot",
"cttmean",
"cttpeak",
"cttrise",
"ctXmuInGC",
"ctXmuInIV",
"ctXmuOuIV",
"date",
"fido_buflen",
"fido_chi2",
"fido_didfit",
"fido_endx",
"fido_endy",
"fido_endz",
"fido_entrx",
"fido_entry",
"fido_entrz",
"fido_gclen",
"fido_ivlen",
"fido_minuit_happiness",
"fido_phi",
"fido_targlen",
"fido_th",
"fido_used_ov",
"hamphi",
"hamth",
"hamx",
"hamxe",
"HEMuDeltaT",
"IVX",
"lilike",
"nev",
"nhit",
"nhitIV",
"novhit",
"novloxy",
"novtrk",
"novupxy",
"ovbadtrk",
"ovloxylike",
"ovloxyx",
"ovloxyy",
"ovloxyz",
"ovtightloxy",
"ovtighttrk",
"ovtightupxy",
"ovtrigid",
"ovtrklike",
"ovtrkphi",
"ovtrkth",
"ovtrkx",
"ovtrky",
"ovupxylike",
"ovupxyx",
"ovupxyy",
"ovupxyz",
"pmtmultpe",
"pmtmultpe_IV",
"timeid",
"timeiv",
"tref",
"trefextIV",
"trefIV",
"trgId",
"trgWord",
"ttovtrig",
"vctnpulse",
"vctnpulseIV",
"vctq",
"vctqIV",
"vcttime",
"vcttimeIV",
};

static bool lightnoise(const float qrms, const float mqtq,
                       const float rmsts, const float qdiff)
{
  if(mqtq > 0.12 || mqtq < 0) return true;
  if(qdiff > 30e3) return true;
  if(rmsts >= 36 && qrms >= 464 - 8*rmsts) return true;
  return false;
}

static void stopper_search(dataparts & parts, TTree * const ctree,
                           TTree * const ftree, const int prompt)
{
  if(prompt >= ctree->GetEntries()){
    fprintf(stderr,
            "Prompt event %d given, but file has only %ld events\n",
            prompt, (long int)ctree->GetEntries());
    return;
  }

  ctree->GetEntry(prompt);
  const double prompttime = parts.trgtime;
  const double li9x=parts.ctX[0],
               li9y=parts.ctX[1],
               li9z=parts.ctX[2];

  for(int back = 1; ; back++){
    if(prompt-back < 0) return;
    ctree->GetEvent(prompt-back);
    ftree->GetEvent(prompt-back);

    // Must be possible that it's a stopper
    if(!parts.ids_didfit) continue;

    // Within possible muons, may want to clean up, but not sure how yet
    /*if(parts.fido_qiv < 5000) continue;
      if(parts.fido_qid/8300 > 700) continue;
      if(parts.fido_chi2/(parts.fido_nidtubes+parts.fido_nivtubes-6)>10)
        continue;
    */

    // Open up a big window
    if(prompttime - parts.trgtime > 10e9) return;

    const double mux = parts.ids_end_x, muy = parts.ids_end_y,
                 muz = parts.ids_end_z;

    const double li9tomu = sqrt(pow(mux-li9x, 2)
                               +pow(muy-li9y, 2)
                               +pow(muz-li9z, 2));

    // And they must be near each other.
    if(li9tomu > 400) continue;

    const double mutime = parts.trgtime;

    ctree->GetEvent(prompt-back+1);

    const float miche = parts.deltaT < 5500? parts.ctEvisID:0.;

    unsigned int nneutron = 0;
    for(back--; back > 0; back--){
      ctree->GetEvent(prompt-back);

      // Not more than ~4 nH lifetimes
      if(parts.trgtime - mutime > 800e3) break;

      // pass light noise
      if(lightnoise(parts.qrms, parts.ctmqtqall,
                    parts.ctrmsts, parts.qdiff)) continue;

      // right energy for a neutron capture
      if(!((parts.ctEvisID > 1.8 && parts.ctEvisID < 2.6) ||
           (parts.ctEvisID > 4.0 && parts.ctEvisID < 10 )))
        continue;

      // near the point the muon stopped (~97% efficient - doc4450)
      if(sqrt(pow(mux - parts.ctX[0], 2)
             +pow(muy - parts.ctX[1], 2)
             +pow(muz - parts.ctX[2], 2)) > 1000) continue;

      nneutron++;
    }

    printf("Muon for %d %d is %d dt %lf ms distance %f with "
           "%u neutrons at %f %f %f michel_e %f\n",
           parts.run, prompt, prompt-back,
           (prompttime - parts.trgtime)/1e6, li9tomu, nneutron,
           li9x, li9y, li9z, miche);
    fflush(stdout);

    break;
  }
}

int main()
{
  std::string line;
  while(std::getline(std::cin, line)){
    std::stringstream ss(line);
    int run;
    if(!(ss >> run)){
      fprintf(stderr, "Could not parse line: %s\n", line.c_str());
      break;
    }

    TFile * const cinfile = new TFile(Form("/cp/s4/dchooz/cheetah/"
      "prod-08-05_p01_v2/reduced.Run%07d_Seq010.root", run), "read");

    TFile * const finfile = new TFile(Form("/cp/s4/strait/fido_seq010/"
      "fido.%07d.root", run), "read");

    TTree * ctree = NULL, * ftree = NULL;

    if(!cinfile || cinfile->IsZombie()){
      fprintf(stderr, "I couldn't read cheetah run %d\n", run);
      goto cleanup;
    }

    if(!finfile || finfile->IsZombie()){
      fprintf(stderr, "I couldn't read fido run %d\n", run);
      goto cleanup;
    }

    ctree = (TTree *)cinfile->Get("data");
    if(!ctree){
      fprintf(stderr, "Run %d lacks a data tree!\n", run);
      goto cleanup;
    }

    ftree = (TTree *)finfile->Get("RecoMuonFIDOInfoTree");
    if(!ftree){
      fprintf(stderr, "Run %d lacks a RecoMuonFIDOInfoTree!\n", run);
      goto cleanup;
    }

    if(ctree->GetEntries() != ftree->GetEntries()){
      fprintf(stderr,
              "Run %d: cheetah has %ld entries, but fido has %ld\n",
              run,
              long(ctree->GetEntries()), long(ftree->GetEntries()));
      goto cleanup;
    }

    ftree->SetMakeClass(1);

    for(unsigned int i = 0; i < noff; i++)
      ctree->SetBranchStatus(turnoff[i], 0);

    dataparts parts;
    #define SBA(x) ctree->SetBranchAddress(#x, &parts.x);
    #define fSBA(x) ftree->SetBranchStatus(#x, 1); \
                    ftree->SetBranchAddress(#x, &parts.x);
    fSBA(ids_didfit);
    fSBA(ids_end_x);
    fSBA(ids_end_y);
    fSBA(ids_end_z);
    fSBA(ids_entr_x);
    fSBA(ids_entr_y);
    fSBA(ids_entr_z);
    fSBA(ids_gclen);
    fSBA(ids_chi2);
    fSBA(ids_ivlen);
    fSBA(ids_buflen);

    SBA(qdiff);
    SBA(ctEvisID);
    SBA(fido_qiv);
    SBA(fido_qid);
    SBA(fido_nidtubes);
    SBA(fido_nivtubes);
    SBA(deltaT);
    SBA(trgtime);
    SBA(ctmqtqall);
    SBA(ctrmsts);
    SBA(qrms);
    SBA(run);
    if(ctree->SetBranchAddress("ctX", parts.ctX) < 0)
      ctree->SetBranchAddress("ctX[3]", parts.ctX);

    int prompt;
    while(ss >> prompt) stopper_search(parts, ctree, ftree, prompt);

    cleanup:

    if(ctree) delete ctree;
    if(ftree) delete ftree;
    if(cinfile) delete cinfile;
    if(finfile) delete finfile;

  }

}
