#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <sstream>
#include <iostream>

#include "TTree.h"
#include "TFile.h"
#include "TError.h"

#include "search.h"

static const double window = 30e9;

static bool earlymich(const int run, const int mutrig)
{
  static TFile *f=TFile::Open("/cp/s4/strait/feb25.t.withmichonly.root");
  static TTree *pulse = (TTree *)f->Get("t");

  // Fantastically inefficient since we know we'll be reading in order,
  // but much easier to code than going sequentially.
  return
    pulse->GetEntries(Form("run == %d && trig == %d", run, mutrig));
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

  // Can't look back all the way
  if(prompttime < window) return;

  const double li9x=bamacorrxy(parts.ctX[0], parts.ctEvisID),
               li9y=bamacorrxy(parts.ctX[1], parts.ctEvisID),
               li9z=bamacorrz( parts.ctX[2], parts.ctEvisID);

  for(int back = 1; ; back++){
    if(prompt-back < 0) return;
    ctree->GetEvent(prompt-back);
    ftree->GetEvent(prompt-back);

    // Must be possible that it's a stopper
    if(!parts.ids_didfit) continue;

    if(parts.fido_qiv < 5000) continue;
    if(parts.fido_qid/8300 > 700) continue;
    if(parts.ids_chi2/(parts.fido_nidtubes+parts.fido_nivtubes-6)>10)
      continue;

    // This correctly uses the *uncorrected* position
    if(pow(parts.ids_end_x, 2)+pow(parts.ids_end_y, 2) > pow(1708-35,2))
      continue;

    // This correctly uses the *uncorrected* position
    if(parts.ids_end_z < -1786+35) continue;

    const double dedxslant = parts.ids_entr_z
      + 62*parts.fido_qiv/(parts.id_ivlen-parts.id_buflen);

    if(dedxslant > 11500) continue;

    const double chi2qual = parts.ids_chi2-parts.id_chi2;
    if(chi2qual > 800) continue;

    if(parts.fido_nidtubes+parts.fido_nivtubes < 6) continue;

    if(pow(parts.id_entr_x, 2)+pow(parts.id_entr_y, 2) < pow(1000, 2) &&
       pow(parts.ids_entr_x,2)+pow(parts.ids_entr_y,2) > pow(2758, 2))
      continue;


    const double mux = fidocorrxy(parts.ids_end_x),
                 muy = fidocorrxy(parts.ids_end_y),
                 muz = fidocorrz( parts.ids_end_z);
    const double imux = parts.ids_entr_x,
                 imuy = parts.ids_entr_y,
                 imuz = parts.ids_entr_z;

    const double li9tomu = sqrt(pow(mux-li9x, 2)
                               +pow(muy-li9y, 2)
                               +pow(muz-li9z, 2));

    // And they must be near each other.
    if(li9tomu > 800) continue;

    // Ok, we've almost accepted this muon.  Now we have to do the hard
    // work of figuring out if it happened right after another muon,
    // because if so the neutron count we're about to do is unreliable.


    // First save the structure we're about to tromp all over
    dataparts selectedmuon;
    memcpy(&selectedmuon, &parts, sizeof(dataparts));

    // Require at least 500us since the last muon so we don't count
    // neutrons from the previous one as belonging to this one.
    bool previousmuon = false;
    for(int moreback = 1; ; moreback++){
      if(prompt-back-moreback < 1) break; // BOR and didn't find any
      ctree->GetEvent(prompt-back-moreback);

      // Out of 500us window and didn't find any
      if(selectedmuon.trgtime - parts.trgtime > 500e3) break;

      if(parts.coinov || parts.fido_qiv > 5000 || parts.ctEvisID > 60)
        previousmuon = true;
    }

    if(previousmuon) continue;

    // restore our muon after the look back for other muons
    memcpy(&parts, &selectedmuon, sizeof(dataparts));

    // Open up a big window
    if(prompttime - parts.trgtime > window) return;

    const double mutime = parts.trgtime;
    const double fidoqid = parts.fido_qid;

    ctree->GetEvent(prompt-back+1);

    const float miche = parts.deltaT < 5500? parts.ctEvisID:0.;
    const float micht = parts.deltaT < 5500? parts.deltaT:0.;

    // XXX really need to do more careful counting like in b12search
    unsigned int nneutron = 0, nneutronnotmichel = 0;
    for(int forward = 1; back-forward > 0; forward++){
      ctree->GetEvent(prompt-back+forward);

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
      if(sqrt(pow(mux - bamacorrxy(parts.ctX[0], parts.ctEvisID), 2)
             +pow(muy - bamacorrxy(parts.ctX[1], parts.ctEvisID), 2)
             +pow(muz - bamacorrz( parts.ctX[2], parts.ctEvisID), 2))
         > 1000) continue;

      nneutron++;
      if(parts.trgtime - mutime > 5500) nneutronnotmichel++;
    }

    printf("%d %d %d %lf %f %u %u %f %f %f %f %f %f %f %f %f "
           "%f %f %f %f %f %d\n",
           parts.run, prompt, prompt-back,
           (prompttime - mutime)/1e6, li9tomu, nneutron, nneutronnotmichel,
           li9x, li9y, li9z, mux, muy, muz, imux, imuy, imuz,
           miche, micht, fidoqid, chi2qual, dedxslant,
           earlymich(parts.run, prompt-back));
    fflush(stdout);
  }
}

int main()
{
  gErrorIgnoreLevel = kError;
  unsigned int errcode = 0;
  std::string line;

  printf("run:trig:mutrig:dt:dist:n:nlate:dx:dy:dz:mx:my:mz:"
         "imx:imy:imz:miche:micht:fidoqid:chi2qual:dedxslant:"
         "earlymich\n");

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
      errcode |= 0x01;
      goto cleanup;
    }

    if(!finfile || finfile->IsZombie()){
      fprintf(stderr, "I couldn't read fido run %d\n", run);
      errcode |= 0x02;
      goto cleanup;
    }

    ctree = (TTree *)cinfile->Get("data");
    if(!ctree){
      fprintf(stderr, "Run %d lacks a data tree!\n", run);
      errcode |= 0x04;
      goto cleanup;
    }

    ftree = (TTree *)finfile->Get("RecoMuonFIDOInfoTree");
    if(!ftree){
      fprintf(stderr, "Run %d lacks a RecoMuonFIDOInfoTree!\n", run);
      errcode |= 0x08;
      goto cleanup;
    }

    if(ctree->GetEntries() != ftree->GetEntries()){
      fprintf(stderr,
              "Run %d: cheetah has %ld entries, but fido has %ld\n",
              run,
              long(ctree->GetEntries()), long(ftree->GetEntries()));
      errcode |= 0x10;
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
    fSBA(id_entr_x);
    fSBA(id_entr_y);
    fSBA(ids_entr_x);
    fSBA(ids_entr_y);
    fSBA(ids_entr_z);
    fSBA(ids_chi2);
    fSBA(ids_ivlen);
    fSBA(ids_buflen);
    fSBA(id_chi2);
    fSBA(id_ivlen);
    fSBA(id_buflen);

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

  if(errcode != 0) fprintf(stderr, "Errors encountered, see above\n");

  return errcode;
}
