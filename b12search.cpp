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
  bool fido_stop;
  float fido_qiv, fido_qid;
  float fido_chi2;
  float fido_endx, fido_endy, fido_endz;
  int fido_nidtubes, fido_nivtubes;
  double deltaT;
  double trgtime;
  float ctX[3];
  float ctEvisID;
  float qrms, qdiff;
};

static const unsigned int noff = 90;
static const char * turnoff[noff] = { 
  "nev", "trgWord", "date", "tref", "trefIV", "trefextIV",
  "nhitIV", "nhit", "ctq", "ctnpe", "ctqIV", "ctnpeIV", "ctnpulse",
  "ctnbadch", "ctrho", "ctphi", "ctR", "ctmqtq", "ctmqtqflag",
  "ctgoodness", "ctnbadchIV", "ctnpulseIV",
  "vctnpulseIV", "vctqIV", "vcttimeIV", "vctnpulse", "vctq",
  "vcttime", "pmtmultpe", "pmtmultpe_IV", "IVX", 
  "cttrise", "ctfwhm", "ctt2tot", "cttpeak", "cttmean",
  "ctIDMuDeltaT", "ctIVMuDeltaT", "HEMuDeltaT", "ctlightflux",
  "ctFlagMu", "ctXmuInGC", "ctXmuInIV", "ctXmuOuIV", "ctqtot",
  "ctqtotIV", "ctsphericity", "ctaplanarity", "lilike", 
  "timeid", "timeiv", "ovtrigid", "ttovtrig",
  "novhit", "novupxy", "novloxy", "novtrk", "ovupxyx",
  "ovupxyy", "ovupxyz", "ovtightupxy", "ovupxylike", "ovloxyx",
  "ovloxyy", "ovloxyz", "ovtightloxy", "ovloxylike", "ovtrkx",
  "ovtrky", "ovtrkth", "ovtrkphi", "ovtrklike", "ovbadtrk",
  "ovtighttrk", "hamx", "hamxe", "hamth", "hamphi", "fido_didfit", 
  "fido_used_ov", "fido_minuit_happiness", "fido_ivlen",
  "fido_buflen", "fido_gclen", "fido_targlen", "fido_entrx",
  "fido_entry", "fido_entrz", "fido_th", "fido_phi"
};

static bool lightnoise(const float qrms, const float mqtq,
                       const float rmsts, const float qdiff)
{
  if(mqtq > 0.12 || mqtq < 0) return true;
  if(qdiff > 30e3) return true;
  if(rmsts >= 36 && qrms >= 464 - 8*rmsts) return true;
  return false;
}


static void searchfrommuon(dataparts & parts, TTree * const data,
                           const unsigned int muoni)
{
  const double mutime = parts.trgtime;
  const double mux = parts.fido_endx,
               muy = parts.fido_endy,
               muz = parts.fido_endz;
  const int mutrgid = parts.trgId;
  const int murun = parts.run;
  const bool mucoinov = parts.coinov;
  const float murchi2 =
    parts.fido_chi2/(parts.fido_nidtubes+parts.fido_nivtubes-6);
  const float mufidoqiv = parts.fido_qiv;

  unsigned int got = 0;
  for(unsigned int i = muoni+1; i < data->GetEntries(); i++){
    data->GetEntry(i);  

    if(parts.run != murun) break; // Stop at run boundaries

    const double b12time = parts.trgtime; 
    const double dt_ms = (b12time - mutime)/1e6;
    
    if(dt_ms < 1) continue; // Go past all H neutron captures

    if(dt_ms > 1000) break; // Stop looking after lots of b12 lifetimes

    // Ignore low energy accidentals and H-neutrons
    if(parts.ctEvisID < 3) continue;

    // No IV energy
    if(parts.fido_qiv > 1000) continue;
    if(parts.coinov) continue;

    // pass light noise
    if(lightnoise(parts.qrms, parts.ctmqtqall,
                  parts.ctrmsts, parts.qdiff)) continue;

    // Go up high enough to also get B-8, Li-8 end points of 17 and 16.
    if(parts.ctEvisID > 20) continue;

    const double b12x = parts.ctX[0]*(0.970+0.013*parts.ctEvisID),
                 b12y = parts.ctX[1]*(0.970+0.013*parts.ctEvisID),
                 b12z = parts.ctX[2]*(0.985+0.005*parts.ctEvisID);

    const double b12tomu = sqrt(pow(mux-b12x, 2)
                               +pow(muy-b12y, 2)
                               +pow(muz-b12z, 2));

    // And they must be near each other.
    if(b12tomu > 800) continue;

    printf("B-12_for %d %d %d is %d dt %lf distance %f energy %f b12_at %f %f %f chi2,qiv: %f %f ",
        parts.run, mutrgid, mucoinov, parts.trgId, dt_ms, b12tomu, parts.ctEvisID, b12x, b12y, b12z,
        murchi2, mufidoqiv);
    got++;

    if(got >= 2) break;
  }
  if(got){
    printf("\n");
    fflush(stdout);
  }
}

static void search(dataparts & parts, TTree * const data)
{
  for(unsigned int mi = 0; mi < data->GetEntries()-1; mi++){
    data->GetBranch("fido_stop")->GetEntry(mi);

    // Must be reasonably sure that it's a stopper
    if(!parts.fido_stop) continue;

    data->GetEntry(mi);

    if(parts.fido_qiv < 5000) continue;
    if(parts.fido_qid/8300 > 700) continue;
    if(parts.fido_chi2/(parts.fido_nidtubes+parts.fido_nivtubes-6) > 10)
      continue;

    if(parts.fido_nidtubes+parts.fido_nivtubes < 30) continue;

    data->GetEntry(mi+1);

    // And must not have a Michel, loose cuts
    if(parts.deltaT < 5500 && parts.ctEvisID > 12) continue;

    data->GetEntry(mi);
    searchfrommuon(parts, data, mi);
  }
}

int main(int argc, char ** argv)
{
  gErrorIgnoreLevel = kFatal;
  int errcode = 0;
  for(int i = 1; i < argc; i++){
    TFile * const infile = new TFile(argv[i], "read");

    if(!infile || infile->IsZombie()){
      fprintf(stderr, "I couldn't read %s\n", argv[i]);
      errcode |= 1;
      continue;
    }

    TTree * const data = (TTree *)infile->Get("data");
    if(!data){
      fprintf(stderr, "%s doesn't have a data tree!\n", argv[i]);
      errcode |= 2;
      continue;
    }

    for(unsigned int i = 0; i < noff; i++)
      data->SetBranchStatus(turnoff[i], 0);

    dataparts parts;
    memset(&parts, 0, sizeof(parts));
    #define SBA(x) data->SetBranchAddress(#x, &parts.x);
    SBA(fido_stop);
    SBA(fido_endx);
    SBA(fido_endy);
    SBA(fido_endz);
    SBA(fido_qiv);
    SBA(fido_qid);
    SBA(fido_chi2);
    SBA(fido_nidtubes);
    SBA(fido_nivtubes);
    SBA(deltaT);
    SBA(trgtime);
    SBA(run);
    SBA(trgId);
    SBA(ctmqtqall);
    SBA(ctrmsts);
    SBA(coinov);
    SBA(ctEvisID);
    SBA(qrms);
    SBA(qdiff);
    if(data->SetBranchAddress("ctX", parts.ctX) < 0)
      data->SetBranchAddress("ctX[3]", parts.ctX);

    search(parts, data);
  
    delete data;
    delete infile;
  }
  return errcode;
}
