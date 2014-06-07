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
  int run;
  float ctmqtqall, ctrmsts;
  bool fido_stop;
  float fido_qiv, fido_qid;
  float fido_chi2;
  float fido_endx, fido_endy, fido_endz;
  int fido_nidtubes, fido_nivtubes;
  double deltaT;
  double trgtime;
  float ctX[3];
};

static const unsigned int noff = 94;
static const char * turnoff[noff] = { 
  "nev", "trgId", "trgWord", "date", "tref", "trefIV", "trefextIV",
  "nhitIV", "nhit", "ctq", "ctnpe", "ctqIV", "ctnpeIV", "ctnpulse",
  "ctnbadch", "ctrho", "ctphi", "ctR", "ctmqtq", "ctmqtqflag",
  "ctgoodness", "ctnbadchIV", "ctnpulseIV",
  "vctnpulseIV", "vctqIV", "vcttimeIV", "vctnpulse", "vctq",
  "vcttime", "pmtmultpe", "pmtmultpe_IV", "IVX", 
  "cttrise", "ctfwhm", "ctt2tot", "cttpeak", "cttmean",
  "ctIDMuDeltaT", "ctIVMuDeltaT", "HEMuDeltaT", "ctlightflux",
  "ctFlagMu", "ctXmuInGC", "ctXmuInIV", "ctXmuOuIV", "ctqtot",
  "ctqtotIV", "ctsphericity", "ctaplanarity", "lilike", "qrms", 
  "qdiff", "timeid", "timeiv", "ctEvisID", "ovtrigid", "ttovtrig",
  "novhit", "novupxy", "novloxy", "novtrk", "ovupxyx",
  "ovupxyy", "ovupxyz", "ovtightupxy", "ovupxylike", "ovloxyx",
  "ovloxyy", "ovloxyz", "ovtightloxy", "ovloxylike", "ovtrkx",
  "ovtrky", "ovtrkth", "ovtrkphi", "ovtrklike", "ovbadtrk",
  "ovtighttrk", "hamx", "hamxe", "hamth", "hamphi", "fido_didfit", 
  "fido_used_ov", "fido_minuit_happiness", "fido_ivlen",
  "fido_buflen", "fido_gclen", "fido_targlen", "fido_entrx",
  "fido_entry", "fido_entrz", "fido_th", "fido_phi"
};


static bool searchfrommuon(dataparts & parts, TTree * const data,
                           const unsigned int muoni)
{
  const double mutime = parts.trgtime;
  const double mux = parts.fido_endx,
               muy = parts.fido_endy,
               muz = parts.fido_endz;

  for(unsigned int i = muoni+1; i < data->GetEntries(); i++){
    data->GetEntry(i);  

    const double b12time = parts.trgtime; 

    const double dt_ms = (b12time - mutime)/1e6;
    
    // Go past all H neutron captures
    if(dt_ms < 1) continue;

    // Stop looking after several b12 lifetimes
    if(dt_ms > 100) return false;

    const double b12x = parts.ctX[0],
                 b12y = parts.ctX[1],
                 b12z = parts.ctX[2];

    const double b12tomu = sqrt(pow(mux-b12x, 2)
                               +pow(muy-b12y, 2)
                               +pow(muz-b12z, 2));

    // And they must be near each other.
    if(b12tomu > 500) continue;

    // Ignore low energy accidentals and H-neutrons
    if(parts.fido_qid/8300 < 3) continue;

    // No IV energy
    if(parts.fido_qiv > 1000) continue;

    if(parts.coinov) continue;

    if(parts.ctmqtqall > 0.12) continue;
    if(parts.ctrmsts > 40) continue;

    // Must be below b12 endpoint 
    if(parts.fido_qid/8300 > 14) continue;

    printf("B-12 for %d %d is %d dt %lf ms distance %f energy %f\n",
        parts.run, muoni, i, dt_ms, b12tomu, parts.fido_qid/8300);

    return true;
  }
  return false;
}

static void search(dataparts & parts, TTree * const data)
{
  for(unsigned int mi = 0; mi < data->GetEntries()-1; mi++){
    data->GetEntry(mi);

    // Must be reasonably sure that it's a stopper
    if(!parts.fido_stop) continue;
    if(parts.fido_qiv < 5000) continue;
    if(parts.fido_qid/8300 > 700) continue;
    if(parts.fido_chi2/(parts.fido_nidtubes+parts.fido_nivtubes-6) > 10)
      continue;

    // No muons near the edges of the GC, both to remove badly
    // reconstructed muons and to avoid accidentals near the edges.
    if(abs(parts.fido_endz) > 1300 ||
       sqrt(pow(parts.fido_endx,2)+pow(parts.fido_endy,2)) > 1300)
      continue;

    data->GetEntry(mi+1);

    // And must not have a Michel, moderate cuts
    if(parts.deltaT < 5500 && parts.fido_qid/8300 > 12 &&
                              parts.fido_qid/8300 < 70) continue;

    data->GetEntry(mi);
    searchfrommuon(parts, data, mi);
  }
}

int main(int argc, char ** argv)
{
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
    SBA(ctmqtqall);
    SBA(ctrmsts);
    SBA(coinov);
    data->SetBranchAddress("ctX", parts.ctX);

    search(parts, data);
  
    delete data;
    delete infile;
  }
  return errcode;
}
