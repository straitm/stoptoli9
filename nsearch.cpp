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
  float ctEvisID;
  float qdiff, qrms;
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

static const unsigned int noff = 91;
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


static bool searchfrommuon(dataparts & parts, TTree * const data,
                           const unsigned int muoni)
{
  const double mutime = parts.trgtime;
  const double mux = parts.fido_endx,
               muy = parts.fido_endy,
               muz = parts.fido_endz;

  unsigned int nneutron = 0;
  for(unsigned int i = muoni+1; i < data->GetEntries(); i++){
    data->GetEntry(i);  

    const double ntime = parts.trgtime; 

    const double dt = ntime - mutime;

    // Exclude the first 10 us to get past low energy
    // Michels
    if(dt < 10e3) continue;
    
    // Stop looking after ~4 nH lifetimes
    if(dt > 800e3) return false;

    // Stop counting if there's another muon
    if(parts.fido_qiv > 1000 || parts.ctEvisID > 50 || parts.coinov == 1)
      break;

    // right energy for a neutron capture
    if(!((parts.ctEvisID > 1.8 && parts.ctEvisID < 2.6) ||
         (parts.ctEvisID > 4.0 && parts.ctEvisID < 10 )))
      continue;
      
    const double
      nx = parts.ctX[0], ny = parts.ctX[1], nz = parts.ctX[2];

    const double ntomu =
      sqrt(pow(mux-nx, 2) + pow(muy-ny, 2) + pow(muz-nz, 2));

    // And they must be near each other.
    if(ntomu > 1000) continue;

    // pass light noise
    if(lightnoise(parts.qrms, parts.ctmqtqall,
                  parts.ctrmsts, parts.qdiff)) continue;

    printf("n %d for %d %d is %d dt %lf us distance %f energy %f "
           "muon_at %f %f %f\n", ++nneutron,
           parts.run, muoni, i, dt/1000, ntomu, parts.ctEvisID,
           mux, muy, muz);
    fflush(stdout);
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
    SBA(qdiff);
    SBA(qrms);
    SBA(ctEvisID);
    data->SetBranchAddress("ctX", parts.ctX);

    search(parts, data);
  
    delete data;
    delete infile;
  }
  return errcode;
}
