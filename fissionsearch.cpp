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
  float ctEvisID;
  int run;
  bool fido_stop;
  float fido_qiv, fido_qid;
  float fido_chi2;
  float fido_endx, fido_endy, fido_endz;
  float qdiff, qrms;
  int fido_nidtubes, fido_nivtubes;
  double deltaT;
  double trgtime;
  float ctmqtqall, ctrmsts;
  float ctX[3];
};

static bool lightnoise(const float qrms, const float mqtq,
                       const float rmsts, const float qdiff)
{
  if(mqtq > 0.12 || mqtq < 0) return true;
  if(qdiff > 30e3) return true;
  if(rmsts >= 36 && qrms >= 464 - 8*rmsts) return true;
  return false;
}

static void fission_search(dataparts & parts, TTree * const data,
                           const int prompt, const int delayed)
{
  if(prompt >= data->GetEntries() ||
     delayed >= data->GetEntries()){
    fprintf(stderr,
            "Events %d+%d given, but file has only %ld events\n",
            prompt, delayed, (long int)data->GetEntries());
    return;
  }

  data->GetEntry(prompt);
  const double prompttime = parts.trgtime; 
  const double promptE = parts.ctEvisID;
  const double px=parts.ctX[0],
               py=parts.ctX[1],
               pz=parts.ctX[2];

  for(int e = delayed+1; e < data->GetEntries(); e++){
    data->GetEvent(e);

    // Open up a big window
    if(parts.trgtime - prompttime > 1e9) return;

    const double candx = parts.ctX[0], candy = parts.ctX[1],
                 candz = parts.ctX[2];

    const double dist = sqrt(pow(candx-px, 2)
                            +pow(candy-py, 2)
                            +pow(candz-pz, 2));

    if(dist > 500) continue;

    if(parts.fido_qiv > 1000 || parts.ctEvisID > 50 || parts.coinov)
      continue;

    if(parts.ctEvisID < 0.5 || parts.ctEvisID > 6) continue;

    if(lightnoise(parts.qrms, parts.ctmqtqall,
                  parts.ctrmsts, parts.qdiff)) continue;

    printf("Late event for %d %d is %d dt %lf ms distance %f with "
           "E %f promt %f at %f %f %f\n", parts.run, prompt, e,
           (parts.trgtime - prompttime)/1e6, dist, parts.ctEvisID, 
           promptE, px, py, pz);

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

    TFile * const infile = new TFile(Form("/cp/s4/dchooz/cheetah/"
      "prod-08-05_p01_v2/reduced.Run00%05d_Seq010.root", run), "read");

    if(!infile || infile->IsZombie()){
      fprintf(stderr, "I couldn't read run %d\n", run);
      exit(1);
    }

    TTree * const data = (TTree *)infile->Get("data");
    if(!data){
      fprintf(stderr, "Run %d doesn't have a data tree!\n", run);
      exit(1);
    }

    const unsigned int noff = 91;
    const char * turnoff[noff] = { 
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

    for(unsigned int i = 0; i < noff; i++)
      data->SetBranchStatus(turnoff[i], 0);

    dataparts parts;
    #define SBA(x) data->SetBranchAddress(#x, &parts.x);
    SBA(qdiff);
    SBA(ctEvisID);
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
    SBA(ctmqtqall);
    SBA(ctrmsts);
    SBA(qrms);
    SBA(run);
    SBA(coinov);
    data->SetBranchAddress("ctX", parts.ctX);

    int prompt, delayed;
    while(ss >> prompt >> delayed)
      fission_search(parts, data, prompt, delayed);
  
    delete data;
    delete infile;
  }

}
