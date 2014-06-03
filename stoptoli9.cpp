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
  int run;
  bool fido_stop;
  float fido_qiv, fido_qid;
  float fido_chi2;
  float fido_endx, fido_endy, fido_endz;
  int fido_nidtubes, fido_nivtubes;
  double deltaT;
  double trgtime;
  float ctX[3];
};

void stopper_search(dataparts & parts, TTree * const data, 
                    const int prompt)
{
  if(prompt >= data->GetEntries()){
    fprintf(stderr,
            "Prompt event %d given, but file has only %ld events\n",
            prompt, (long int)data->GetEntries());
    return;
  }

  data->GetEntry(prompt);
  const double prompttime = parts.trgtime; 
  const double li9x=parts.ctX[0],
               li9y=parts.ctX[1],
               li9z=parts.ctX[2];

  for(int back = 1; ; back++){
    if(prompt-back < 0) return;
    data->GetEvent(prompt-back);

    // Must be reasonably sure that it's a stopper
    if(!parts.fido_stop) continue;
    if(parts.fido_qiv < 5000) continue;
    if(parts.fido_qid/8300 > 700) continue;
    if(parts.fido_chi2/(parts.fido_nidtubes+parts.fido_nivtubes-6) > 10)
      continue;

    // and not too long before the Li-9 candidate
    if(prompttime - parts.trgtime > 5e9) return;

    const double mux = parts.fido_endx, muy = parts.fido_endy,
                 muz = parts.fido_endz;

    const double li9tomu = sqrt(pow(mux-li9x, 2)
                               +pow(muy-li9y, 2)
                               +pow(muz-li9z, 2));

    // And they must be near each other.
    if(li9tomu > 1000) continue;

    data->GetEvent(prompt-back+1);

    // And must not have a Michel, with loose cuts.
    if(parts.deltaT < 5500 && parts.fido_qid/8300 > 12) continue;

    printf("Muon for %d %d is %d dt %lf ms distance %f at %f %f %f\n",
           parts.run, prompt, prompt-back,
           (prompttime - parts.trgtime)/1e6, li9tomu, li9x, li9y, li9z);

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

    const unsigned int noff = 97;
    const char * turnoff[noff] = { 
      "nev", "trgId", "trgWord", "date", "tref", "trefIV", "trefextIV",
      "nhitIV", "nhit", "ctq", "ctnpe", "ctqIV", "ctnpeIV", "ctnpulse",
      "ctnbadch", "ctrho", "ctphi", "ctR", "ctmqtq", "ctmqtqflag",
      "ctmqtqall", "ctgoodness", "ctnbadchIV", "ctnpulseIV",
      "vctnpulseIV", "vctqIV", "vcttimeIV", "vctnpulse", "vctq",
      "vcttime", "pmtmultpe", "pmtmultpe_IV", "IVX", "ctrmsts",
      "cttrise", "ctfwhm", "ctt2tot", "cttpeak", "cttmean",
      "ctIDMuDeltaT", "ctIVMuDeltaT", "HEMuDeltaT", "ctlightflux",
      "ctFlagMu", "ctXmuInGC", "ctXmuInIV", "ctXmuOuIV", "ctqtot",
      "ctqtotIV", "ctsphericity", "ctaplanarity", "lilike", "qrms", 
      "qdiff", "timeid", "timeiv", "ctEvisID", "ovtrigid", "ttovtrig",
      "coinov", "novhit", "novupxy", "novloxy", "novtrk", "ovupxyx",
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
    data->SetBranchAddress("ctX", parts.ctX);

    int prompt;
    while(ss >> prompt) stopper_search(parts, data, prompt);
  
    delete data;
    delete infile;
  }

}
