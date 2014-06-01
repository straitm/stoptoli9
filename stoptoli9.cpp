#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "TTree.h"
#include "TFile.h"
#include "TError.h"

struct dataparts{
  int run;
  bool fido_stop;
  float fido_qiv, fido_qid;
  float fido_chi2;
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

    const double mux=parts.ctX[0], muy=parts.ctX[1], muz=parts.ctX[2];

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

int main(int argc, char ** argv)
{
  if(argc < 3){
    fprintf(stderr,"Syntax: reduced_file prompt_ev [more prompts]\n");
    exit(1);
  }

  TFile * const infile = new TFile(argv[1], "read");
  if(!infile || infile->IsZombie()){
    fprintf(stderr, "I couldn't read %s\n", argv[1]);
    exit(1);
  }

  TTree * const data = (TTree *)infile->Get("data");
  if(!data){
    fprintf(stderr, "%s doesn't have a data tree!\n", argv[1]);
    exit(1);
  }

  dataparts parts;
  #define SBA(x) data->SetBranchAddress(#x, &parts.x);
  SBA(fido_stop);
  SBA(fido_qiv);
  SBA(fido_qid);
  SBA(fido_chi2);
  SBA(fido_nidtubes);
  SBA(fido_nivtubes);
  SBA(deltaT);
  SBA(trgtime);
  SBA(run);
  data->SetBranchAddress("ctX", parts.ctX);

  for(int i = 2; i < argc; i++){
    stopper_search(parts, data, atoi(argv[i]));
  }
  
  string line;
  while(std::getline(line)){
    printf("%s\n", line.c_str());
  }

}
