#include "TFile.h"
#include "TTree.h"
#include <stdio.h>

void preselectforfido(const char * const jpfilename)
{
  TFile * jpfile = new TFile(jpfilename, "read");
  if(!jpfile || jpfile->IsZombie()){
    fprintf(stderr, "Could not open %s\n", jpfilename);
    exit(1);
  }
  TTree * GI = dynamic_cast<TTree *>(jpfile->Get("GI"));
  if(!GI){
    fprintf(stderr, "Could not get GI tree from %s\n", jpfilename);
    exit(1);
  }
  double TrigTime, EvisID;
  unsigned int TriggerID;
  GI->SetBranchAddress("TrigTime", &TrigTime);
  GI->SetBranchAddress("TriggerID", &TriggerID);
  GI->SetBranchAddress("EvisID", &EvisID);

  for(int i = 0; i < GI->GetEntries(); i++){
    GI->GetEntry(i);
    
    if(EvisID > 70){
      const unsigned int muontrigid = TriggerID;
      const double muontrigtime = TrigTime;
      bool ok = true;

      // Require no Michel
      for(int j = i+1; j < GI->GetEntries(); j++){
        GI->GetEntry(j);
        if(TrigTime - muontrigtime <= 512) continue;
        if(TrigTime - muontrigtime > 10*2197.) break;
        if(EvisID > 10){
          ok = false;
          break;
        }
      }
      if(ok) printf("%d\n", muontrigid);
    }
  }
}
