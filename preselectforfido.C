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
  double TrigTime, EvisIDg[3];
  unsigned int TriggerID;
  GI->SetBranchAddress("TrigTime", &TrigTime);
  GI->SetBranchAddress("TriggerID", &TriggerID);
  GI->SetBranchAddress("EvisIDg", EvisIDg);

  for(int i = 0; i < GI->GetEntries(); i++){
    GI->GetEntry(i);
    
    if(EvisID[0] > 70){
      const unsigned int muontrigid = TriggerID;
      const double muontrigtime = TrigTime;
      bool ok = true;

      // Require no Michel
      // Contrariwise, *do* allow Michels, because
      // we need decaying muons for finding the
      // accidental neutron background. Now, do we
      // really need them *reconstructed*? Yes, because
      // the accidental rate is a function of
      // position (even if the neutron efficiency
      // is not).
      //
      // Do we need to reconstruct *all* decaying
      // muons to get adequate statistics? Maybe not,
      // so come back here and reenable this test 
      // modified so it lets through a tenth or
      // whatever if it turns out we are wasting time.
      //
      // But since there may also be gammas over 10 MeV
      // or interesting features of the Michel spectrum
      // I bet we want to stick with doing them all.
      /*for(int j = i+1; j < GI->GetEntries(); j++){
        GI->GetEntry(j);
        if(TrigTime - muontrigtime <= 512) continue;
        if(TrigTime - muontrigtime > 10*2197.) break;
        if(EvisIDg[0] > 10){
          ok = false;
          break;
        }
      }*/
      if(ok) printf("%d\n", muontrigid);
    }
  }
}
