#include "TTree.h"
#include "TFile.h"
#include "TBranch.h"
#include "TObject.h"
#include <iostream>
using namespace std;

void earlymich()
{
  TFile *_file1 = TFile::Open("/cp/s4/strait/feb25.t.withmichonly.root");
  TTree * pulse = (TTree *)_file1->Get("t");
  pulse->SetName("pulse");

  TFile *_file0 = TFile::Open("/cp/s4/strait/fullfido-100s-0-25MeV-20140925-earlymich.root", "update");
  TTree * b12 = (TTree *)_file0->Get("t");
  b12->SetName("b12");
  cout << b12->GetEntries() << " entries" << endl;

  int run, trig;
  b12->SetBranchAddress("run", &run);
  b12->SetBranchAddress("mutrig", &trig);
  TBranch * runbr   = b12->GetBranch("run");
  TBranch * trigbr  = b12->GetBranch("mutrig");

  float prun, ptrig;
  pulse->SetBranchAddress("run", &prun);
  pulse->SetBranchAddress("trig", &ptrig);
  TBranch * prunbr = pulse->GetBranch("run");
  TBranch * ptrigbr = pulse->GetBranch("trig");

  int lastpulse = 0;
  bool earlymich;
  TBranch * earlymichbr = b12->Branch("earlymich", &earlymich, "earlymich/O");

  const int ent = b12->GetEntries();
  for(int i = 0; i < ent; i++){
    if(i%100000==0) cout << 100*float(i)/ent << "%" << endl;

    runbr->GetEntry(i);
    trigbr->GetEntry(i);
    earlymich = false;
    for(int j = lastpulse; j < pulse->GetEntries(); j++){
      prunbr->GetEntry(j);
      if(prun > run) break;
      ptrigbr->GetEntry(j);
      if(prun == run && ptrig > trig) break;

      if(prun == run && ptrig == trig){
        lastpulse = j; // start exactly where we stopped because
                       // we may be about to answer for the same
                       // muon again.
        earlymich = true;
        goto end;
      }
    }
    end: 
    earlymichbr->Fill(); 
  }
  b12->SetName("t");
  b12->Write("", TObject::kOverwrite);
}
