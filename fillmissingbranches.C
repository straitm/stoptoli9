/* Creates a few branches and fills them with zeros for compatibility
 * with old ntuples */

#include "TTree.h"
#include "TFile.h"
#include "TBranch.h"
#include "TObject.h"
#include <iostream>
using namespace std;

void fillmissingbranches()
{
  TFile *_file0 = TFile::Open("/cp/s4/strait/fullfido-300s-3-25MeV-20150924-post3rdpub.root", "update");
  TTree * b12 = (TTree *)_file0->Get("t");
  b12->SetName("b12");
  cout << b12->GetEntries() << " entries" << endl;

  bool earlymich = false;
  float deadt = 0, deade = 0;
  TBranch * earlymichbr = b12->Branch("earlymich", &earlymich, "earlymich/O");
  TBranch * deadtbr = b12->Branch("deadt", &deadt, "deadt/F");
  TBranch * deadebr = b12->Branch("deade", &deade, "deade/F");

  const int ent = b12->GetEntries();
  for(int i = 0; i < ent; i++){
    earlymichbr->Fill(); 
    deadtbr->Fill(); 
    deadebr->Fill(); 
  }
  b12->SetName("t");
  b12->Write("", TObject::kOverwrite);
}
