#include "TTree.h"
#include "TFile.h"
#include <iostream>

void doit(const int run, const int mutrig, TTree * ts)
{
  cout<<ts->GetEntries(Form("run==%d&&mutrig==%d&&earlymich",run,mutrig))<<endl;
}

void getearlymich()
{
  TFile fs("/cp/s4/strait/fullfido-100s-3-25MeV-20150204-earlymich.root", "read");
  TTree * ts = (TTree *)fs.Get("t"); 

  ts->SetBranchStatus("*", 0);
  ts->SetBranchStatus("run", 1);
  ts->SetBranchStatus("mutrig", 1);
  ts->SetBranchStatus("earlymich", 1);

  int run, mutrig;
  while(cin >> run >> mutrig) doit(run, mutrig, ts);
}
