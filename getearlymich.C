#include "TTree.h"
#include "TFile.h"
#include <iostream>

void doit(const int run, const int mutrig, TTree * ts, TTree * tl)
{
  if(ts->GetEntries(Form("run == %d && mutrig == %d", run, mutrig)))
    cout<<ts->GetEntries(Form("run==%d&&mutrig==%d&&earlymich",run,mutrig))<<endl;
  else
    cout<<tl->GetEntries(Form("run==%d&&mutrig==%d&&earlymich",run,mutrig))<<endl;
}

void getearlymich()
{
  TFile fs("/cp/s4/strait/fullfido-100s-3-25MeV-20140925-earlymich.root", "read");
  TTree * ts = (TTree *)fs.Get("t"); 

  TFile fl("/cp/s4/strait/fullfido-100s-3-25MeV-20140925-earlymich.root", "read");
  TTree * tl = (TTree *)fs.Get("t"); 
  
  tl->SetBranchStatus("*", 0);
  tl->SetBranchStatus("run", 1);
  tl->SetBranchStatus("mutrig", 1);
  tl->SetBranchStatus("earlymich", 1);
  ts->SetBranchStatus("*", 0);
  ts->SetBranchStatus("run", 1);
  ts->SetBranchStatus("mutrig", 1);
  ts->SetBranchStatus("earlymich", 1);

  int run, mutrig;
  while(cin >> run >> mutrig) doit(run, mutrig, ts, tl);
}
