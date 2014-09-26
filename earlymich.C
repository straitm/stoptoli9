#include "TCanvas.h"
#include "TH1.h"
#include "TTree.h"
#include "TFile.h"
#include "TBranch.h"
#include <iostream>
using namespace std;

void earlymich()
{
  TFile *_file0 = TFile::Open("/cp/s4/strait/fullfido-100s-3-25MeV-20140925.root");
  TTree * b12 = (TTree *)_file0->Get("t");
  cout << b12->GetEntries() << " entries" << endl;
  b12->SetName("b12");
  TFile *_file1 = TFile::Open("/cp/s4/strait/feb25.t.withmichonly.root");
  TTree * pulse = (TTree *)_file1->Get("t");
  pulse->SetName("pulse");
  float dt;
  float miche;
  b12->SetBranchAddress("dt", &dt);
  b12->SetBranchAddress("miche", &miche);
  int run, trig;
  b12->SetBranchAddress("run", &run);
  b12->SetBranchAddress("mutrig", &trig);
  TCanvas * c1 = new TCanvas();
  c1->SetLogy();
  TH1D * all = new TH1D("all", "", 200, 1, 2010);
  TH1D * vetoed = new TH1D("vetoed", "", 200, 1, 2010);
  vetoed->SetLineColor(kRed);
  vetoed->SetMarkerColor(kRed);
  all->Draw();
  vetoed->Draw("same");

  TBranch * dtbr    = b12->GetBranch("dt");
  TBranch * michebr = b12->GetBranch("miche");
  TBranch * runbr   = b12->GetBranch("run");
  TBranch * trigbr  = b12->GetBranch("mutrig");

  const int ent = b12->GetEntries();

  pulse->SetBranchStatus("*", 0);
  pulse->SetBranchStatus("run", 1);
  pulse->SetBranchStatus("trig", 1);

  float prun, ptrig;
  pulse->SetBranchAddress("run", &prun);
  pulse->SetBranchAddress("trig", &ptrig);
  int lastpulse = 0;

  TBranch * prunbr = pulse->GetBranch("run");
  TBranch * ptrigbr = pulse->GetBranch("trig");

  for(int i = 0; i < ent; i++){
    if((i&0xfffff)==0){
      cout << i << "\t" << vetoed->GetEntries() << "/" << all->GetEntries() << endl;
      c1->Update();
      c1->Modified();
    }
    dtbr->GetEntry(i);
    if(dt == 0) continue;
    michebr->GetEntry(i);
    if(miche > 12) continue;

    all->Fill(dt);
    runbr->GetEntry(i);
    trigbr->GetEntry(i);
    for(int j = lastpulse; j < pulse->GetEntries(); j++){
      prunbr->GetEntry(j);
      if(prun > run) break;
      ptrigbr->GetEntry(j);
      if(prun == run && ptrig > trig) break;

      if(prun == run && ptrig == trig){
        vetoed->Fill(dt);
        lastpulse = j;
        goto end;
      }
    }
    end: 1;
  }
}
