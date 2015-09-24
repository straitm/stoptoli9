#include "time.h"
#include "TTree.h"
#include "TFile.h"

void getrunlivetime(const int run)
{
  TFile rf(Form("/cp/s4/dchooz/cheetah/prod-08-05_p01_v2/reduced.Run%07d_Seq010.root", run), "read");
  TFile rf11(Form("/cp/s4/dchooz/cheetah/seq11/reduced.Run%07d_Seq011.root", run), "read");
  TTree * data = (TTree*)rf.Get("data");
  if(!data || data->IsZombie()) data = (TTree*)rf11.Get("data");

  double trgtime;
  data->SetBranchAddress("trgtime", &trgtime);
  data->GetEntry(data->GetEntries()-1);

  printf("%.9lf\n", trgtime/1e9);

  delete data;
  rf.Close();
  rf11.Close();
}
