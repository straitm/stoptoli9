#include "TTree.h"
#include "TFile.h"

void run2month(const int run)
{
  TFile rf(Form("/cp/s4/dchooz/cheetah/prod-08-05_p01_v2/reduced.Run%07d_Seq010.root", run), "read");
  TFile rf11(Form("/cp/s4/dchooz/cheetah/seq11/reduced.Run%07d_Seq011.root", run), "read");
  TTree * data = (TTree*)rf.Get("data");
  if(!data || data->IsZombie()) data = (TTree*)rf11.Get("data");

  int date[6];
  data->SetBranchAddress("date", date);
  double trgtime;
  data->SetBranchAddress("trgtime", &trgtime);
  data->GetEntry(0);

  printf("%04d%02d\n", date[0], date[1]);
  t.tm_year = date[0]-1900;

  delete data;
  rf.Close();
  rf11.Close();
}
