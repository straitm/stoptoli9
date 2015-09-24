#include "time.h"
#include "TTree.h"
#include "TFile.h"

void trig2unixtime(const int run, const int trig)
{
  TFile rf(Form("/cp/s4/dchooz/cheetah/prod-08-05_p01_v2/reduced.Run%07d_Seq010.root", run), "read");
  TTree * data = (TTree*)rf.Get("data");

  int date[6];
  data->SetBranchAddress("date", date);
  double trgtime;
  data->SetBranchAddress("trgtime", &trgtime);
  data->GetEntry(trig);

  tm t;
  memset(&t, 0, sizeof(t));
  t.tm_year = date[0]-1900;
  t.tm_mon = date[1]-1;
  t.tm_mday = date[2];
  t.tm_hour = date[3];
  t.tm_min = date[4];
  t.tm_sec = date[5];

  printf("%d\n",  mktime(&t) - 6 * 3600 + int(trgtime/1e9));

  delete data;
  rf.Close();
}
