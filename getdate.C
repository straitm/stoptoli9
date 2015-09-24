int getdate(const int run)
{
  TChain t("data");
  t.Add(Form("/cp/s4/dchooz/cheetah/prod-08-05_p01_v2/reduced.Run%07d_Seq010.root", run));
  int date[6];
  t.SetBranchAddress("date", date);
  t.GetEntry(0);
  return date[0]*100 + date[1];
}
