float getqiv(const int run, const int trig)
{
  TChain t("data");
  t.Add(Form("/cp/s4/dchooz/cheetah/prod-08-05_p01_v2/reduced.Run%07d_Seq010.root", run));
  float ans;
  t.SetBranchAddress("fido_qiv", &ans);
  t.GetEntry(trig);
  return ans;
}
