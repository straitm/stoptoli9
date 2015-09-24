bool lightnoise(const float qrms, const float mqtq, const float rmsts,
                const float qdiff)
{
  if(mqtq > 0.12 || mqtq < 0) return true;
  if(qdiff > 30e3) return true;
  if(rmsts >= 36 && qrms >= 464 - 8*rmsts) return true;
  return false;
}


double getdelaye(int run, int prompt)
{
  TFile rf(Form("/cp/s4/dchooz/cheetah/prod-08-05_p01_v2/reduced.Run%07d_Seq010.root", run), "read");
  TTree * data = (TTree*)rf.Get("data");

  float qrms, mqtq, rmsts, qdiff, ctEvisID;
  double trgtime, ctqIV;
  data->SetBranchAddress("qrms", &qrms);
  data->SetBranchAddress("qdiff", &qdiff);
  data->SetBranchAddress("ctmqtqall", &mqtq);
  data->SetBranchAddress("ctrmsts", &rmsts);
  data->SetBranchAddress("ctEvisID", &ctEvisID);
  data->SetBranchAddress("trgtime", &trgtime);
  data->SetBranchAddress("ctqIV", &ctqIV);

  data->GetEntry(prompt);

  const double prompttime = trgtime;

  int delay = prompt;
  do{
    data->GetEntry(++delay);
  }while(lightnoise(qrms, mqtq, rmsts, qdiff) || ctEvisID <= 0.4 || ctEvisID > 20 || ctqIV > 30e3);
  
  return ctEvisID;
}
