// Given a run number and a prompt trigger number,
// print out the prompt energy and Ps variables, as available
void moreli9vars(int run, int prompt, bool red = false)
{
  if(red){
    TFile rf(Form(
      "/cp/s4/dchooz/cheetah/prod-08-05_p01_v2/reduced.Run%07d_Seq010.root",
      run));

    TTree * r = (TTree *)rf.Get("data");

    r->Scan("ctEvisID:deltaT", Form("trgId == %d", prompt));
  }
  else{
    TFile f(Form(
      "/cp/s4/strait/dceuhcandps/DCRun%07d_RecoPSM_v17_EUppHCand_base.root.prompt.root",
      run));

    TFile fd(Form(
      "/cp/s4/strait/dceuhcandps/DCRun%07d_RecoPSM_v17_EUppHCand_base.root.delay.root",
      run));

    TTree * g = (TTree *)f.Get("GlobalInfoTree");
    TTree * e = (TTree *)f.Get("LightEUppInfoTree");
    TTree * gd = (TTree *)fd.Get("GlobalInfoTree");
    gd->SetName("gd");
    TTree * ed = (TTree *)fd.Get("LightEUppInfoTree");
    ed->SetName("ed");
    g->AddFriend(e);
    g->AddFriend(gd);
    g->AddFriend(ed);

    if(g->GetEntries( Form("fEventID == %d", prompt) ))
      g->Scan("fEvisID:ed.fEvisID:foPsdtCo:foPsdtCs:gd.fTrigTime-fTrigTime",
              Form("fEventID == %d", prompt));
    else
      printf(
          "************************************************\n"
          "*    Row   *   fEvisID *  foPsdtCo *  foPsdtCs *\n"
          "*************************************************\n"
          "*       90 * -1 * -1 * -1 *\n"
          "*************************************************\n"
          "==> 1 selected entry\n");
  }
}
