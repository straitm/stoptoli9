#include "TTree.h"
#include "TFile.h"

double gettimeid(const int run, const int trig)
{
  static TTree * datat = NULL;
  static TFile * dataf = NULL;
  static int lastrun = 0;
  static TBranch * br[1];
  static double timeid;
  TDirectory * old = gDirectory;
  if(lastrun != run){
    fprintf(stderr, ".");
    if(datat) delete datat;
    if(dataf) delete dataf;
    dataf = new TFile(Form("/cp/s4/dchooz/cheetah/prod-08-05_p01_v2/reduced.Run%07d_Seq010.root",run), "read");
    if(!dataf || dataf->IsZombie()) {lastrun = 0, datat = NULL, dataf = NULL; return 0;}
    datat = (TTree *) dataf->Get("data"); 
    if(!datat){ lastrun = 0, datat = NULL, dataf = NULL; return 0;}
    datat->SetBranchAddress("timeid", &timeid);
    br[0] = datat->GetBranch("timeid");
  }

  br[0]->GetEntry(trig);

  lastrun = run;

  old->cd();

  return timeid;
}
