#include "TTree.h"
#include "TFile.h"

float getqid(const int run, const int trig)
{
  static TTree * fidot = NULL;
  static TFile * fidof = NULL;
  static int lastrun = 0;
  static TBranch * br[1];
  static float fido_qid;
  TDirectory * old = gDirectory;
  if(lastrun != run){
    fprintf(stderr, ".");
    if(fidot) delete fidot;
    if(fidof) delete fidof;
    fidof = new TFile(Form("/cp/s4/strait/fido_seq010/fido.%07d.root",run), "read");
    if(!fidof || fidof->IsZombie()) {lastrun = 0, fidot = NULL, fidof = NULL; return 0;}
    fidot = (TTree *) fidof->Get("RecoMuonFIDOInfoTree"); 
    if(!fidot){ lastrun = 0, fidot = NULL, fidof = NULL; return 0;}
    fidot->SetMakeClass(1);
    fidot->SetBranchAddress("fido_qid", &fido_qid);
    br[0] = fidot->GetBranch("fido_qid");
  }

  br[0]->GetEntry(trig);

  lastrun = run;

  old->cd();

  return fido_qid;
}
