#include "TTree.h"
#include "TFile.h"

double getfidotstop(const int run, const int trig)
{
  static TTree * fidot = NULL;
  static TFile * fidof = NULL;
  static int lastrun = 0;
  static TBranch * br[2];
  static float ids_t0, ids_ivlen;
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
    fidot->SetBranchAddress("ids_t0", &ids_t0);
    fidot->SetBranchAddress("ids_ivlen", &ids_ivlen);
    br[0] = fidot->GetBranch("ids_t0");
    br[1] = fidot->GetBranch("ids_ivlen");
  }

  br[0]->GetEntry(trig);
  br[1]->GetEntry(trig);

  lastrun = run;

  old->cd();

  return ids_ivlen/300 + ids_t0;
}
