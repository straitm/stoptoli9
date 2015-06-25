#include "TTree.h"
#include "TFile.h"

float getredchi2(const int run, const int trig)
{
  static TTree * fidot = NULL;
  static TFile * fidof = NULL;
  static int lastrun = 0;
  static TBranch * br[3];
  static float ids_chi2;
  static int nidtubes, nivtubes;
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
    fidot->SetBranchAddress("ids_chi2", &ids_chi2);
    fidot->SetBranchAddress("nidtubes", &nidtubes);
    fidot->SetBranchAddress("nivtubes", &nivtubes);
    br[0] = fidot->GetBranch("ids_chi2");
    br[1] = fidot->GetBranch("nidtubes");
    br[2] = fidot->GetBranch("nivtubes");
  }

  for(int i = 0; i < 4; i++) br[i]->GetEntry(trig);
  fidot->GetEntry(trig);

  lastrun = run;

  old->cd();

  if(nidtubes + nivtubes - 6 < 0) return 0;
  return ids_chi2/(nidtubes + nivtubes - 6);
}
