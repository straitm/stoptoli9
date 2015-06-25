
#include "TTree.h"
#include "TFile.h"

Float_t getivdedxslant(const int run, const int trig)
{
  static TTree * fidot = NULL;
  static TFile * fidof = NULL;
  static int lastrun = 0;
  static TBranch * br[4];
  static float ids_entr_z, id_ivlen, id_buflen, fido_qiv;
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
    fidot->SetBranchAddress("ids_entr_z", &ids_entr_z);
    fidot->SetBranchAddress("id_ivlen", &id_ivlen);
    fidot->SetBranchAddress("id_buflen", &id_buflen);
    fidot->SetBranchAddress("fido_qiv", &fido_qiv);
    br[0] = fidot->GetBranch("ids_entr_z");
    br[1] = fidot->GetBranch("id_ivlen");
    br[2] = fidot->GetBranch("id_buflen");
    br[3] = fidot->GetBranch("fido_qiv");
  }

  for(int i = 0; i < 4; i++) br[i]->GetEntry(trig);
  fidot->GetEntry(trig);

  lastrun = run;

  old->cd();

  if(id_ivlen <= 0) return 0;
  return ids_entr_z + 62*fido_qiv/(id_ivlen - id_buflen);
}
