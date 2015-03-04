#include "TTree.h"
#include "TFile.h"

static TTree * t = NULL;
static TFile * f = NULL;

Float_t getm(const int run, const int trig, const int c)
{
  static int lastrun = 0;
  static Float_t ans[3];
  static TBranch * br[3];
  if(lastrun != run){
    fprintf(stderr, ".");
    if(t) delete t;
    if(f) delete f;
    f = new TFile(Form("/cp/s4/strait/fido_seq010/fido.%07d.root",run), "read");
    t = (TTree *) f->Get("RecoMuonFIDOInfoTree"); 
    t->SetMakeClass(1);
    t->SetBranchAddress("ids_end_x", ans + 0);
    t->SetBranchAddress("ids_end_y", ans + 1);
    t->SetBranchAddress("ids_end_z", ans + 2);
    br[0] = t->GetBranch("ids_end_x");
    br[1] = t->GetBranch("ids_end_y");
    br[2] = t->GetBranch("ids_end_z");
  }

  br[c]->GetEntry(trig);

  lastrun = run;

  return ans[c];
}

Float_t getmx(const int run, const int trig){ return getm(run,trig,0); }
Float_t getmy(const int run, const int trig){ return getm(run,trig,1); }
Float_t getmz(const int run, const int trig){ return getm(run,trig,2); }
