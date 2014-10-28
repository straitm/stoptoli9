#include "TGraph.h"
#include "TRandom.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TMinuit.h"
#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "stdio.h"

double bgerr = 0;
double abg[60], ali8[60], ahe6[60], asig[60];

TH1D * ehistsig = NULL, * li8spec = NULL,
     * he6spec = NULL,  * ehistbg = NULL;

void fcn(int & npar, double * gin, double & like, double *par, int flag)
{
  const double bgnorm = par[0];
  const double li8norm = par[1];
  const double he6norm = par[2];

  like = 0;
  for(int i = 7; i < 60; i++){
    const double model = bgnorm*abg[i] + li8norm*ali8[i] + he6norm*ahe6[i];
    const double data = asig[i];
    like += model - data;
    if(data > 0) like += data*log(data/model);
  } 
  like *= 2;
  like += pow((1-bgnorm)/bgerr, 2);
}

void he6finalfit()
{
  const char * const RED     = "\033[31;1m"; // bold red
  const char * const CLR      = "\033[m"    ; // clear

  TFile fiel("/cp/s4/strait/fullfido-100s-1.5-25MeV-20141022.root");
  TTree * t = (TTree *) fiel.Get("t");

  const double r2cut = 1708, zcut = 1786;

  char cut[1000];
  snprintf(cut, 999, 
    "dx**2+dy**2 < %f**2 && abs(dz) < %f &&"
    "b12like < 0.4 && !earlymich && ttlastvalid > 0.1 && ttlastmuon > 1"
    "&& miche < 12 && dist < 400 && timeleft > 100e3", r2cut, zcut);

  const double bglow = 7.13*5, bghigh = 100, siglow = 0.5, sighigh = 0.801*3;
  const double eff = 1
    * 0.981 // subsequent muons
    * 0.977 // previous muons
    * 0.948 // delta r
    * 0.9709 // 100s from end of run
    * 0.986 // ttlastvalid
    * 0.96 // ttlastmuon
    * 0.9994 // b12like
    * (exp(-siglow*log(2)/0.801) - exp(-sighigh*log(2)/0.801))
  ;

  const double li8eff = 1
    * 0.981 // subsequent muons
    * 0.977 // previous muons
    * 0.948 // delta r
    * 0.9709 // 100s from end of run
    * 0.986 // ttlastvalid
    * 0.96 // ttlastmuon
    * 0.9994 // b12like
    * (exp(-siglow*log(2)/0.8399) - exp(-sighigh*log(2)/0.8399))
  ;
  TCanvas * c2 = new TCanvas;

//#include "ehistbg.C"
//#include "ehistsig.C"
 
  if(!(ehistbg = (TH1D*)gROOT->FindObject("ehistbg"))){ 
    ehistbg = new TH1D("ehistbg", "", 60, 0, 15);
    fprintf(stderr, "regenerating ehistbg\n");
    t->Draw("e >> ehistbg ", Form("%s && dt > %f", cut, bglow*1e3), "hist");

    ehistbg->Scale((sighigh-siglow)/(bghigh-bglow));
  }

  if((bghigh-bglow)/15 > 1) bgerr *= sqrt((bghigh-bglow)/15); // XXX fudge
  bgerr = 1/sqrt(ehistbg->Integral());

  if(!(ehistsig = (TH1D*)gROOT->FindObject("ehistsig"))){ 
    ehistsig = new TH1D("ehistsig", "", 60, 0, 15);
    fprintf(stderr, "regenerating ehistsig\n");
    t->Draw("e >> ehistsig", Form("%s && dt > %f && dt < %f", cut, siglow*1e3, sighigh*1e3), "esame");
  }
  
  TF1 betahe6("betahe6", "(x+0.511)*sqrt((x+0.511)**2 - 0.511**2) * (3.5076 - x)**2", 0, 3.51);
  TF1 betali8("betali8", "(x+0.511)*sqrt((x+0.511)**2 - 0.511**2) * (16.00517 - x)**2", 0, 16.00517);
  TF1 betan16l("betan16l", "(x+0.511)*sqrt((x+0.511)**2 - 0.511**2) * (10.4191 - x)**2", 0, 10.4191);
  TF1 betan16h1("betan16h1", "(x+0.511)*sqrt((x+0.511)**2 - 0.511**2) * (4.370 - x)**2", 0, 4.370);
  TF1 betan16h2("betan16h2", "(x+0.511)*sqrt((x+0.511)**2 - 0.511**2) * (3.3022- x)**2", 0,3.3022);
  TF1 betan16h3("betan16h3", "(x+0.511)*sqrt((x+0.511)**2 - 0.511**2) * (1.5472- x)**2", 0,1.5472);
  he6spec = new TH1D("he6spec", "", 60, 0, 15);
  li8spec = new TH1D("li8spec", "", 60, 0, 15);
  for(int i = 0; i < 1e6; i++){
    const double et=betahe6.GetRandom();
    he6spec->Fill(gRandom->Gaus(et, sqrt(et)*0.077));
  }
  for(int i = 0; i < 1e6; i++){
    const double et=betali8.GetRandom();
    li8spec->Fill(gRandom->Gaus(et, sqrt(et)*0.077));
  }

  he6spec->Scale(ehistsig->Integral()/he6spec->Integral());
  li8spec->Scale(ehistsig->Integral()/li8spec->Integral());

  for(int i = 0; i < 60; i++){
    abg[i]  =ehistbg ->GetBinContent(i+1);
    asig[i] =ehistsig->GetBinContent(i+1);
    ali8[i] = li8spec->GetBinContent(i+1);
    ahe6[i] = he6spec->GetBinContent(i+1);
  }

  TMinuit * mn = new TMinuit(3);
  mn->SetFCN(fcn);
  int err;
  mn->mnparm(0, "bg",  1, 0.01, 0, 10, err);
  mn->mnparm(1, "li8", 1, 0.01, 0, 10, err);
  mn->mnparm(2, "he6", 1, 0.01, 0, 10, err);
  printf("err? %d\n", err);

  mn->Command("SET PAR 3 0");
  mn->Command("FIX 3");
  mn->Command("MIGRAD");

  const double chi2nohe6 = mn->fAmin;

  mn->Command("RELEASE 3");
  mn->Command("MIGRAD");
  mn->Command("MINOS 10000 3");

  const double chi2he6 = mn->fAmin;

  if(chi2nohe6 - chi2he6 > 0)
    printf("%sSignificance of He-6: %f%s\n", RED, sqrt(chi2nohe6 - chi2he6), CLR);

  const double upforlim = 6.25 + chi2nohe6 - chi2he6;

  double dum, bgnorm, li8norm, he6norm;
  mn->GetParameter(0, bgnorm, dum);
  mn->GetParameter(1, li8norm, dum);
  mn->GetParameter(2, he6norm, dum);

  double he6normerr;
  mn->mnerrs(2, dum, dum, he6normerr, dum);

  double li8normerr;
  mn->mnerrs(1, dum, dum, li8normerr, dum);

  ehistbg->Scale(bgnorm);
  li8spec->Scale(li8norm);
  he6spec->Scale(he6norm);

  ehistsig->Draw("e");
  ehistbg->Draw("same");
  li8spec->Draw("same");
  he6spec->Draw("same");

  TH1D * all = new TH1D("all", "", 60, 0, 15);
  all->Add(ehistbg);
  all->Add(li8spec);
  all->Add(he6spec);
  all->SetLineColor(kRed);
  all->Draw("same");

  const double captures = 139. * 489.509 * pow(r2cut/1150., 2) * zcut/1229.;

  const double toprob = 1./captures/eff;

  const double nhe6 = he6spec->Integral();
  const double signhe6 = he6spec->Integral() * he6normerr/he6norm;

  printf("%s efficiency: %f%s\n", RED, eff, CLR);
  printf("%s raw n He-6: %f +- %f%s\n", RED, nhe6, signhe6, CLR);
  printf("%s He-6 prob: %f +- %f%s\n", RED, nhe6*toprob, signhe6*toprob, CLR);

  const double nli8 = li8spec->Integral();
  const double signli8 = li8spec->Integral() * li8normerr/li8norm;

  printf("%s Li-8 prob: %f +- %f%s\n", RED, nli8/captures/li8eff, signli8/captures/li8eff, CLR);

  mn->fUp = upforlim;
  mn->Command("MINOS 10000 3");
 
  double he6up90;
  mn->GetParameter(2, he6norm, dum);
  mn->mnerrs(2, he6up90, dum, dum, dum);

  printf("%s 90%% upper limit: %f+%f = %f%s\n", RED, he6norm, he6up90, he6norm+he6up90, CLR);
  
  TH1D * he6demo = (TH1D*) he6spec->Clone("he6demo");
  he6demo->Scale((he6norm+he6up90)/he6norm);
  he6demo->SetLineStyle(kDashed);
  he6demo->Draw("same");

  const double n90he6 = he6demo->Integral();

  printf("%s 90%% upper limit raw n He-6: %f%s\n", RED, n90he6, CLR);
  printf("%s 90%% upper limit He-6 prob: %f%s\n", RED, n90he6*toprob, CLR);

  c2->SaveAs("he6.pdf");
  c2->SaveAs("he6.C");
}

