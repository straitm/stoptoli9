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
double abg[80], ali8[80], ahe6[80], asig[80];

TH1D * ehistsig, * li8spec, * he6spec, * ehistbg;

void fcn(int & npar, double * gin, double & like, double *par, int flag)
{
  const double bgnorm = par[0];
  const double li8norm = par[1];
  const double he6norm = par[2];

  like = 0;
  for(int i = 7; i < 80; i++){
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

  const double r2cut = 1150, zcut = 1229;

  const char * const cut = string(Form(
    "dx**2+dy**2 < 1150**2 && abs(dz) < 1229 &&"
    "b12like < 0.4 && !earlymich && ttlastvalid > 0.1 && ttlastmuon > 1 && miche < 12 && dist < 400 && timeleft > 100e3", r2cut, zcut)).c_str();


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

  TCanvas * c2 = new TCanvas;
 
  TF1 betahe6("betahe6", "(x+0.511)*sqrt((x+0.511)**2 - 0.511**2) * (3.5076 - x)**2", 0, 3.51);
  TF1 betali8("betali8", "(x+0.511)*sqrt((x+0.511)**2 - 0.511**2) * (16.00517 - x)**2", 0, 16.00517);
  TF1 betan16l("betan16l", "(x+0.511)*sqrt((x+0.511)**2 - 0.511**2) * (10.4191 - x)**2", 0, 10.4191);
  TF1 betan16h1("betan16h1", "(x+0.511)*sqrt((x+0.511)**2 - 0.511**2) * (4.370 - x)**2", 0, 4.370);
  TF1 betan16h2("betan16h2", "(x+0.511)*sqrt((x+0.511)**2 - 0.511**2) * (3.3022- x)**2", 0,3.3022);
  TF1 betan16h3("betan16h3", "(x+0.511)*sqrt((x+0.511)**2 - 0.511**2) * (1.5472- x)**2", 0,1.5472);
  he6spec = new TH1D("he6spec", "", 80, 0, 20);
  li8spec = new TH1D("li8spec", "", 80, 0, 20);
  for(int i = 0; i < 1e6; i++){ const double et=betahe6.GetRandom(); he6spec->Fill(gRandom->Gaus(et, sqrt(et)*0.077));}
  for(int i = 0; i < 1e6; i++){ const double et=betali8.GetRandom(); li8spec->Fill(gRandom->Gaus(et, sqrt(et)*0.077));}


  ehistbg = new TH1D("ehistbg", "", 80, 0, 20);
  ehistsig = new TH1D("ehistsig", "", 80, 0, 20);
  t->Draw("e >> ehistbg ", Form("%s && dt > %f", cut, bglow*1e3), "hist");
  t->Draw("e >> ehistsig", Form("%s && dt > %f && dt < %f", cut, siglow*1e3, sighigh*1e3), "esame");
  
  bgerr = 1/sqrt(ehistbg->Integral());
  if((bghigh-bglow)/20 > 1) bgerr *= sqrt((bghigh-bglow)/20); // XXX fudge
  ehistbg->Scale((sighigh-siglow)/(bghigh-bglow));

  he6spec->Scale(ehistsig->Integral()/he6spec->Integral());
  li8spec->Scale(ehistsig->Integral()/li8spec->Integral());

  for(int i = 0; i < 80; i++){
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

  mn->Command("FIX 1");
  mn->Command("MIGRAD");
  mn->Command("MINOS 10000 3");
  mn->Command("RELEASE 1");
  mn->Command("MINOS 10000 3");

  double dum, bgnorm, li8norm, he6norm;
  mn->GetParameter(0, bgnorm, dum);
  mn->GetParameter(1, li8norm, dum);
  mn->GetParameter(2, he6norm, dum);

  double he6normerr;
  mn->mnerrs(2, dum, dum, he6normerr, dum);


  ehistbg->Scale(bgnorm);
  li8spec->Scale(li8norm);
  he6spec->Scale(he6norm);


  ehistsig->Draw();
  ehistbg->Draw("esame");
  li8spec->Draw("esame");
  he6spec->Draw("esame");

  const double captures = 139. * 489.509;

  const double toprob = 1./captures/eff;

  const double nhe6 = he6spec->Integral();
  const double signhe6 = he6spec->Integral() * he6normerr/he6norm;

  printf("%s efficiency: %f%s\n", RED, eff, CLR);
  printf("%s raw n He-6: %f +- %f%s\n", RED, nhe6, signhe6, CLR);
  printf("%s He-6 prob: %f +- %f%s\n", RED, nhe6*toprob, signhe6*toprob, CLR);
}
