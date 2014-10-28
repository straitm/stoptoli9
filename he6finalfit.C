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

/*
  //TCanvas c0;
  //t->Draw("dy**2+dx**2 >> p1(40, 0, 3240000)", Form("%s && dz > -1000 && dz < 800 && e > 1.5 && e < 2.0 && dt/1000 > 1", cut ), "e");
  //t->Draw("dy**2+dx**2 >> p2(40, 0, 3240000)", Form("%s && dz > -1000 && dz < 800 && e > 2.0 && e < 2.5 && dt/1000 > 1", cut ), "samee");
  //t->Draw("dy**2+dx**2 >> p3(40, 0, 3240000)", Form("%s && dz > -1000 && dz < 800 && e > 2.5 && e < 3.1 && dt/1000 > 1", cut ), "samee");
  //t->Draw("dy**2+dx**2 >> p5(40, 0, 3240000)", Form("%s && dz > -1000 && dz < 800 && e > 3.1 && e < 3.5 && dt/1000 > 1", cut ), "samee");
  //return;

  TCanvas * c1 = new TCanvas;
  c1->SetLogy();

  t->Draw("dt/1000 >> hfit(3000, 0.001, 100)", Form("%s && e > 2.8 && e < 10", cut));
  TH1D * hfit = (TH1D*)gROOT->FindObject("hfit");

  TF1 ee("ee", Form("%f*([0]*exp(-x*log(2)/0.0202) + "
               "%f*[1]/[2]*log(2)*exp(-x*log(2)/[2]) + "
               "[3]*exp(-x*log(2)/7.13) + "
               "[4])", hfit->GetBinWidth(1), eff), 0, 100);

  ee.SetParameters(1, 1, 0.801, 1, 1);
  ee.FixParameter(2, 0.801);

  ee.FixParameter(3, 0);
  
  ee.FixParameter(1, 0);
  hfit->Fit("ee", "l");
  const double nohe8chi2 = gMinuit->fAmin;

  ee.ReleaseParameter(1);
  ee.SetParLimits(1, 0, 100000);
  hfit->Fit("ee", "le");
  const double whe8chi2 = gMinuit->fAmin;

  printf("%sEfficiency: %f%s\n", RED, eff, CLR);

  printf("%ssignificance of something: %f%s\n", RED,
    nohe8chi2 < whe8chi2? 0 : sqrt(2)*sqrt(nohe8chi2 - whe8chi2), CLR);

  t->Draw("dt/1000 >> hdisp(100, 0.1, 100)", Form("%s && e > 2.8 && e < 10", cut), "e");
  TH1D * hdisp = (TH1D*)gROOT->FindObject("hdisp");

  TF1 * eedisp = (TF1 *)ee.Clone("eedisp");
  eedisp->SetNpx(400);
  eedisp->SetLineColor(kRed);

  int tomult[4] = { 0, 1, 3, 4};
  const double mult = hdisp->GetBinWidth(1)/hfit->GetBinWidth(1);
  for(int i = 0; i < 4; i++)
    eedisp->SetParameter(tomult[i], eedisp->GetParameter(tomult[i])*mult);
  eedisp->Draw("same");

  TF1 b12("b12", Form("%f*[0]*exp(-x*log(2)/0.0202)", hfit->GetBinWidth(1)), 0, 100);
  TF1 he6("he6", Form("%f*%f*[0]/[2]*log(2)*exp(-x*log(2)/0.801)", hfit->GetBinWidth(1), eff), 0, 100);
  TF1 n16("n16", Form("%f*[0]*exp(-x*log(2)/7.13)", hfit->GetBinWidth(1))  , 0, 100);
  TF1 acc("acc", Form("%f*[0]", hfit->GetBinWidth(1)), 0, 100);

  b12.SetNpx(400);

  TF1 * parts[4] = { &b12, &he6, &n16, &acc };


  b12.SetParameter(0, eedisp->GetParameter(0));
  he6.SetParameter(0, eedisp->GetParameter(1));
  n16.SetParameter(0, eedisp->GetParameter(3));
  acc.SetParameter(0, eedisp->GetParameter(4));

  hdisp->GetYaxis()->SetRangeUser(acc.GetParameter(0)/50,
                                  acc.GetParameter(0)*10);

  for(int i = 0; i < 4; i++){
    parts[i]->SetLineStyle(7);
    parts[i]->SetLineWidth(2);
    parts[i]->Draw("Same");
  } 

  const double Nfound = ee.GetParameter(1);
  const double Nerrup = Nfound * gMinuit->fErp[1]/ee.GetParameter(1);
  const double Nerrlo = Nfound * gMinuit->fErn[1]/ee.GetParameter(1);

  printf("%sN found: %f +%f %f%s\n", RED, Nfound, Nerrup, Nerrlo, CLR);

  const double captures = 139. * 489.509;

  const double toprob = 1./captures/eff;

  printf("%sProb: %g +%g %g%s\n", 
      RED, toprob*Nfound, toprob*Nerrup, toprob*Nerrlo, CLR);
*/
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

/*
  for(int i = 0; i < 28.0e5; i++){ const double et=betan16l->GetRandom()        ; spec->Fill(gRandom->Gaus(et, sqrt(et)*0.077));}
  for(int i = 0; i < 66.2e5; i++){ const double et=(gRandom->Rndm() < 0.6)*betan16h1->GetRandom()+6.0477; spec->Fill(gRandom->Gaus(et, sqrt(et)*0.077));}
  for(int i = 0; i <  4.8e5; i++){ const double et=(gRandom->Rndm() < 0.6)*betan16h2->GetRandom()+7.11515;spec->Fill(gRandom->Gaus(et, sqrt(et)*0.077));}
  for(int i = 0; i < 1.06e5; i++){ const double et=(gRandom->Rndm() < 0.6)*betan16h3->GetRandom()+8.8719; spec->Fill(gRandom->Gaus(et, sqrt(et)*0.077));}
*/

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

/*
  ehistbg->Sumw2();
  ehistsig->Sumw2();

  ehistbg ->Draw("hist");
  ehistsig->Draw("esame");

  TH1D * ehistsub = (TH1D*) ehistsig->Clone("ehistsub");
  ehistsub->Add(ehistbg, -1);
  ehistsub->SetLineColor(kRed);
  ehistsub->Draw("samee");
  ehistbg ->GetYaxis()->SetRangeUser(ehistsub->GetMinimum()-1, ehistbg->GetMaximum());
*/
/*
  TCanvas c3;
  gMinuit->fUp = nohe8chi2 - whe8chi2 + 2.3/2;

  gMinuit->Command("Set print 0");
  gMinuit->Command("mncont 2 5 200");
  TGraph * ninty_1d = (TGraph*)((TGraph*)gMinuit->GetPlot())->Clone();
  ninty_1d->Draw("al");
  double lim = 0;
  for(int i = 0; i < ninty_1d->GetN(); i++)
    if(ninty_1d->GetX()[i] > lim)
      lim = ninty_1d->GetX()[i];
  printf("%sHe-6 90%% CL upper limit events: %.3lf %s\n", RED, lim, CLR);
  printf("%sHe-6 90%% CL upper limit prob: %.3lg %s\n", RED, lim*toprob, CLR);
*/
}
