#include "TGraph.h"
#include "TChain.h"
#include "TRandom.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TMinuit.h"
#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "stdio.h"

double bgerr[4] = {0};
double abg[4][120], asig[4][120], an16[4][120], ali8[120], ahe6[120];

TH2D * ehistsig = NULL, * ehistbg = NULL;
TH1D * li8spec = NULL, * he6spec = NULL, * n16spec[4] = {NULL};

const double mus[4] = {0.45120, 0.44157, 0.49416, 2.2830};

double bamacorrxy(const double xy, const double e)
{
  return xy * (0.970 + 0.013*e);
}

double bamacorrz(const double z, const double e)
{
  return z  * (0.985 + 0.005*e);
}

void fcn(int & npar, double * gin, double & like, double *par, int flag)
{
  double mnbgnorm[4] = { par[0], par[1], par[2], par[3] };
  double mnli8norm = par[4];
  double mnhe6norm = par[5];
  double mnn16norm[4] = { par[6], par[7], par[8], par[9] };

  like = 0;
  for(int j = 0; j < 4; j++){
    bool datahasstarted = false;
    for(int i = 6; i < 120; i++){
      double model = mnbgnorm[j]*abg[j][i]
                   + mnn16norm[j]*an16[j][i]
                   + mnli8norm*mus[j]*ali8[i]
                   + mnhe6norm*mus[j]*ahe6[i];
      double data = asig[j][i];
      if(datahasstarted || data > 0) datahasstarted = true;
      else continue;
      like += model - data;
      if(data > 0 && model > 0) like += data*log(data/model);
    } 
    like += pow((1-mnbgnorm[j])/bgerr[j], 2);
  }
  // pull term for Li-8.  However, this doesn't let the Li-8
  // part of the fit account for N-16, which is useful.
  //like += pow((0.00229277 - mnli8norm)/1.79122e-04, 2);
  like *= 2;
}

int classi(const double x, const double y, const double z)
{
  const double r2 = x*x+y*y;
  const double r = sqrt(r2);
  const double az = abs(z);
  if(r2 > 1150*1150 || az > 1229 + 0.03*(1150-r)) return 3;
  if(r2 > 1000*1000 || az > 1000) return 2;
  if(r2 > 800*800 || az > 800) return 1;
  return 0;
}

void he6finalfit()
{
  TH1::SetDefaultSumw2();

  const char * const RED     = "\033[31;1m"; // bold red
  const char * const CLR      = "\033[m"    ; // clear

  TCanvas * c2 = new TCanvas();
  c2->Divide(1, 4);
  c2->cd(1);
  c2->ToggleEventStatus();

  TFile fiel("/cp/s4/strait/fullfido-100s-0-25MeV-20141022.root");
  TTree * t = (TTree *) fiel.Get("t");

  const double r2cut = 800, zcut = 800;

  const bool tgc = true, targ = false;

  if(tgc && targ){
    fprintf(stderr, "Danger Will Robinson!\n");
    return;
  }

  const double distcut = 200;

  char cut[1000];
  if(tgc)
    snprintf(cut, 999, 
      "b12like < 0.4 && !earlymich && ttlastvalid > 0.1 && ttlastmuon>1"
      "&& miche < 12 && dist < %f && timeleft > 100e3", distcut);
  else if(targ)
    snprintf(cut, 999, 
      "dx**2+dy**2 < %f**2 && abs(dz) < %f &&"
      "b12like < 0.4 && !earlymich && ttlastvalid > 0.1 && ttlastmuon>1"
      "&& miche < 12 && dist < %f && timeleft > 100e3", 1150., 1229., distcut);
  else
    snprintf(cut, 999, 
      "dx**2+dy**2 < %f**2 && abs(dz) < %f &&"
      "b12like < 0.4 && !earlymich && ttlastvalid > 0.1 && ttlastmuon>1"
      "&& miche < 12 && dist < %f && timeleft > 100e3", r2cut, zcut, distcut);

  const double distcuteff = (distcut == 400?0.948:distcut == 300?0.852:distcut == 200?0.565:distcut==159?0.376:100000);
    
  const double bglow = 7.13*5, bghigh = 100, siglow = 0.3, sighigh = 0.801*2;
  const double eff = 1
    * 0.981 // subsequent muons
    * 0.977 // previous muons
    * distcuteff // delta r
    * 0.9709 // 100s from end of run
    * 0.986 // ttlastvalid
    * 0.96 // ttlastmuon
    * 0.9994 // b12like
    * (exp(-siglow*log(2)/0.801) - exp(-sighigh*log(2)/0.801))
  ;

  const double li8eff = 1
    * 0.981 // subsequent muons
    * 0.977 // previous muons
    * distcuteff // delta r
    * 0.9709 // 100s from end of run
    * 0.986 // ttlastvalid
    * 0.96 // ttlastmuon
    * 0.9994 // b12like
    * (exp(-siglow*log(2)/0.8399) - exp(-sighigh*log(2)/0.8399))
  ;

#include "ehistbg.C"
#include "ehistsig.C"
 
  if(!ehistbg){ 
    ehistbg = new TH2D("ehistbg", "", 4, 0, 4, 120, 0, 15);
    fprintf(stderr, "regenerating ehistbg\n");
    t->Draw("e:classi(dx, dy, dz) >> ehistbg ", Form("%s && dt > %f", cut, bglow*1e3));

    ehistbg->Scale((sighigh-siglow)/(bghigh-bglow));
  }

  for(int j = 0; j < 4; j++){
    bgerr[j] = 1/sqrt(ehistbg->Integral(j+1,j+1));
    if((bghigh-bglow)/15 > 1) bgerr[j] *= sqrt((bghigh-bglow)/15); // XXX fudge
  }

  if(!ehistsig){ 
    ehistsig = new TH2D("ehistsig", "", 4, 0, 4, 120, 0, 15);
    fprintf(stderr, "regenerating ehistsig\n");
    t->Draw("e:classi(dx, dy, dz) >> ehistsig", Form("%s && dt > %f && dt < %f", cut, siglow*1e3, sighigh*1e3));
  }
  
  TF1 betahe6("betahe6", "(x+0.511)*sqrt((x+0.511)**2 - 0.511**2) * (3.5076 - x)**2", 0, 3.51);
  he6spec = new TH1D("he6spec", "", 120, 0, 15);
  li8spec = new TH1D("li8spec", "", 120, 0, 15);
  for(int j = 0; j < 4; j++)
    n16spec[j] = new TH1D(Form("n16spec%d", j), "", 120, 0, 15);
  for(int i = 0; i < 1e7; i++){
    const double et=betahe6.GetRandom();
    he6spec->Fill(gRandom->Gaus(et, sqrt(et*pow(0.077, 2)
                                     +et*et*pow(0.018, 2)
                                          + pow(0.017, 2)) ));
  }
  
  TChain li8mc("data");
  li8mc.Add("/cp/s4/strait/li8mc/*.root");
  li8mc.Draw("ctEvisID >> li8spec", "bamacorrxy(ctX[0], ctEvisID)**2 + bamacorrxy(ctX[1], ctEvisID)**2 < (1150*1.2)**2 && abs(bamacorrz(ctX[2], ctEvisID)) < 1229*1.2");

  TChain n16mc("data");
  n16mc.Add("/cp/s4/strait/n16mc/*.root");
  n16mc.Draw("ctEvisID >> n16spec0", "bamacorrxy(ctX[0], ctEvisID)**2 + bamacorrxy(ctX[1], ctEvisID)**2 < 800**2 && abs(bamacorrz(ctX[2], ctEvisID)) < 800");
  n16mc.Draw("ctEvisID >> n16spec1", "bamacorrxy(ctX[0], ctEvisID)**2 + bamacorrxy(ctX[1], ctEvisID)**2 > 800**2 && "
                                     "bamacorrxy(ctX[0], ctEvisID)**2 + bamacorrxy(ctX[1], ctEvisID)**2 < 1000**2 && "
                                     "abs(bamacorrz(ctX[2], ctEvisID)) > 800 && abs(bamacorrz(ctX[2], ctEvisID)) < 1000");
  n16mc.Draw("ctEvisID >> n16spec2", "bamacorrxy(ctX[0], ctEvisID)**2 + bamacorrxy(ctX[1], ctEvisID)**2 > 1000**2 && "
                                     "bamacorrxy(ctX[0], ctEvisID)**2 + bamacorrxy(ctX[1], ctEvisID)**2 < 1150**2 && "
                                     "abs(bamacorrz(ctX[2], ctEvisID)) > 800 && abs(bamacorrz(ctX[2], ctEvisID)) < 1229");
  n16mc.Draw("ctEvisID >> n16spec3", "bamacorrxy(ctX[0], ctEvisID)**2 + bamacorrxy(ctX[1], ctEvisID)**2 > 1150**2 && abs(bamacorrz(ctX[2], ctEvisID)) > 1229");

  he6spec->Scale(ehistsig->Integral()/he6spec->Integral());
  li8spec->Scale(ehistsig->Integral()/li8spec->Integral());
  for(int j = 0; j < 4; j++)
    n16spec[j]->Scale(ehistsig->Integral()/n16spec[j]->Integral());

  for(int i = 0; i < 120; i++){
    ali8[i] = li8spec->GetBinContent(i+1);
    ahe6[i] = he6spec->GetBinContent(i+1);
    for(int j = 0; j < 4; j++){
      abg[j][i]  =ehistbg ->GetBinContent(j+1, i+1);
      asig[j][i] =ehistsig->GetBinContent(j+1, i+1);
      an16[j][i] = n16spec[j]->GetBinContent(i+1);
    }
  }


  TMinuit * mn = new TMinuit(9);
  mn->SetFCN(fcn);
  int err;
  mn->mnparm(0, "bgtargin",   1, 0.01, 0,  0, err);
  mn->mnparm(1, "bgtargmid",  1, 0.01, 0,  0, err);
  mn->mnparm(2, "bgtargout",  1, 0.01, 0,  0, err);
  mn->mnparm(3, "bggc",       1, 0.01, 0,  0, err);
  mn->mnparm(4, "li8",        1, 0.001, 0.0001, 1, err);
  mn->mnparm(5, "he6",        1, 0.01, 0, 1, err);
  mn->mnparm(6, "n16targin",  1, 0.001, 0, 0.1, err);
  mn->mnparm(7, "n16targmid", 1, 0.001, 0, 0.1, err);
  mn->mnparm(8, "n16targout", 1, 0.001, 0, 0.1, err);
  mn->mnparm(9, "n16gc",      1, 0.001, 0, 0.1, err);
  printf("err? %d\n", err);

  mn->Command("SET PAR 6 0");
  mn->Command("FIX 6");
  mn->Command("MIGRAD");

  const double chi2nohe6 = mn->fAmin;

  mn->Command("RELEASE 6");
  mn->Command("SET PAR 6 0.1");
  mn->Command("MIGRAD");
  mn->Command("MINOS 10000 6");

  const double chi2he6 = mn->fAmin;

  if(chi2nohe6 - chi2he6 > 0)
    printf("%sSignificance of He-6: %f%s\n",
           RED, sqrt(chi2nohe6 - chi2he6), CLR);

  const double upforlim = 2.71 + chi2nohe6 - chi2he6;

  double dum, bgnorm[4], li8norm, he6norm, n16norm[4];

  for(int j = 0; j < 4; j++) mn->GetParameter(j, bgnorm[j], dum);
  mn->GetParameter(4, li8norm, dum);
  mn->GetParameter(5, he6norm, dum);
  for(int j = 0; j < 4; j++) mn->GetParameter(j+6, n16norm[j], dum);

  double li8normerr, he6normerr;
  mn->mnerrs(4, dum, dum, li8normerr, dum);
  mn->mnerrs(5, dum, dum, he6normerr, dum);

  const double captures = 489.509 *
                  (tgc ?367.:
                   targ?139.:
                        139. * pow(r2cut/1150., 2) * zcut/1229.);

  const double toprob = 1./captures/eff;

  const double nhe6 = he6spec->Integral() * he6norm *
                      (mus[0] + mus[1] + mus[2] + mus[3]);
  const double signhe6 = he6spec->Integral() * he6normerr *
                      (mus[0] + mus[1] + mus[2] + mus[3]);

  printf("%sefficiency: %f%s\n", RED, eff, CLR);
  printf("%sraw n He-6: %f +- %f%s\n", RED, nhe6, signhe6, CLR);
  printf("%sHe-6 prob: %f +- %f%s\n", RED, nhe6*toprob,signhe6*toprob, CLR);

  const double nli8 = li8spec->Integral() * li8norm *
                      (mus[0] + mus[1] + mus[2] + mus[3]);
  const double signli8 = li8spec->Integral() * li8normerr *
                      (mus[0] + mus[1] + mus[2] + mus[3]);

  printf("%sLi-8 prob: %f +- %f%s\n",RED,nli8/captures/li8eff,
                                     signli8/captures/li8eff, CLR);

  mn->fUp = upforlim;
  mn->Command("MINOS 10000 6");
 
  double he6up90;
  mn->GetParameter(5, he6norm, dum);
  mn->mnerrs(5, he6up90, dum, dum, dum);

  printf("%s 90%% upper limit: %f+%f = %f%s\n", RED, he6norm, he6up90, he6norm+he6up90, CLR);
  
  mn->Command("set print -3");
  mn->Command("Fix 6");
  double sump = 0;
  for(int i =0;i<1000; i++){
    mn->Command(Form("set  par 6 %f", i*0.01*he6norm));
    mn->Command("Migrad");
    double p = exp(-mn->fAmin);
    sump += p;
  }

  double sump2 = 0;
  int i;
  for(i =0;i<1000; i++){
    mn->Command(Form("set par 6 %f", i*0.01*he6norm));
    mn->Command("Migrad");
    double p = exp(-mn->fAmin);
    sump2 += p;
    if(sump2>sump*0.9){ printf("%s 90% bays lim: %f (best = %f)%s\n", RED, i*0.1*he6norm, he6norm, CLR); break;}
  }
  mn->Command("set print 0");

  double he6norm90 = i*0.01*he6norm;

  double n90he6 = 0;
  for(int j = 1; j <= 4; j++){
    c2->cd(j);

    TH1D * ehistbg_p=ehistbg->ProjectionY(Form("ehistbg_%d", j), j, j);
    TH1D * ehistsig_p=ehistsig->ProjectionY(Form("ehistsig_%d", j), j, j);
    TH1D * li8spec_p = (TH1D * )li8spec->Clone(Form("li8spec_%d", j));
    TH1D * he6spec_p = (TH1D * )he6spec->Clone(Form("he6spec_%d", j));
    TH1D * n16spec_p = (TH1D * )n16spec[j-1]->Clone(Form("n16specp_%d", j));

    ehistbg_p->Scale(bgnorm[j-1]);
    n16spec_p->Scale(n16norm[j-1]);
    li8spec_p->Scale(li8norm*mus[j-1]);
    he6spec_p->Scale(he6norm*mus[j-1]);


    ehistsig_p->SetLineColor(kRed);
    ehistsig_p->SetMarkerColor(kRed);
    ehistsig_p->GetXaxis()->SetRange(6, 120);
    ehistsig_p->Draw("e");

    ehistbg_p->Draw("same");

    n16spec_p->SetLineColor(kOrange);
    n16spec_p->Draw("samehist");

    TH1D * bg_plus_li8 = (TH1D*)ehistbg_p->Clone(Form("bg_plus_li8_%d", j));
    bg_plus_li8->Add(li8spec_p);
    bg_plus_li8->SetLineColor(kBlue);
    bg_plus_li8->Draw("same");

    TH1D * bg_plus_li8_n16 = (TH1D*)bg_plus_li8->Clone(Form("bg_plus_li8_n16_%d", j));
    bg_plus_li8_n16->Add(n16spec_p);
    bg_plus_li8_n16->SetLineColor(kViolet);
    bg_plus_li8_n16->Draw("same");

    TH1D * bg_plus_li8_n16_he6 = (TH1D*)bg_plus_li8_n16->Clone(Form("bg_plus_li8_n16_he6_%d", j));
    bg_plus_li8_n16_he6->Add(he6spec_p);
    bg_plus_li8_n16_he6->SetLineColor(kGreen+2);
    bg_plus_li8_n16_he6->Draw("same");

    TH1D * he6demo = (TH1D*) he6spec_p->Clone(Form("he6demo%d", j));
    he6demo->Scale(he6norm90/he6norm);
    he6demo->SetLineStyle(kDashed);
    he6demo->Draw("same");

    TH1D * alldemo = new TH1D(Form("alldemo%d", j), "", 120, 0, 15);
    alldemo->Add(ehistbg_p);
    alldemo->Add(li8spec_p);
    alldemo->Add(he6demo);
    alldemo->SetLineColor(kGreen+2);
    alldemo->SetLineStyle(kDashed);
    alldemo->Draw("same");

    n90he6 += he6demo->Integral(1, he6demo->GetNbinsX());

    /*alldemo->Rebin(4);
    he6demo->Rebin(4);
    bg_plus_li8_plus_he6->Rebin(4);
    bg_plus_li8->Rebin(4);
    ehistbg_p->Rebin(4);
    ehistsig_p->Rebin(4); */
  }

  printf("%s 90%% upper limit raw n He-6: %f%s\n", RED, n90he6, CLR);
  printf("%s 90%% upper limit He-6 prob: %f%s\n", RED, n90he6*toprob, CLR);

  c2->SaveAs("he6.pdf");
  c2->SaveAs("he6.C");

  ehistsig->SaveAs("tmp-ehistsig.C");
  ehistbg->SaveAs("tmp-ehistbg.C");
}

