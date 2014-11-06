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

int ncut = 0;

double teff = 0, gceff = 0;

const int nrbins = 5;

double bgerr[nrbins] = {0};
double abg[nrbins][120], asig[nrbins][120], an16[nrbins][120], ali8[120], ahe6[120];

TH2D * ehistsig = NULL, * ehistbg = NULL;
TH1D * li8spec = NULL, * he6spec = NULL, * n16spec[5] = {NULL};

TH2D * mdispbg = NULL;

const double mus[nrbins] =
  {0.36274982, 0.3592311, 0.35758143, 0.30827676, 2.2830};

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
  double mnbgnorm[nrbins] = { par[0], par[1], par[2], par[3], par[4] };
  double mnli8norm = par[nrbins];
  double mnhe6norm = par[nrbins+1];
  double mnn16norm[nrbins] = { par[nrbins+2], par[nrbins+3], par[nrbins+4], par[nrbins+5], par[nrbins+6] };

  like = 0;
  for(int j = 0; j < nrbins; j++){
    bool datahasstarted = false;
    for(int i = 6; i < 120; i++){
      double model = mnbgnorm[j]*abg[j][i]
                   + mnn16norm[j]*an16[j][i]
                   + mnli8norm*mus[j]*ali8[i]
                   + mnhe6norm*mus[j]*ahe6[i]*(j < nrbins-1?teff:gceff);
      double data = asig[j][i];
      if(datahasstarted || data > 0) datahasstarted = true;
      else continue;
      like += model - data;
      if(data > 0 && model > 0) like += data*log(data/model);
    } 
    like += pow((1-mnbgnorm[j])/bgerr[j], 2);
  }
  // pull term for Li-8.
  //like += pow((0.00229277 - mnli8norm)/1.79122e-04, 2);
  like *= 2;
}

int classi(const double x, const double y, const double z)
{
  const double r2 = x*x+y*y;
  const double r = sqrt(r2);
  const double az = abs(z);
  if(r2 > 1150*1150 || az > 1229 + 0.03*(1150-r)) return 4;
  if(r2 > 1068.*1068. || az > 1068.) return 3;
  if(r2 > 933.*933. || az > 933.) return 2;
  if(r2 > 740*740 || az > 740) return 1;
  return 0;
}

void he6finalfit(const int ncutin = 0)
{
  ncut = ncutin;
  
  TH1::SetDefaultSumw2();

  const char * const RED     = "\033[31;1m"; // bold red
  const char * const CLR      = "\033[m"    ; // clear

  TCanvas * c2 = new TCanvas();
  c2->Divide(1, nrbins);
  c2->cd(1);
  c2->ToggleEventStatus();

  TFile fiel("/cp/s4/strait/fullfido-100s-0-25MeV-20141022.root");
  TTree * t = (TTree *) fiel.Get("t");

  const double distcut = 200;

  char cutnoncut[1000];
  snprintf(cutnoncut, 999, 
    "b12like < 0.4 && !earlymich && ttlastvalid > 0.1 && ttlastmuon>1"
    "&& miche < 12 && dist < %f && timeleft > 100e3", distcut);
  char cut[1000];
  snprintf(cut, 999, "n >= %d && n <= 3 && %s", ncut, distcut, cutnoncut);

  const double distcuteff = (distcut == 400?0.948:distcut == 300?0.852:distcut == 200?0.565:distcut==159?0.376:100000);

  const double bglow = 7.13*3, bghigh = 100, siglow = 0.3, sighigh = 0.801*2;
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
    
  teff  = eff * pow(0.64, ncut);
  gceff = eff * pow(0.93, ncut);

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

//#include "ehistbg.C"
//#include "ehistsig.C"

  double nfrac[nrbins];
  double nfracerr[nrbins];

  TH2D tmp("tmp", "", 2, 0, 2, nrbins, 0, nrbins);
  t->Draw(Form("classi(mx, my, mz):n>=%d && n <= 3 >> tmp", ncut), "ndecay == 0 && miche < 12 && !earlymich");
  for(int i = 0; i < nrbins; i++){
    double muonswithn = tmp.GetBinContent(2, i+1);
    double allmuons   = tmp.GetBinContent(1, i+1) + muonswithn;

    nfrac[i] = muonswithn/allmuons;
    nfracerr[i] = sqrt(nfrac[i]*(1-nfrac[i]))/sqrt(allmuons);
    printf("frac with neutrons in %d is %f\n", i, nfrac[i]);
  }

  if(!ehistbg){ 
    ehistbg = new TH2D("ehistbg", "", 5, 0, 5, 120, 0, 15);
    fprintf(stderr, "regenerating ehistbg\n");
   
    ehistbg->Sumw2();
    t->Draw("e:classi(mx, my, mz) >> ehistbg ", Form("%s && dt > %f", cutnoncut, bglow*1e3));

    ehistbg->Scale((sighigh-siglow)/(bghigh-bglow));

    for(int j = 1; j <= nrbins; j++)
      for(int i = 1; i <= ehistbg->GetNbinsY(); i++){
        ehistbg->SetBinContent(j, i, ehistbg->GetBinContent(j, i)*nfrac[j-1]); 
        ehistbg->SetBinError  (j, i, ehistbg->GetBinError  (j, i)*nfrac[j-1]); 
      }
  }

  for(int j = 0; j < nrbins; j++){
    const double floorlead = (bghigh-bglow)/(sighigh-siglow)*ehistbg->Integral(j+1, j+1, 1, ehistbg->GetNbinsY());
    bgerr[j] = sqrt(1/floorlead + nfracerr[j]*nfracerr[j]);
    printf("background fractional error for %d: %f\n", j, bgerr[j]);
  }

  if(!ehistsig){ 
    ehistsig = new TH2D("ehistsig", "", 5, 0, 5, 120, 0, 15);
    fprintf(stderr, "regenerating ehistsig\n");
    t->Draw("e:classi(mx, my, mz) >> ehistsig", Form("%s && dt > %f && dt < %f", cut, siglow*1e3, sighigh*1e3));
  }
  
  TF1 betahe6("betahe6", "(x+0.511)*sqrt((x+0.511)**2 - 0.511**2) * (3.5076 - x)**2", 0, 3.51);
  he6spec = new TH1D("he6spec", "", 120, 0, 15);
  li8spec = new TH1D("li8spec", "", 120, 0, 15);
  for(int j = 0; j < 5; j++)
    n16spec[j] = new TH1D(Form("n16spec%d", j), "", 120, 0, 15);
  for(int i = 0; i < 1e8; i++){
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
  for(int j = 0; j < nrbins; j++)
    n16mc.Draw(Form("ctEvisID >> n16spec%d", j),
               Form("classi(bamacorrxy(ctX[0], ctEvisID),"
                           "bamacorrxy(ctX[1], ctEvisID),"
                           "bamacorrz (ctX[2], ctEvisID)) == %d", j));

  he6spec->Scale(ehistsig->Integral()/he6spec->Integral());
  li8spec->Scale(ehistsig->Integral()/li8spec->Integral());
  for(int j = 0; j < nrbins; j++)
    n16spec[j]->Scale(ehistsig->Integral()/n16spec[j]->Integral());

  for(int i = 0; i < 120; i++){
    ali8[i] = li8spec->GetBinContent(i+1);
    ahe6[i] = he6spec->GetBinContent(i+1);
    for(int j = 0; j < nrbins; j++){
      abg[j][i]  =ehistbg ->GetBinContent(j+1, i+1);
      asig[j][i] =ehistsig->GetBinContent(j+1, i+1);
      an16[j][i] = n16spec[j]->GetBinContent(i+1);
    }
  }


  TMinuit * mn = new TMinuit(11);
  mn->SetFCN(fcn);
  int err;
  mn->mnparm(0, "bgtargin",   1, 0.01,  0,  0, err);
  mn->mnparm(1, "bgtargmid1", 1, 0.01,  0,  0, err);
  mn->mnparm(2, "bgtargmid2", 1, 0.01,  0,  0, err);
  mn->mnparm(3, "bgtargout",  1, 0.01,  0,  0, err);
  mn->mnparm(4, "bggc",       1, 0.01,  0,  0, err);
  mn->mnparm(5, "li8",        1, 0.001, 0, 1, err);
  mn->mnparm(6, "he6",        1, 0.01,  0, 1*(ncut == 3? 100:1), err);
  mn->mnparm(7, "n16targin",  1, 0.001, 0, 0.1, err);
  mn->mnparm(8, "n16targmid", 1, 0.001, 0, 0.1, err);
  mn->mnparm(9, "n16targout", 1, 0.001, 0, 0.1, err);
  mn->mnparm(10, "n16gc",     1, 0.001, 0, 0.1, err);
  printf("err? %d\n", err);

  if(ncut > 2){ // no n16 
    for(int i = nrbins+2+1; i < nrbins+2+1+4; i++){
      mn->Command(Form("SET PAR %d 0", i));
      mn->Command(Form("FIX %d", i));
    }
  }

  mn->Command(Form("SET PAR %d 0", nrbins+1+1));
  mn->Command(Form("FIX %d", nrbins+1+1));
  mn->Command("MIGRAD");

  const double chi2nohe6 = mn->fAmin;

  mn->Command(Form("RELEASE %d", nrbins+1+1));
  mn->Command(Form("SET PAR %d 0.05", nrbins+1+1));
  mn->Command(Form("MIGRAD"));
  mn->Command(Form("MINOS 10000 %d", nrbins+1+1));

  const double chi2he6 = mn->fAmin;

  if(chi2nohe6 - chi2he6 > 0)
    printf("%sSignificance of He-6: %f%s\n",
           RED, sqrt(chi2nohe6 - chi2he6), CLR);

  const double upforlim = 2.71 + chi2nohe6 - chi2he6;

  double dum, bgnorm[nrbins], li8norm, he6norm, n16norm[nrbins];

  for(int j = 0; j < nrbins; j++) mn->GetParameter(j, bgnorm[j], dum);
  mn->GetParameter(nrbins, li8norm, dum);
  mn->GetParameter(nrbins+1, he6norm, dum);
  for(int j = 0; j < nrbins; j++) mn->GetParameter(j+nrbins+2, n16norm[j], dum);

  double li8normerr, he6normerr;
  mn->mnerrs(nrbins, dum, dum, li8normerr, dum);
  mn->mnerrs(nrbins+1, dum, dum, he6normerr, dum);

  const double captures = 489.509 * 367.;

  const double raw_nhe6 = he6spec->Integral() * he6norm *
                      ((mus[0] + mus[1] + mus[2] + mus[3])*teff + mus[4]*gceff);
  const double raw_signhe6 = he6spec->Integral() * he6normerr *
                      ((mus[0] + mus[1] + mus[2] + mus[3])*teff + mus[4]*gceff);

  const double cooked_nhe6 = he6spec->Integral() * he6norm *
                      (mus[0] + mus[1] + mus[2] + mus[3] + mus[4]);
  const double cooked_signhe6 = he6spec->Integral() * he6normerr *
                      (mus[0] + mus[1] + mus[2] + mus[3] + mus[4]);

  printf("%sefficiency: %f %f%s\n", RED, teff, gceff, CLR);
  printf("%sraw    He-6: %f +- %f%s\n", RED, raw_nhe6, raw_signhe6, CLR);
  printf("%scooked He-6: %f +- %f%s\n", RED, cooked_nhe6, cooked_signhe6, CLR);
  printf("%sHe-6 prob: %f +- %f%s\n", RED, cooked_nhe6/captures, cooked_signhe6/captures, CLR);

  const double defaultfup = mn->fUp;
  mn->fUp = upforlim;
  mn->Command(Form("MINOS 10000 %d", nrbins+1+1));
  mn->fUp = defaultfup;
 
  double he6up90;
  mn->GetParameter(nrbins+1, he6norm, dum);
  mn->mnerrs(nrbins+1, he6up90, dum, dum, dum);

  printf("%s 90%% upper limit: %f+%f = %f%s\n", RED, he6norm, he6up90, he6norm+he6up90, CLR);
  
  mn->Command("set print -3");
  mn->Command(Form("Fix %d", nrbins+1+1));

  const double scan = he6up90 > 0.0001? he6up90: 0.001;

  double sump = 0;
  for(int i =0;i<1000; i++){
    mn->Command(Form("set  par %d %f", nrbins+1+1, i*0.01*scan));
    mn->Command("Migrad");
    const double p = exp(-mn->fAmin+chi2he6);
    printf("%f\t%g\t%g\n", i*0.01*scan, -mn->fAmin+chi2he6, p);
    sump += p;
  }

  double sump2 = 0;
  int i;
  for(i =0;i<1000; i++){
    mn->Command(Form("set par %d %f", nrbins+1+1, i*0.01*scan));
    mn->Command("Migrad");
    const double p = exp(-mn->fAmin+chi2he6);
    sump2 += p;
    printf("%f\t%g\n", i*0.01*scan, sump2/sump);
    if(sump2>sump*0.9){ printf("%s 90% bays lim: %f (best = %f)%s\n", RED, i*0.01*scan, scan, CLR); break;}
  }
  mn->Command("set print 0");

  const double he6norm90 = i*0.01*scan;

  double n90he6 = 0;
  for(int j = 1; j <= nrbins; j++){
    c2->cd(j);

    TH1D * ehistbg_p=ehistbg->ProjectionY(Form("ehistbg_%d", j), j, j);
    TH1D * ehistsig_p=ehistsig->ProjectionY(Form("ehistsig_%d", j), j, j);
    TH1D * li8spec_p = (TH1D * )li8spec->Clone(Form("li8spec_%d", j));
    TH1D * he6spec_p = (TH1D * )he6spec->Clone(Form("he6spec_%d", j));
    TH1D * he6demo   = (TH1D * )he6spec->Clone(Form("he6demo_%d", j));
    TH1D * n16spec_p = (TH1D * )n16spec[j-1]->Clone(Form("n16specp_%d", j));

    const double thiseff = j < nrbins?teff:gceff;

    ehistbg_p->Scale(bgnorm[j-1]);
    n16spec_p->Scale(n16norm[j-1]);
    li8spec_p->Scale(li8norm*mus[j-1]);
    he6spec_p->Scale(he6norm*mus[j-1]*thiseff);
    he6demo->Scale(he6norm90*mus[j-1]*thiseff);

    ehistsig_p->SetLineColor(kRed);
    ehistsig_p->SetMarkerColor(kRed);
    ehistsig_p->GetXaxis()->SetRange(6, 120);
    ehistsig_p->Draw("e");

    ehistbg_p->Draw("samehist");

    if(ncut < 3){
      n16spec_p->SetLineColor(kOrange);
      n16spec_p->Draw("samehist");
    }

    TH1D * bg_plus_li8 = (TH1D*)ehistbg_p->Clone(Form("bg_plus_li8_%d", j));
    bg_plus_li8->Add(li8spec_p);
    bg_plus_li8->SetLineColor(kBlue);
    bg_plus_li8->Draw("samehist");

    TH1D * bg_plus_li8_n16 = (TH1D*)bg_plus_li8->Clone(Form("bg_plus_li8_n16_%d", j));
    bg_plus_li8_n16->Add(n16spec_p);
    bg_plus_li8_n16->SetLineColor(kViolet);
    if(ncut < 3) bg_plus_li8_n16->Draw("samehist");

    TH1D * bg_plus_li8_n16_he6 = (TH1D*)bg_plus_li8_n16->Clone(Form("bg_plus_li8_n16_he6_%d", j));
    bg_plus_li8_n16_he6->Add(he6spec_p);
    bg_plus_li8_n16_he6->SetLineColor(kGreen+2);
    bg_plus_li8_n16_he6->Draw("samehist");

    he6spec_p->SetLineWidth(1);
    he6spec_p->Draw("samehist");

    he6demo->SetLineStyle(kDashed);
    he6demo->Draw("samehist");

    TH1D * alldemo = new TH1D(Form("alldemo%d", j), "", 120, 0, 15);
    alldemo->Add(ehistbg_p);
    alldemo->Add(li8spec_p);
    alldemo->Add(n16spec_p);
    alldemo->Add(he6demo);
    alldemo->SetLineColor(kGreen+2);
    alldemo->SetLineStyle(kDashed);
    alldemo->Draw("samehist");

    n90he6 += he6demo->Integral(1, he6demo->GetNbinsX())/thiseff;

    /*alldemo->Rebin(4);
    he6demo->Rebin(4);
    bg_plus_li8_plus_he6->Rebin(4);
    bg_plus_li8->Rebin(4);
    ehistbg_p->Rebin(4);
    ehistsig_p->Rebin(4); */
  }

  printf("%s 90%% upper limit n He-6, eff corrected: %f%s\n", RED, n90he6, CLR);
  printf("%s 90%% upper limit He-6 prob: %f%s\n", RED, n90he6/captures, CLR);

  c2->SaveAs(Form("he6-ncut%d.pdf", ncut));
  c2->SaveAs(Form("he6-ncut%d.C", ncut));

  ehistsig->SaveAs(Form("tmp-ehistsig-ncut%d.C", ncut));
  ehistbg->SaveAs(Form("tmp-ehistbg-ncut%d.C", ncut));
}

