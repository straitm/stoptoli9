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

int ncutmin = 0, ncutmax = 100;

const int nrbins = 5;

double teff[nrbins] = {0};

bool rbinson[nrbins] = { true, true, true, true, true };

double bgerr[nrbins] = {0};
double abg[nrbins][120], asig[nrbins][120], an16[nrbins][120], ali8[nrbins][120], ahe6[nrbins][120];

TH2D * ehistsig = NULL, * ehistbg = NULL;
TH1D * li8spec[5] = {NULL}, * he6spec[5] = {NULL}, * n16spec[5] = {NULL};

TH2D * mdispbg = NULL;

const double mus[nrbins] =
  {0.36335419, 0.374389, 0.357561, 0.292535, 2.2830};

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
  double mnbgnorm[nrbins] = { par[0], par[1], par[2],
                              par[3], par[4] };
  double mnhe6norm = par[nrbins];
  double mnli8norm[nrbins] = { par[nrbins+1], par[nrbins+2],
                               par[nrbins+3], par[nrbins+4],
                               par[nrbins+5] };
  double mnn16norm[nrbins] = { par[nrbins+6], par[nrbins+7],
                               par[nrbins+8], par[nrbins+9],
                               par[nrbins+10] };

  like = 0;
  for(int j = 0; j < nrbins; j++){
    if(!rbinson[j]) continue;
    bool datahasstarted = false;
    for(int i = 3 /* 0.375 MeV */; i < 120; i++){
      double model = mnbgnorm[j]*abg[j][i]
                   + mnn16norm[j]*an16[j][i]
                   + mnli8norm[j]*ali8[j][i]
                   + mnhe6norm*mus[j]*ahe6[j][i]*teff[j];
      if(model < 0) model = 0;
      double data = asig[j][i];
      if(datahasstarted || data > 0) datahasstarted = true;
      else continue;
      like += model - data;
      if(data > 0 && model > 0) like += data*log(data/model);
    } 
    like += pow((1-mnbgnorm[j])/bgerr[j], 2);
  }

  like *= 2;
}

int classi(const double x, const double y, const double z)
{
  const double r2 = x*x+y*y;
  const double r = sqrt(r2);
  const double az = abs(z);
  if(r2 > 1150*1150 || az > 1229 + 0.03*(1150-r)) return 4;
  if(r2 > 1045.*1045. || az > 1127.) return 3;
  if(r2 > 913.*913. || az > 985.) return 2;
  if(r2 > 724.*724. || az > 781) return 1;
  return 0;
}

void he6finalfit(const int ncutmin_ = 0, const int ncutmax_ = 100,
                 const bool r0 = true, const bool r1 = true,
                 const bool r2 = true, const bool r3 = true,
                 const bool r4 = true)
{
  if(ncutmax_ < 0 || ncutmin_ < 0){
    printf("negative input means that I will just get compiled\n");
    return;
  }

  rbinson[0] = r0;
  rbinson[1] = r1;
  rbinson[2] = r2;
  rbinson[3] = r3;
  rbinson[4] = r4;

  ncutmin = ncutmin_;
  ncutmax = ncutmax_;
  
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
    "miche < 12 && dist < %f && timeleft > 100e3 && " // XXX miche > 0.8MeV?
    "b12like < 0.4 && !earlymich && ttlastvalid > 0.1 && ttlastmuon>1"
    , distcut);
  char cut[1000];
  snprintf(cut, 999, "latennear >= %d && latennear <= %d && %s", ncutmin, ncutmax, distcut, cutnoncut);

  const double distcuteff = (distcut == 400?0.944:distcut == 300?0.852:distcut == 200?0.565:distcut==159?0.376:100000);

  const double bglow = 7.13*3, bghigh = 100, siglow = 0.3, sighigh = 0.801*2;
  const double eff = 1
    * 0.981 // subsequent muons
    * 0.977 // previous muons
    * distcuteff // delta r
    * 0.9709 // 100s from end of run
    * 0.986 // ttlastvalid
    * 0.96 // ttlastmuon
    * 0.994 // b12like
    * (exp(-siglow*log(2)/0.801) - exp(-sighigh*log(2)/0.801))
  ;
    
  // must cut on late neutrons for these to be valid
  teff[0] = eff * pow(0.5685, ncutmin);
  teff[1] = eff * pow(0.6437, ncutmin);
  teff[2] = eff * pow(0.6634, ncutmin);
  teff[3] = eff * pow(0.6906, ncutmin);
  teff[4] = eff * pow(0.90,   ncutmin);

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
  t->Draw(Form("classi(dx, dy, dz):latennear>=%d && latennear <= %d >> tmp", ncutmin, ncutmax),
          "ndecay == 0 && miche < 12 && !earlymich && timeleft > 100e3");
  for(int i = 0; i < nrbins; i++){
    double muonswithn = tmp.GetBinContent(2, i+1);
    double allmuons   = tmp.GetBinContent(1, i+1) + muonswithn;

    nfrac[i] = muonswithn/allmuons;
    nfracerr[i] = sqrt(nfrac[i]*(1-nfrac[i]))/sqrt(allmuons);
    printf("frac with %d-%d neutrons in %d is %f\n", ncutmin, ncutmax, i, nfrac[i]);
  }

  if(!ehistbg){ 
    ehistbg = new TH2D("ehistbg", "", 5, 0, 5, 120, 0, 15);
    fprintf(stderr, "regenerating ehistbg\n");
   
    ehistbg->Sumw2();
    t->Draw("e:classi(dx, dy, dz) >> ehistbg ", Form("%s && dt > %f", cutnoncut, bglow*1e3));

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
    t->Draw("e:classi(dx, dy, dz) >> ehistsig", Form("%s && dt > %f && dt < %f", cut, siglow*1e3, sighigh*1e3));
  }
  
  for(int j = 0; j < 5; j++){
    n16spec[j] = new TH1D(Form("n16spec%d", j), "", 120, 0, 15);
    li8spec[j] = new TH1D(Form("li8spec%d", j), "", 120, 0, 15);
    he6spec[j] = new TH1D(Form("he6spec%d", j), "", 120, 0, 15);
  }


  TChain he6mc("data");
  he6mc.Add("/cp/s4/strait/he6mc/*.root");
  for(int j = 0; j < nrbins; j++)
    he6mc.Draw(Form("ctEvisID >> he6spec%d", j),
               Form("classi(bamacorrxy(ctX[0], ctEvisID),"
                           "bamacorrxy(ctX[1], ctEvisID),"
                           "bamacorrz (ctX[2], ctEvisID)) == %d", j));

  TChain li8mc("data");
  li8mc.Add("/cp/s4/strait/li8mc/*.root");
  for(int j = 0; j < nrbins; j++)
    li8mc.Draw(Form("ctEvisID >> li8spec%d", j),
               Form("classi(bamacorrxy(ctX[0], ctEvisID),"
                           "bamacorrxy(ctX[1], ctEvisID),"
                           "bamacorrz (ctX[2], ctEvisID)) == %d", j));

  TChain n16mc("data");
  n16mc.Add("/cp/s4/strait/n16mc/*.root");
  for(int j = 0; j < nrbins; j++)
    n16mc.Draw(Form("ctEvisID >> n16spec%d", j),
               Form("classi(bamacorrxy(ctX[0], ctEvisID),"
                           "bamacorrxy(ctX[1], ctEvisID),"
                           "bamacorrz (ctX[2], ctEvisID)) == %d", j));

  for(int j = 0; j < nrbins; j++){
    li8spec[j]->Scale(ehistsig->Integral()/li8spec[j]->Integral());
    n16spec[j]->Scale(ehistsig->Integral()/n16spec[j]->Integral());
    he6spec[j]->Scale(ehistsig->Integral()/he6spec[j]->Integral());
  }

  for(int i = 0; i < 120; i++){
    for(int j = 0; j < nrbins; j++){
      abg[j][i]  =ehistbg ->GetBinContent(j+1, i+1);
      asig[j][i] =ehistsig->GetBinContent(j+1, i+1);
      an16[j][i] = n16spec[j]->GetBinContent(i+1);
      ali8[j][i] = li8spec[j]->GetBinContent(i+1);
      ahe6[j][i] = he6spec[j]->GetBinContent(i+1);
    }
  }


  TMinuit * mn = new TMinuit(11);
  mn->SetFCN(fcn);
  int err;

  const int li8fpn = nrbins+1+1, n16fpn = 2*nrbins+1+1,
            he6fpn = nrbins+1;

  mn->mnparm(0, "bgtargin",    1,     0.03,  0, 0, err);
  mn->mnparm(1, "bgtargmid1",  1,     0.03,  0, 0, err);
  mn->mnparm(2, "bgtargmid2",  1,     0.03,  0, 0, err);
  mn->mnparm(3, "bgtargout",   1,     0.03,  0, 0, err);
  mn->mnparm(4, "bggc",        1,     0.03,  0, 0, err);
  mn->mnparm(5, "he6",         1,     0.01,  0, 1*(ncutmin>2?100:ncutmin>1?10:10), err);
  mn->mnparm(6, "li8targin",   0.003, 0.001, 0, 0.1, err);
  mn->mnparm(7, "li8targmid1", 0.003, 0.001, 0, 0.1, err);
  mn->mnparm(8, "li8targmid2", 0.003, 0.001, 0, 0.1, err);
  mn->mnparm(9, "li8targout",  0.003, 0.001, 0, 0.1, err);
  mn->mnparm(10, "li8gc",      0.003, 0.001, 0, 0.1, err);
  mn->mnparm(11, "n16targin",  0.001, 0.001, 0, 0.1, err);
  mn->mnparm(12, "n16targmid1",0.001, 0.001, 0, 0.1, err);
  mn->mnparm(13, "n16targmid2",0.001, 0.001, 0, 0.1, err);
  mn->mnparm(14,"n16targout",  0.001, 0.001, 0, 0.1, err);
  mn->mnparm(15,"n16gc",           1, 0.001, 0,  1,  err);
  printf("err? %d\n", err);

  // Don't float accidentals, Li8 or N16 in regions that we ignore.
  for(int i = 0; i < nrbins; i++){
    if(!rbinson[i]){
      mn->Command(Form("FIX %d", i+1));
      mn->Command(Form("FIX %d", n16fpn+i));
      mn->Command(Form("FIX %d", li8fpn+i));
    }
  }

  if(ncutmin > 2){ // no n16, no li8
    for(int i = nrbins+2+1; i < nrbins+2+1+4; i++){
      mn->Command(Form("SET PAR %d 0", li8fpn+i));
      mn->Command(Form("FIX %d", li8fpn+i));
      mn->Command(Form("SET PAR %d 0", n16fpn+i));
      mn->Command(Form("FIX %d", n16fpn+i));
    }
  }

  mn->Command(Form("SET PAR %d 0", he6fpn));
  mn->Command(Form("FIX %d", he6fpn));
  mn->Command("MIGRAD");

  const double chi2nohe6 = mn->fAmin;

  mn->Command(Form("RELEASE %d", he6fpn));
  mn->Command(Form("SET PAR %d 0.05", he6fpn));
  mn->Command("MIGRAD");
  mn->Command(Form("MINOS 10000 %d", he6fpn));

  const double chi2he6 = mn->fAmin;

  if(chi2nohe6 - chi2he6 > 0)
    printf("%sSignificance of He-6: %f%s\n",
           RED, sqrt(chi2nohe6 - chi2he6), CLR);

  const double upforlim = 2.71 + chi2nohe6 - chi2he6;

  double dum, bgnorm[nrbins], li8norm[nrbins], he6norm, n16norm[nrbins];

  mn->GetParameter(he6fpn-1, he6norm, dum);
  for(int j = 0; j < nrbins; j++){
    mn->GetParameter(j+n16fpn-1, n16norm[j], dum);
    mn->GetParameter(j+li8fpn-1, li8norm[j], dum);
    mn->GetParameter(j, bgnorm[j], dum);
  }

  double he6normerrup, he6normerrlo;

  int iforhe = 0;
  for(int i = 0; i < 11; i++){
    printf("mn->fNexofi[%d] = %d\n", i, mn->fNexofi[i]);
    if(mn->fNexofi[i] == he6fpn) iforhe = i;
  }

  he6normerrup = mn->fErp[iforhe];
  he6normerrlo = mn->fErn[iforhe];

  const double captures = 489.509 * 367.;

  double raw_nhe6 = 0, raw_signhe6up = 0, raw_signhe6lo = 0;
  double cooked_nhe6 = 0, cooked_signhe6up = 0, cooked_signhe6lo = 0;
  for(int i = 0; i < nrbins; i++){
    raw_nhe6 += he6spec[i]->Integral() * he6norm * mus[i] * teff[i];
    raw_signhe6up+=he6spec[i]->Integral()*he6normerrup*mus[i]* teff[i];
    raw_signhe6lo+=he6spec[i]->Integral()*he6normerrlo*mus[i]* teff[i];

    cooked_nhe6    += he6spec[i]->Integral() * he6norm  * mus[i];
    cooked_signhe6up+=he6spec[i]->Integral()*he6normerrup*mus[i];
    cooked_signhe6lo+=he6spec[i]->Integral()*he6normerrlo*mus[i];
  }

  printf("%sraw    He-6: %f %f +%f%s\n", RED, raw_nhe6, raw_signhe6lo, raw_signhe6up, CLR);
  printf("%scooked He-6: %f %f +%f%s\n", RED, cooked_nhe6, cooked_signhe6lo, cooked_signhe6up, CLR);
  printf("%sHe-6 prob: %f %f +%f%s\n", RED, cooked_nhe6/captures, cooked_signhe6lo/captures, cooked_signhe6up/captures, CLR);

  const double defaultfup = mn->fUp;
  mn->fUp = upforlim;
  //mn->Command(Form("MINOS 10000 %d", he6fpn));
  mn->fUp = defaultfup;
 
  double he6up90 = mn->fErp[he6fpn-1];

  printf("%s 90%% upper limit: %f+%f = %f%s\n", RED, he6norm, he6up90, he6norm+he6up90, CLR);

  mn->Command(Form("Fix %d", he6fpn));

  const double scan = he6up90 > 0.0001? he6up90: 0.001;
/*
  double sump = 0;
  for(int i =0;i<1000; i++){
    mn->Command(Form("set  par %d %f", he6fpn, i*0.01*scan));
    mn->Command("Migrad");
    const double p = exp(-mn->fAmin+chi2he6);
    printf("%f\t%g\t%g\n", i*0.01*scan, -mn->fAmin+chi2he6, p);
    sump += p;
  }

  double sump2 = 0;
  int i;
  for(i =0;i<1000; i++){
    mn->Command(Form("set par %d %f", he6fpn, i*0.01*scan));
    mn->Command("Migrad");
    const double p = exp(-mn->fAmin+chi2he6);
    sump2 += p;
    printf("%f\t%g\n", i*0.01*scan, sump2/sump);
    if(sump2>sump*0.9){ printf("%s 90% bays lim: %f (best = %f)%s\n", RED, i*0.01*scan, scan, CLR); break;}
  }
  const double he6norm90 = i*0.01*scan;
*/

  const double he6norm90 = he6up90; // XXX?
  double n90he6 = 0;
  for(int j = 1; j <= nrbins; j++){
    c2->cd(j);

    TH1D * ehistbg_p=ehistbg->ProjectionY(Form("ehistbg_%d", j), j, j);
    TH1D * ehistsig_p=ehistsig->ProjectionY(Form("ehistsig_%d", j), j, j);
    TH1D * li8spec_p = (TH1D * )li8spec[j-1]->Clone(Form("li8specp_%d", j));
    TH1D * he6spec_p = (TH1D * )he6spec[j-1]->Clone(Form("he6specp_%d", j));
    TH1D * he6demo   = (TH1D * )he6spec[j-1]->Clone(Form("he6demo_%d", j));
    TH1D * n16spec_p = (TH1D * )n16spec[j-1]->Clone(Form("n16specp_%d", j));

    const double thiseff = teff[j-1];

    ehistbg_p->Scale(bgnorm[j-1]);
    n16spec_p->Scale(n16norm[j-1]);
    li8spec_p->Scale(li8norm[j-1]);
    he6spec_p->Scale(he6norm*mus[j-1]*thiseff);
    he6demo->Scale(he6norm90*mus[j-1]*thiseff);

    ehistsig_p->SetLineColor(kRed);
    ehistsig_p->SetMarkerColor(kRed);
    ehistsig_p->GetXaxis()->SetRange(6, 120);
    ehistsig_p->Draw("e");

    ehistbg_p->Draw("samehist");

    if(ncutmin < 3){
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
    if(ncutmin < 3) bg_plus_li8_n16->Draw("samehist");

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
  }

  printf("%s 90%% upper limit n He-6, eff corrected: %f%s\n", RED, n90he6, CLR);
  printf("%s 90%% upper limit He-6 prob: %f%s\n", RED, n90he6/captures, CLR);

  const char * const name = "";

  c2->SaveAs(Form("he6-%s-ncut%d-%d-%d%d%d%d%d.pdf", name, ncutmin, ncutmax, rbinson[0], rbinson[1], rbinson[2], rbinson[3], rbinson[4]));
  c2->SaveAs(Form("he6-%s-ncut%d-%d-%d%d%d%d%d.C", name, ncutmin, ncutmax, rbinson[0], rbinson[1], rbinson[2], rbinson[3], rbinson[4]));

  ehistsig->SaveAs(Form("tmp-%s-ehistsig-ncut%d-%d-%d%d%d%d%d.C", name, ncutmin, ncutmax, rbinson[0], rbinson[1], rbinson[2], rbinson[3], rbinson[4]));
  ehistbg->SaveAs(Form("tmp-%s-ehistbg-ncut%d-%d-%d%d%d%d%d.C", name, ncutmin, ncutmax, rbinson[0], rbinson[1], rbinson[2], rbinson[3], rbinson[4]));
}
