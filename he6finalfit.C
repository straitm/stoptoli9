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
#include <stdio.h>
#include "consts.h"
#include "carbondenominators_finalfit_out.h"

// Turn on to use the high-purity muon sample.  I found that this 
// gave similar, but slightly worse, results
//#define HP

double bamacorrz(const double z, const double e)
{
  // Romain's thesis's correction (eq. 7.24):
  return z + 7.466 
       + (0.008475 + 0.01029*e)*z
       - 1.053e-5*z*z
       + 2.694e-8*z*z*z;
}

double bamacorrxy(const double xy, const double e)
{
  return (1.013 - 7.0e-3*e)*xy + 0.0795e-3*xy*fabs(xy);
}

int nreq = 0;

const int nrbins = 5;

double teff[nrbins] = {0};

bool rbinson[nrbins] = { true, true, true, true, true };

double bgerr[nrbins] = {0};
double abg[nrbins][120], asig[nrbins][120], an16[nrbins][120], ali8[nrbins][120], ahe6[nrbins][120];

TH2D * ehistsig = NULL, * ehistbg = NULL;
TH1D * li8spec[5] = {NULL}, * he6spec[5] = {NULL}, * n16spec[5] = {NULL};

TH2D * mdispbg = NULL;

const double disteff[nrbins] = {
0.7444, // +-0.86e-2
0.6645, // +-0.96e-2
0.6280, // +-1.02e-2
0.5932, // +-1.11e-2
0.4606  // +-0.40e-2
};

// The number of stopping muons in each region, from the inside out,
// with an arbitrary overall normalization
const double mus[nrbins] = {
#ifdef HP
2.201618962,
1.960695653,
1.882183636,
2.580928391,
1.139839258
#else
0.36808601,
0.34753130,
0.31973898,
0.27913184,
2.25442232
#endif
};


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
  if(r2 > 1154*1154 || az > 1233 + 0.03*(1154-r)) return 4;
  if(r2 > 1068.*1068. || az > 1068.) return 3;
  if(r2 > 933.*933. || az > 933.) return 2;
  if(r2 > 740.*740. || az > 740) return 1;
  return 0;
}

void he6finalfit(const int nreq_ = 0,
                 const bool r0 = true, const bool r1 = true,
                 const bool r2 = true, const bool r3 = true,
                 const bool r4 = true, const bool jefferys = false)
{
  if(nreq_ < 0){
    printf("negative input means that I will just get compiled\n");
    return;
  }

  rbinson[0] = r0;
  rbinson[1] = r1;
  rbinson[2] = r2;
  rbinson[3] = r3;
  rbinson[4] = r4;

  nreq = nreq_;
  
  TH1::SetDefaultSumw2();

  TCanvas * c2 = new TCanvas();
  c2->Divide(1, nrbins);
  c2->cd(1);
  c2->ToggleEventStatus();

  TFile fiel(rootfile0up);
  TTree * t = (TTree *) fiel.Get("t");

  const double distcut = 200;

#define BASICMUONCUT "miche < 12 && !earlymich && timeleft > 100e3 && "

  char cutnoncut[1000];
  snprintf(cutnoncut, 999, 
#ifdef HP
  "mx**2+my**2 < 1050**2 && mz > -1175 && "
  "abs(fez + 62*ivdedx/2 - 8847.2) < 1000 && rchi2 < 2 && "
#endif
    BASICMUONCUT
    "dist < %f && b12like < 0.4 && ttlastvalid > 0.1 && ttlastmuon>1"
    , distcut);
  char cut[1000];
  snprintf(cut, 999, "latennear == %d && %s", nreq, distcut, cutnoncut);


  const double bglow = 7.13*3, bghigh = 100,
              siglow = 0.3,   sighigh = 0.801*2;

  // delta r cut is by region below
  const double eff = 1
    * light_noise_eff
    * mich_eff
    * sub_muon_eff10 // subsequent muons with 1ms veto
    * 0.9709 // 100s from end of run
    * 0.986 // ttlastvalid
    * 0.994 // b12like
    * (exp(-siglow*log(2)/0.801) - exp(-sighigh*log(2)/0.801))
  ;
    
  // Must cut on late neutrons for these to be valid. Neutron efficiency
  // varies by region within the target due to different muon track
  // lengths.
  //
  // Delta-t neutron efficiencies by region (first index, T0-4/GC)
  // and number of neutrons (second index, 0-3). Because the neutron
  // efficiency is a function of muon energy, the mean efficiency for,
  // e.g., 2 neutrons is not the efficiency of 1 neutron squared, but
  // rather somewhat better than that.
  const double neffdtby_pos_n[5][4] = {
  {1.0, 0.5160,     0.2784, 0.1567},
  {1.0, 0.5655,     0.3468, 0.2273},
  {1.0, 0.6005,     0.3953, 0.2779},
  {1.0, 0.6063,     0.4079, 0.2948},
  {1.0, neff_dt_gc, 0.8075, 0.7297},
  };

  teff[0]=eff*disteff[0]*neffdtby_pos_n[0][nreq]*pow(neff_dr_800_targ,nreq);
  teff[1]=eff*disteff[1]*neffdtby_pos_n[1][nreq]*pow(neff_dr_800_targ,nreq);
  teff[2]=eff*disteff[2]*neffdtby_pos_n[2][nreq]*pow(neff_dr_800_targ,nreq);
  teff[3]=eff*disteff[3]*neffdtby_pos_n[3][nreq]*pow(neff_dr_800_targ,nreq);
  teff[4]=eff*disteff[4]*neffdtby_pos_n[4][nreq]*pow(neff_dr_800_h,   nreq);

  for(int i = 0; i < nrbins; i++)
    printf("%sEfficiency in region %d: %.1f%s\n",
           RED, i, 100*teff[i], CLR);

/*
  const double li8eff = 1
    * light_noise_eff
    * mich_eff
    * sub_muon_eff10 // subsequent muons with 1ms veto
    * 0.9709 // 100s from end of run
    * 0.986 // ttlastvalid
    * 0.9994 // b12like
    * (exp(-siglow*log(2)/0.8399) - exp(-sighigh*log(2)/0.8399))
  ; */

//#include "ehistbg.C"
//#include "ehistsig.C"

  double nfrac[nrbins];
  double nfracerr[nrbins];

  TH2D tmp("tmp", "", 2, 0, 2, nrbins, 0, nrbins);
  t->Draw(Form("classi(dx, dy, dz):latennear==%d >> tmp", nreq),
        // This has to be the same cuts as used for the decay, 
        // except that it can't make any mention of decays
        #ifdef HP
          "mx**2+my**2 < 1050**2 && mz > -1175 && "
          "abs(fez + 62*ivdedx/2 - 8847.2) < 1000 && rchi2 < 2 && "
        #endif
        BASICMUONCUT
        "ndecay == 0");
  for(int i = 0; i < nrbins; i++){
    double muonswithn = tmp.GetBinContent(2, i+1);
    double allmuons   = tmp.GetBinContent(1, i+1) + muonswithn;

    nfrac[i] = muonswithn/allmuons;
    nfracerr[i] = sqrt(nfrac[i]*(1-nfrac[i]))/sqrt(allmuons);
    printf("frac with %d neutrons in %d is %f\n", nreq, i, nfrac[i]);
  }

  if(!ehistbg){ 
    ehistbg = new TH2D("ehistbg", "", 5, 0, 5, 120, 0, 15);
    fprintf(stderr, "regenerating ehistbg\n");
   
    ehistbg->Sumw2();
    t->Draw("e:classi(dx, dy, dz) >> ehistbg ",
            Form("%s && dt > %f", cutnoncut, bglow*1e3));

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
    t->Draw("e:classi(dx, dy, dz) >> ehistsig",
      Form("%s && dt > %f && dt < %f", cut, siglow*1e3, sighigh*1e3));
  }

#ifdef HP
  const double captures = livetime * n_c12cap_hp;
#else
  const double captures = livetime * n_c12cap;
#endif

  if(ehistsig->Integral(4, 32) == 0){
    printf("NO signal events.  Assuming zero background, \n"
           "limit is roughly < %f\n", 2.30258509299404590/captures/teff[1]);
    return;
  }
  else{
    printf("%f signal events\n", ehistsig->Integral(1,32));
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
  mn->mnparm(5, "he6",         1,     0.01,  0, 1*(nreq>2?100:nreq>1?10:10), err);
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
    printf("%sTECHNOTE Table 2: Significance of He-6 with %d neutrons: %f%s\n",
           RED, nreq, sqrt(chi2nohe6 - chi2he6), CLR);
  else
    printf("%sTECHNOTE Table 2: Significance of He-6 with %d neutrons: 0%s\n",
           RED, nreq, CLR);

  const double upforlim = 1.64237441514981608 + chi2nohe6 - chi2he6;

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
    //printf("mn->fNexofi[%d] = %d\n", i, mn->fNexofi[i]);
    if(mn->fNexofi[i] == he6fpn) iforhe = i;
  }

  he6normerrup = mn->fErp[iforhe];
  he6normerrlo = mn->fErn[iforhe];

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
  printf("%s%sHe-6 prob: (%.2f %.2f +%.2f)%%%s\n", RED,
        nreq == 1? "TECHNOTE results.tex probSixHeFromTwelveCnval: ":"",
        100*cooked_nhe6/captures, 100*cooked_signhe6lo/captures, 100*cooked_signhe6up/captures, CLR);

  const double defaultfup = mn->fUp;
  mn->fUp = upforlim;
  printf("Using fUp of %.2f\n", mn->fUp);
  mn->Command(Form("MINOS 10000 %d", he6fpn));
  mn->fUp = defaultfup;
 
  double he6up90 = mn->fErp[he6fpn-1];

  const double he6norm90 = he6norm+he6up90; // XXX?
  printf("%s 90%% upper limit: %f+%f = %f (fit parameter)%s\n",
         "", he6norm, he6up90, he6norm90, "");


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

    n16spec_p->SetLineColor(kOrange);
    n16spec_p->Draw("samehist");

    TH1D * bg_plus_li8 = (TH1D*)ehistbg_p->Clone(Form("bg_plus_li8_%d", j));
    bg_plus_li8->Add(li8spec_p);
    bg_plus_li8->SetLineColor(kBlue);
    bg_plus_li8->Draw("samehist");

    TH1D * bg_plus_li8_n16 = (TH1D*)bg_plus_li8->Clone(Form("bg_plus_li8_n16_%d", j));
    bg_plus_li8_n16->Add(n16spec_p);
    bg_plus_li8_n16->SetLineColor(kViolet);
    bg_plus_li8_n16->Draw("samehist");

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

  c2->SaveAs(Form("he6out/he6-%s-ncut%d-%d%d%d%d%d.pdf", name, nreq,
             rbinson[0], rbinson[1], rbinson[2], rbinson[3], rbinson[4]));
  c2->SaveAs(Form("he6out/he6-%s-ncut%d-%d%d%d%d%d.C", name, nreq,
             rbinson[0], rbinson[1], rbinson[2], rbinson[3], rbinson[4]));

  ehistsig->SaveAs(Form("he6out/tmp-%s-ehistsig-ncut%d-%d%d%d%d%d.C", name,
    nreq, rbinson[0], rbinson[1], rbinson[2], rbinson[3], rbinson[4]));
  ehistbg->SaveAs(Form("he6out/tmp-%s-ehistbg-ncut%d-%d%d%d%d%d.C", name,
    nreq, rbinson[0], rbinson[1], rbinson[2], rbinson[3], rbinson[4]));


  // conversion between the parameter value and production prob 
  const double rat = (n90he6/captures)/he6norm90;

  mn->SetPrintLevel(-1);
  const double scan = 0;
  double sump = 0;
  const int N = 20000;
  double ps[N];

  unsigned int smallcount = 0;
  const double increment = 1e-6;

  for(int i = 1; i < N; i++){
    const double prob = i*increment;
    mn->Command(Form("rel %d\n", he6fpn));
    mn->Command(Form("set par %d %f", he6fpn, prob / rat));
    mn->Command(Form("fix %d\n", he6fpn));
    mn->Command("Migrad");
    const double p = exp(chi2he6-mn->fAmin) * (jefferys? 1/sqrt(prob): 1);
    printf("%8.6f %8.3g %8.3g ", prob, mn->fAmin-chi2he6, p);
    for(int j = 0; j < p*10 - 1; j++) printf("#");
    if     (p*10 - int(p*10) > 0.67) printf("+");
    else if(p*10 - int(p*10) > 0.33) printf("|");
    printf("\n");
    sump += p;
    ps[i] = p;
    if(p < 1e-6 && ++smallcount > 3)
      break;
  }
  
  printf("Norm: %f\n", sump);

  double sump2 = 0;
  for(int i =1;i<N; i++){
    const double prob = i*increment;
    sump2 += ps[i]/sump;
    if(sump2 > 0.9){
       // conservatively use the high side of the interval.
       printf("TECHNOTE results.tex %s: Bayes limit = %.2f%%\n",
              nreq == 0? "probSixHeFromTwelveCZeron":
              nreq == 1? "probSixHeFromTwelveCnlim":
              nreq == 2? "probSixHeFromTwelveCnn":
              nreq == 3? "probSixHeFromTwelveCnnn": "???",
              prob*100);
       break;
    }
  }

  if(smallcount <= 3)
    printf("Not sure you integrated out far enough\n");
}
