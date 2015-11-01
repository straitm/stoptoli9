/* Actually the fitter for C-15 and Be-11, for which N-16 is the major 
 * background */

#include "TMarker.h"
#include "TH1.h"
#include "TROOT.h"
#include "TGraph.h"
#include "TF1.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TMinuit.h"
#include "TGaxis.h"
#include "consts.h"
#include "noncarbondenominators_finalfit.out.h"
#include "carbondenominators_finalfit.out.h"
#include <string>

using std::string;
using std::cout;

// High-purity sample
//#define HP

// Low-purity sample, i.e. the Loose sample minus the HP sample
//#define LP

// Without either of the above #defines, you get the Loose sample.

#if defined(HP)
  const double Ccaptures = n_c12cap_hp*livetime;
  const double Ocaptures = n_o16cap_beta_hp*livetime;
#elif defined(LP)
  const double Ccaptures = (n_c12cap-n_c12cap_hp)*livetime;
  const double Ocaptures = (n_o16cap_beta-n_o16cap_beta_hp)*livetime;
#else
  const double Ccaptures = n_c12cap*livetime;
  const double Ocaptures = n_o16cap_beta*livetime;
#endif
 
static const int nbins = 10000;
static double lowtime = 0.001;
static double hightime = 300;
static double binwidth = (hightime - lowtime)/nbins;
static double sig[nbins] = {0};


static const double n16eff = 1
  * light_noise_eff
  * mich_eff
  * sub_muon_eff10 // subsequent muons, 1ms
  * 0.4836 * 0.8 // delta r for 200mm near acrlyic, fudged down for 6.1MeV gammas
  * 0.9709 // 100s from end of run
  * 0.798 // energy, from MC, a bit rough
  * 0.986 // from ttlastvalid cut, very naive
  * 0.906 // b12
;

static const double be11eff = 1
  * sub_muon_eff10 // subsequent muons, 1ms
  * wholedet_dist200eff // delta r for 200mm
  * 0.9709 // 100s from end of run
  * 0.705 // energy, estimated from scaled b12 mc
  * 0.986 // from ttlastvalid cut, very naive
  * 0.906 // b12
;

static const double c15eff = 1
  * sub_muon_eff10 // subsequent muons, 1ms
  * 0.4836 * 0.82 // delta r for 200mm near acrlyic, fudged down for 5.X MeV gammas
  * 0.9709 // 100s from end of run
  * 0.789 // energy, estimated from scaled n16 MC
  * 0.986 // from ttlastvalid cut, very naive
  * 0.906 // b12
;

static int dopull;

void command(TMinuit * mn, const char * const cmd)
{
  char cmp[7];
  strncpy(cmp, cmd, 6);
  cmp[6] = 0;
  if(!strcmp(cmp, "mncont")) mn->SetPrintLevel(0);
  if(mn->Command(cmd) != 0){
    printf("Trouble with %s\n", cmd);
    mn->SetPrintLevel(1);
    mn->Command(cmd);
    mn->SetPrintLevel(-1);
  }
}

TGraph * getplot(TMinuit * mn)
{
  if(!mn->GetPlot()) return NULL;
  return (TGraph*)((TGraph*)mn->GetPlot())->Clone();
}

void fcn(int & npar, double * gin, double & like, double *par, int flag)
{
  const double acc = par[0], b12 = par[1], li8 = par[2],
               c15 = par[3], n16 = par[4], be11 = par[5];

  like = 0;
  for(int i = 0; i < nbins; i++){
    const double tl = lowtime + binwidth*i;
    const double th = lowtime + binwidth*(i+1);

    const double data = sig[i];
    const double model = acc*binwidth
        + b12* b12life*(exp(-tl/b12life) - exp(-th/b12life))
        + li8* li8life*(exp(-tl/li8life) - exp(-th/li8life))
        + c15* c15life*(exp(-tl/c15life) - exp(-th/c15life))
        + n16* n16life*(exp(-tl/n16life) - exp(-th/n16life))
        +be11*be11life*(exp(-tl/be11life)- exp(-th/be11life));

    like += model - data;
    if(data > 0 && model > 0) like += data*log(data/model);
  }

  like *= 2;

  // pull term for N-16, which is already measured by other people to be
  // 11+-1%. I take this to be 10.7+-1, where I recompute the central
  // value from Measday table 5.13, Kane column, but take Measday's 1%
  // error rather than the ~0.5% that Kane's errors sum to since he
  // probably knows better than me how much to trust Kane's errors.
  like += dopull * pow((n16 - n16eff*Ocaptures*0.107/n16life)/
                             (n16eff*Ocaptures*0.025/n16life), 2);
}

static void scalemarker(TMarker * m)
{
  const double newx = m->GetX() * c15life/c15eff / Ocaptures * 100;
  const double newy = m->GetY() * be11life/be11eff / Ccaptures * 100;
  m->SetX(newx);
  m->SetY(newy);
}

static void scalegraph(TGraph * g)
{
  if(!g) return;
  for(int i = 0; i < g->GetN(); i++){
    const double newx = g->GetX()[i] * c15life/c15eff / Ocaptures * 100;
    const double newy = g->GetY()[i] * be11life/be11eff / Ccaptures * 100;
    g->SetPoint(i, newx, newy);
  }
}

void n16_finalfit()
{
  TFile * fiel= new TFile(rootfile3up,"read");
  TTree * t = (TTree *) fiel->Get("t");

  TCanvas * c1 = new TCanvas("c1", "c1", 0, 0, 160, 160);

  const char * const scutwoe =
#ifdef LP
  "!"
#endif
#if defined(LP) || defined(HP)
  "(mx**2+my**2 < 1050**2 && mz > -1175 && "
  "abs(fez + 62*ivdedx/2 - 8847.2) < 1000 && rchi2 < 2) && "
#endif
    "!earlymich && latennear == 0 && miche < 12 && dist < 200 && "
    "ttlastvalid > 0.1 && ttlastmuon > 1 && timeleft > 300e3 ";

  const char * const cute = "e > 4 && e < 10";

  const char * const cutb12 = "b12like < 0.02";

  const string scut = Form("%s && %s && %s", scutwoe, cutb12, cute);

  printf("Drawing...\n");
  t->Draw(Form("dt/1000 >> hfit(%d, %f, %f)", nbins, lowtime, hightime),
          scut.c_str());
  TH1D * hfit = (TH1D*)gROOT->FindObject("hfit");

  printf("Drawing...\n");
  t->Draw("dt/1000 >> hdisp(50, 3, 203)", scut.c_str(), "e");
  TH1D * hdisp = (TH1D*)gROOT->FindObject("hdisp");

  for(int i = 0; i < nbins; i++) sig[i] = hfit->GetBinContent(i+1);
  printf("%sN-16 eff: %f\nBe-11 eff: %f\nC-15 eff: %f\n%s",
        RED, n16eff, be11eff, c15eff, CLR);

  TMinuit * mn = new TMinuit(6);
  mn->SetPrintLevel(-1);
  mn->SetFCN(fcn);
  int err;
  mn->mnparm(1 -1, "acc",  10, 0.01, 0, 100, err);
  mn->mnparm(2 -1, "b12", 4e5,    1, 0, 1e6, err);
  mn->mnparm(3 -1, "li8", 200,  0.1, 0, 1000, err);
  mn->mnparm(4 -1, "c15",   1, 0.01, 0, 1000, err);
  mn->mnparm(5 -1, "n16",   5, 0.01, 0, 100, err);
  mn->mnparm(6 -1, "be11",  1, 0.01, 0, 1000, err);

  command(mn, "SET PAR 6 0");
  command(mn, "SET PAR 4 0");
  command(mn, "FIX 6");
  command(mn, "FIX 4");

  command(mn, "SET PAR 5 0");
  command(mn, "FIX 5");
  dopull = 0;
  command(mn, "MIGRAD");

  const double nothingchi2 = mn->fAmin;

  command(mn, "SET PAR 5 1");
  command(mn, "REL 5");
  command(mn, "MIGRAD");
  command(mn, "MINOS 10000 5");

  {
    TF1* n16= new TF1("n16i", Form("[0]*exp(-x/%f)", n16life), 0, hightime);
    double val, derr;
    mn->GetParameter(4, val, derr);
    n16->SetParameter(0, val);
    
    const double n16found = n16->Integral(0, 200)/ n16eff;
    const double n16founderrup = n16found / val * mn->fErp[3];
    const double n16founderrlo = n16found / val * mn->fErn[3];
    
    printf("%sN-16 found, no pull term, assuming no Be-11 or "
           "C-15 = %.0f +%.0f %.0f%s\n",
           RED, n16found, n16founderrup, n16founderrlo, CLR);
  }

  TF1 * eep = new TF1("eep", Form(
    " [0]"
    "+[1]*exp(-x/%f)"
    "+[2]*exp(-x/%f)"
    "+[3]*exp(-x/%f)"
    "+[4]*exp(-x/%f)"
    "+[5]*exp(-x/%f)", b12life, li8life, c15life, n16life, be11life),
    0, hightime);
  eep->SetNpx(200);

  double val[6];

  dopull = 1;
  command(mn, "MIGRAD");
  for(int i = 0; i < 6; i++){
    double derr;
    mn->GetParameter(i, val[i], derr);
    eep->SetParameter(i, val[i]*hdisp->GetBinWidth(1));
  }
  //eep->SetParameter(0, 0);
  eep->Draw("same");
  eep->SavePrimitive(cout);

  command(mn, "REL 4");
  command(mn, "REL 6");
  command(mn, "MIGRAD");
  command(mn, "MINOS 10000 5");

  const double allchi2 = mn->fAmin;

  {
    TF1* n16= new TF1("n16i", Form("[0]*exp(-x/%f)", n16life), 0, hightime);
    double val, derr;
    mn->GetParameter(4, val, derr);
    n16->SetParameter(0, val);
    
    const double n16found = n16->Integral(0, 200)/ n16eff;
    const double n16founderrup = n16found / val * mn->fErp[4];
    const double n16founderrlo = n16found / val * mn->fErn[4];
    
    printf("%sN-16 found, no pull term, Be-11 and C-15 free = "
           "%.0f +%.0f %.0f%s\n",
           RED, n16found, n16founderrup, n16founderrlo, CLR);
  }

  printf("%sAll preferred to none, no pull, by %.1f%s\n", RED,
         sqrt(nothingchi2 - allchi2), CLR);

  dopull = 1;
  mn->SetPrintLevel(0);
  mn->Command("MIGRAD");
  command(mn, "MIGRAD");
  command(mn, "MINOS 10000 5");

  TF1 * ee = new TF1("ee", Form(
    " [0]"
    "+[1]*exp(-x/%f)"
    "+[2]*exp(-x/%f)"
    "+[3]*exp(-x/%f)"
    "+[4]*exp(-x/%f)"
    "+[5]*exp(-x/%f)", b12life, li8life, c15life, n16life, be11life),
    0, hightime);
  ee->SetNpx(200);


  for(int i = 0; i < 6; i++){
    double derr;
    mn->GetParameter(i, val[i], derr);
    ee->SetParameter(i, val[i]*hdisp->GetBinWidth(1));
  }


  //for(int i = 1; i <= hdisp->GetNbinsX(); i++)
    //hdisp->SetBinContent(i, hdisp->GetBinContent(i) - ee->GetParameter(0));

  //ee->SetParameter(0, 0);
  ee->SavePrimitive(cout);
  ee->Draw("same");

  TF1* b12= new TF1("b12", Form("[0]*exp(-x/%f)", b12life), 0, hightime);
  TF1* li8= new TF1("li8", Form("[0]*exp(-x/%f)", li8life), 0, hightime);
  TF1* c15= new TF1("c15", Form("[0]*exp(-x/%f)", c15life), 0, hightime);
  TF1* n16= new TF1("n16", Form("[0]*exp(-x/%f)", n16life), 0, hightime);
  TF1* be11=new TF1("be11",Form("[0]*exp(-x/%f)",be11life), 0, hightime);

  b12->SetParameter(0, val[1]*hdisp->GetBinWidth(1));
  li8->SetParameter(0, val[2]*hdisp->GetBinWidth(1));
  c15->SetParameter(0, val[3]*hdisp->GetBinWidth(1));
  n16->SetParameter(0, val[4]*hdisp->GetBinWidth(1));
  be11->SetParameter(0, val[5]*hdisp->GetBinWidth(1));

  b12->SetNpx(400);
  n16->SetNpx(400);
  be11->SetNpx(400);
  li8->SetNpx(400);

  b12->SetLineStyle(kDashed);
  n16->SetLineStyle(kDashed);
  be11->SetLineStyle(kDashed);
  li8->SetLineStyle(kDashed);

  li8->Draw("same");
  //c15->Draw("same");
  n16->Draw("same");
  //be11->Draw("same");
  c1->Modified();
  c1->Update();

  const double n16found = n16->Integral(0, 200)/ n16eff / hdisp->GetBinWidth(1);
  const double n16founderrup = n16found / val[4] * mn->fErp[4];
  const double n16founderrlo = n16found / val[4] * mn->fErn[4];
  
  printf("%sN-16 found with pull = %.0f +%.0f %.0f%s\n",
         RED, n16found, n16founderrup, n16founderrlo, CLR);

  TMarker * best = new TMarker(val[3], val[5], kStar);
  scalemarker(best);

  mn->fUp = 2.296; // 68% contour in 2D
  mn->Command("set strategy 2");
  command(mn, "mncont 4 6 400");
  TGraph * onesigma_2d = getplot(mn); 
  onesigma_2d->SetNameTitle("onesigma", "onesigma");
  scalegraph(onesigma_2d);

  double be11min = 10000, be11max = 0, c15min = 10000, c15max = 0;
  if(onesigma_2d){
    for(int i = 0; i < onesigma_2d->GetN(); i++){
      if(onesigma_2d->GetY()[i] > be11max) be11max=onesigma_2d->GetY()[i];
      if(onesigma_2d->GetX()[i] > c15max)  c15max =onesigma_2d->GetX()[i];

      if(onesigma_2d->GetY()[i] < be11min) be11min=onesigma_2d->GetY()[i];
      if(onesigma_2d->GetX()[i] < c15min)  c15min =onesigma_2d->GetX()[i];
    }
  }

  printf("%sBe-11 68%% CL events with pull: %.3lf-%.3lf %s\n",
         RED, be11min, be11max, CLR);

  printf("%sC-15 68%% CL upper limit events with pull: %.3lf-%.3lf %s\n",
         RED, c15min, c15max, CLR);

  mn->fUp = 2.71; // 90% in 1D
  command(mn, "mncont 4 6 200");
  TGraph * ninty_1d = getplot(mn);
  scalegraph(ninty_1d);

  double be11lim = 0, c15lim = 0;
  if(ninty_1d){
    for(int i = 0; i < ninty_1d->GetN(); i++){
      // Correction for uncertainty on O-16 denominator
      if(ninty_1d->GetX()[i] >  c15lim)  c15lim = lim_inflation_for_obeta*ninty_1d->GetX()[i];
      if(ninty_1d->GetY()[i] > be11lim) be11lim = ninty_1d->GetY()[i];
    }
  }

  printf("%sTECHNOTE results.tex probElevenBeFromTwelveC: Be-11 90%% CL upper limit prob: %.2lf%% %s\n", RED, be11lim, CLR);

  printf("%sTECHNOTE results.tex probFifteenCFromSixteenO: C-15 90%% CL upper limit prob: %.0lf%% %s\n", RED, c15lim, CLR);

  TCanvas * c2 = new TCanvas("c2", "c2", 1000, 1000);
  c2->cd();

  mn->fUp = 4.6051; // 90% CL contour in 2D
  command(mn, "mncont 4 6 200");
  TGraph * ninty_2d = getplot(mn);
  ninty_2d->SetNameTitle("ninty", "ninty");
  scalegraph(ninty_2d);

  if(ninty_2d){
    ninty_2d->GetXaxis()->SetTitle("Probability of ^{15}C (%)");
    ninty_2d->GetYaxis()->SetTitle("Probability of ^{11}Be (%)");
    ninty_2d->GetXaxis()->CenterTitle();
    ninty_2d->GetYaxis()->CenterTitle();
    ((TGaxis*)(ninty_2d->GetXaxis()))->SetMaxDigits(3);
    ((TGaxis*)(ninty_2d->GetYaxis()))->SetMaxDigits(3);
    ninty_2d->Draw("al");
    
    //ninty_1d->Draw("l");
    onesigma_2d->Draw("l");
  }
  best->Draw();

  onesigma_2d->SavePrimitive(cout);
  ninty_2d   ->SavePrimitive(cout);
  printf("best fit marker at %f %f\n", best->GetX(), best->GetY());
}
