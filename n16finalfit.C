#include "TH1.h"
#include "TROOT.h"
#include "TGraph.h"
#include "TF1.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TMinuit.h"
#include <string>

static const int nbins = 10000;
static double lowtime = 0.001;
static double hightime = 100;
static double binwidth = (hightime - lowtime)/nbins;
static double sig[nbins] = {0};

static double b12life = 0.0202/log(2.),
              li8life = 0.8399/log(2.),
              c15life = 2.449 /log(2.),
              n16life = 7.13  /log(2.),
              be11life=13.81  /log(2.);

static const double n16eff = 1
  * 0.981 // subsequent muons
  * 0.977 // previous muons
  * 0.565 * 0.8 // delta r for 200mm, fudged down for 6.1MeV gammas
  * 0.9709 // 100s from end of run
  * 0.798 // energy, from MC, a bit rough
  * 0.986 // from ttlastvalid cut, very naive
  * 0.96 // from ttlastmuon cut, vary naive
  * (1-0.00504 * 20.20/178.3) // from b12like cut, using li-9 as a guide
;

static const double be11eff = 1
  * 0.981 // subsequent muons
  * 0.977 // previous muons
  * 0.565 // delta r for 200mm
  * 0.9709 // 100s from end of run
  * 0.705 // energy, estimated from scaled b12 mc
  * 0.986 // from ttlastvalid cut, very naive
  * 0.96 // from ttlastmuon cut, vary naive
  * (1-0.00504 * 20.20/178.3) // from b12like cut, using li-9 as a guide
;

static const double c15eff = 1
  * 0.981 // subsequent muons
  * 0.977 // previous muons
  * 0.565 // delta r for 200mm
  * 0.9709 // 100s from end of run
  * 0.789 // energy, estimated from scaled n16 MC
  * 0.986 // from ttlastvalid cut, very naive
  * 0.96 // from ttlastmuon cut, vary naive
  * (1-0.00504 * 20.20/178.3) // from b12like cut, using li-9 as a guide
;

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

//  like += pow((n16 - n16eff*410./n16life)/(n16eff*100./n16life),2);

}

static void scalegraph(TGraph * g)
{
  for(int i = 0; i < g->GetN(); i++){
    const double newx = g->GetX()[i] * c15life/c15eff;
    const double newy = g->GetY()[i] * be11life/be11eff;
    g->SetPoint(i, newx, newy);
  }
}

void n16finalfit()
{
  const bool fit = true;

  TFile * fiel= 
    new TFile("/cp/s4/strait/fullfido-100s-3-25MeV-20141022.root","read");
  TTree * t = (TTree *) fiel->Get("t");

  TCanvas * c = new TCanvas("c1", "c1", 1000, 1000);

  const char * const RED = "\033[31;1m"; // bold red
  const char * const CLR = "\033[m"    ; // clear

  const char * const scutwoe = "!earlymich && nnear == 0 && miche < 12 "
    "&& dist < 200 && ttlastvalid > 0.1 && ttlastmuon > 1 && timeleft > 100e3 ";

  const char * const cute = "e > 4 && e < 10";

  const char * const cutb12 = "b12like < 0.4";

  const string scut = Form("%s && %s && %s", scutwoe, cutb12, cute);

  const char * xscut = string(Form("%s && dt > 4e3 && dt < 20e3 && %s && %s",
                                   scutwoe, cute, cutb12)).c_str();

  printf("Drawing...\n");
  t->Draw(Form("dt/1000 >> hfit(%d, %f, %f)", nbins, lowtime, hightime),
          scut.c_str());
  TH1D * hfit = (TH1D*)gROOT->FindObject("hfit");

  printf("Drawing...\n");
  t->Draw(Form("dt/1000 >> hdisp(%d, %f, %f)", 49, 2., 100.),
          scut.c_str(), "e");
  TH1D * hdisp = (TH1D*)gROOT->FindObject("hdisp");

  for(int i = 0; i < nbins; i++) sig[i] = hfit->GetBinContent(i+1);
  printf("%sN-16 eff: %f\nBe-11 eff: %f\nC-15 eff: %f\n%s",
        RED, n16eff, be11eff, c15eff, CLR);

  TMinuit * mn = new TMinuit(6);
  mn->SetFCN(fcn);
  int err;
  mn->mnparm(1 -1, "acc",   10, 0.01, 0, 100, err);
  mn->mnparm(2 -1, "b12",     4e5, 1, 0, 1e6, err);
  mn->mnparm(3 -1, "li8",   200, 0.1, 0, 1000, err);
  mn->mnparm(4 -1, "c15",   1, 0.01, 0, 1000, err);
  mn->mnparm(5 -1, "n16",   1, 0.01, 0, 100, err);
  mn->mnparm(6 -1, "be11",  1, 0.01, 0, 1000, err);

  mn->Command("MIGRAD");

  TF1 * ee = new TF1("ee", Form(
    " [0]"
    "+[1]*exp(-x/%f)"
    "+[2]*exp(-x/%f)"
    "+[3]*exp(-x/%f)"
    "+[4]*exp(-x/%f)"
    "+[5]*exp(-x/%f)", b12life, li8life, c15life, n16life, be11life),
    0, 100);
  ee->SetNpx(200);

  for(int i = 0; i < 6; i++){
    double val, err;
    mn->GetParameter(i, val, err);
    ee->SetParameter(i, val*hdisp->GetBinWidth(1));
  }

  ee->Draw("same");

  /*
  TF1* acc= new TF1("acc",
    Form("%f*[0]",
    hfit->GetBinWidth(1)), 0, 100);
  TF1* b12= new TF1("b12",
    Form("%f*([0]*exp(-x*log(2)/0.0202) + [1])",
    hfit->GetBinWidth(1)), 0, 100);
  TF1* li8= new TF1("li8",
    Form("%f*([0]*exp(-x*log(2)/0.8399) + [1])",
    hfit->GetBinWidth(1)), 0, 100);
  TF1* n16= new TF1("n16",
    Form("%f*(%f*[0]/7.13*log(2)*exp(-x*log(2)/7.13) + [1])",
    hfit->GetBinWidth(1), n16eff), 0,100);
  TF1* be11=new TF1("be11",
    Form("%f*(%f*[0]/13.81*log(2)*exp(-x*log(2)/13.81) + [1])",
    hfit->GetBinWidth(1), be11eff),0,100);

  b12->SetNpx(400);
  n16->SetNpx(400);
  be11->SetNpx(400);
  li8->SetNpx(400);
*/

  const double Ccaptures = 367*489.509;
  const double Ocaptures = 7.6*489.509;
 
  gMinuit->Command("set print 0");

  gMinuit->fUp = 2.3/2; // 90% in 1D
  gMinuit->Command("mncont 4 6 100");
  TGraph * ninty_1d = (TGraph*)((TGraph*)gMinuit->GetPlot())->Clone();
  scalegraph(ninty_1d);

  double be11lim = 0, c15lim = 0;
  for(int i = 0; i < ninty_1d->GetN(); i++){
    if(ninty_1d->GetY()[i] > be11lim)
      be11lim = ninty_1d->GetY()[i];
    if(ninty_1d->GetX()[i] > c15lim)
      c15lim = ninty_1d->GetX()[i];
  }

  printf("%sBe-11 90%% CL upper limit events: %.3lf %s\n", RED, be11lim, CLR);
  printf("%sBe-11 90%% CL upper limit prob: %.3lg %s\n",
         RED, be11lim/Ccaptures, CLR);

  printf("%sC-15 90%% CL upper limit events: %.3lf %s\n", RED, c15lim, CLR);
  printf("%sC-15 90%% CL upper limit prob: %.3lg %s\n",
         RED, c15lim/Ocaptures, CLR);

  gMinuit->fUp = 4.61/2; // 90% CL contour in 2D
  gMinuit->Command("mncont 4 6 100");
  TGraph * ninty_2d = (TGraph*)((TGraph*)gMinuit->GetPlot())->Clone();
  scalegraph(ninty_2d);

  gMinuit->fUp = 11.83/2; // 99.73% contour in 2D
  gMinuit->Command("mncont 4 6 100");
  TGraph * ninty983_2d = (TGraph*)((TGraph*)gMinuit->GetPlot())->Clone();
  scalegraph(ninty983_2d);

  ninty983_2d->GetXaxis()->SetTitle("N C-15");
  ninty983_2d->GetYaxis()->SetTitle("N Be-11");
  ninty983_2d->Draw("al");
  ninty_2d->Draw("l");
  ninty_1d->Draw("l");
}
