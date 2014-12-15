#include "TMarker.h"
#include "TH1.h"
#include "TROOT.h"
#include "TGraph.h"
#include "TF1.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TMinuit.h"
#include <string>

static const int nbins = 2000;
static double lowtime = 0.001;
static double hightime = 10.001;
static double binwidth = (hightime - lowtime)/nbins;
static double sig[nbins] = {0};

static double nb12life = 0.02020/log(2.), // +- 0.00002 (0.02ms) (0.1%)
              nb13life = 0.01733/log(2.), // +- 0.00017 (0.17ms)
              li8life = 0.8399/log(2.), // +- 0.0009 (0.9ms) (0.1%)
              c15life = 2.449 /log(2.), // +- 0.005 (0.2%)
              n16life = 7.13  /log(2.), // +- 0.02 (0.3%) 
              be11life=13.81  /log(2.); // +- 0.08 (0.6%)

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
  * 0.565 * 0.82 // delta r for 200mm, fudged down for 5.X MeV gammas
  * 0.9709 // 100s from end of run
  * 0.789 // energy, estimated from scaled n16 MC
  * 0.986 // from ttlastvalid cut, very naive
  * 0.96 // from ttlastmuon cut, vary naive
  * (1-0.00504 * 20.20/178.3) // from b12like cut, using li-9 as a guide
;

void fcn(int & npar, double * gin, double & like, double *par, int flag)
{
  const double acc = par[0], b12 = par[1], li8 = par[2],
               c15 = par[3], n16 = par[4], be11 = par[5], b13 = par[6],
               b12life = par[7], b13life = par[8];

  like = 0;
  for(int i = 0; i < nbins; i++){
    const double tl = lowtime + binwidth*i;
    const double th = lowtime + binwidth*(i+1);

    const double data = sig[i];
    const double model = acc*binwidth
        + b13* b13life*(exp(-tl/b13life) - exp(-th/b13life))
        + b12* b12life*(exp(-tl/b12life) - exp(-th/b12life))
        + li8* li8life*(exp(-tl/li8life) - exp(-th/li8life))
        + c15* c15life*(exp(-tl/c15life) - exp(-th/c15life))
        + n16* n16life*(exp(-tl/n16life) - exp(-th/n16life))
        +be11*be11life*(exp(-tl/be11life)- exp(-th/be11life));

    like += model - data;
    if(data > 0 && model > 0) like += data*log(data/model);
  }

  like *= 2;

  like += pow((b12life - 20.20e-3/log(2))/0.02e-3*log(2),2);
  like += pow((b13life - 17.33e-3/log(2))/0.17e-3*log(2),2);

}

static void scalemarker(TMarker * m)
{
  const double newx = m->GetX() * c15life/c15eff;
  const double newy = m->GetY() * be11life/be11eff;
  m->SetX(newx);
  m->SetY(newy);
}

static void scalegraph(TGraph * g)
{
  for(int i = 0; i < g->GetN(); i++){
    const double newx = g->GetX()[i] * c15life/c15eff;
    const double newy = g->GetY()[i] * be11life/be11eff;
    g->SetPoint(i, newx, newy);
  }
}

void b13finalfit()
{
  const bool fit = true;

  TFile * fiel= 
    new TFile("/cp/s4/strait/fullfido-300s-3-25MeV-20141117.root","read");
  TTree * t = (TTree *) fiel->Get("t");

  TCanvas * c = new TCanvas("c1", "c1", 1000, 1000);
  c->SetLogy();

  const char * const RED = "\033[31;1m"; // bold red
  const char * const CLR = "\033[m"    ; // clear

  const char * const scutwoe = "!earlymich && latennear == 0 && miche < 12 "
    "&& dist < 400 && ttlastvalid > 0.1 && ttlastmuon > 1 && "
    "timeleft > 10.001e3 ";

  const char * const cute = "e > 4 && e < 15";

  const char * const cutb12 = "b12like < 0.4";

  const string scut = Form("%s && %s && %s", scutwoe, cutb12, cute);

  printf("Drawing...\n");
  t->Draw(Form("dt/1000 >> hfit(%d, %f, %f)", nbins, lowtime, hightime),
          scut.c_str());
  TH1D * hfit = (TH1D*)gROOT->FindObject("hfit");

  printf("Drawing...\n");
  t->Draw(Form("dt >> hdisp(%d, %f, %f)", 100, 1., 501.),
          scut.c_str(), "e");
  TH1D * hdisp = (TH1D*)gROOT->FindObject("hdisp");

  for(int i = 0; i < nbins; i++) sig[i] = hfit->GetBinContent(i+1);

  TMinuit * mn = new TMinuit(9);
  mn->SetFCN(fcn);
  int err;
  mn->mnparm(1 -1, "acc",   200, 0.01, 0, 1000, err);
  mn->mnparm(2 -1, "b12",    8e5, 1, 0, 1e7, err);
  mn->mnparm(3 -1, "li8",   600, 0.1, 0, 10000, err);
  mn->mnparm(4 -1, "c15",   1, 0.01, 0, 1000, err);
  mn->mnparm(5 -1, "n16",   1, 0.01, 0, 100, err);
  mn->mnparm(6 -1, "be11",  1, 0.01, 0, 1000, err);
  mn->mnparm(7 -1, "b13",   4e3, 0.01, 0, 1e5, err);
  mn->mnparm(8 -1, "b12life", 0.02020/log(2), 1e-5, 0, 0, err);
  mn->mnparm(9 -1, "b13life", 0.01733/log(2), 1e-5, 0, 0, err);

  mn->Command("MIGRAD");
  mn->Command("MIGRAD");

  TF1 * ee = new TF1("ee", Form(
    " [0]"
    "+[1]*exp(-x/[7])"
    "+[2]*exp(-x/%f)"
    "+[3]*exp(-x/%f)"
    "+[4]*exp(-x/%f)"
    "+[5]*exp(-x/%f)"
    "+[6]*exp(-x/[8])",
    li8life*1000, c15life*1000, n16life*1000, be11life*1000),
    0, 501);
  ee->SetNpx(200);

  double val[9];

  for(int i = 0; i < 9; i++){
    double errd;
    mn->GetParameter(i, val[i], errd);
    if(i < 7)
      ee->SetParameter(i, val[i]*hdisp->GetBinWidth(1)/1000);
    else
      ee->SetParameter(i, val[i]*1000);
  }

  ee->Draw("same");

  TF1* b12= new TF1("b12", Form("[0]*exp(-x/[1])"), 0, 501);
  TF1* li8= new TF1("li8", Form("[0]*exp(-x/%f)", li8life*1000), 0, 501);
  TF1* c15= new TF1("c15", Form("[0]*exp(-x/%f)", c15life*1000), 0, 501);
  TF1* n16= new TF1("n16", Form("[0]*exp(-x/%f)", n16life*1000), 0, 501);
  TF1* be11=new TF1("be11",Form("[0]*exp(-x/%f)",be11life*1000), 0, 501);
  TF1* b13= new TF1("b13", Form("[0]*exp(-x/[1])"), 0, 501);

  b12->SetParameter(0, val[1]*hdisp->GetBinWidth(1)/1000);
  b12->SetParameter(1, val[7]*1000);
  li8->SetParameter(0, val[2]*hdisp->GetBinWidth(1)/1000);
  c15->SetParameter(0, val[3]*hdisp->GetBinWidth(1)/1000);
  n16->SetParameter(0, val[4]*hdisp->GetBinWidth(1)/1000);
  be11->SetParameter(0,val[5]*hdisp->GetBinWidth(1)/1000);
  b13->SetParameter(0, val[6]*hdisp->GetBinWidth(1)/1000);
  b13->SetParameter(1, val[8]*1000);

  b12->SetNpx(400);
  n16->SetNpx(400);
  be11->SetNpx(400);
  li8->SetNpx(400);
  b13->SetNpx(400);

  li8->Draw("same");
//  c15->Draw("same");
//  n16->Draw("same");
//  be11->Draw("same");
  b13->Draw("same");
  b12->Draw("same");

  new TCanvas;

  mn->fUp = 2.3; // 90% in 1D
  mn->Command("set print 0");
  mn->Command("mncont 7 9");

  TGraph * ninty_1d =
    mn->GetPlot()?(TGraph*)((TGraph*)mn->GetPlot())->Clone():NULL;

  if(ninty_1d) ninty_1d->Draw("alf");
  
}
