#include "TGraphAsymmErrors.h"
#include "TFile.h"
#include "TMinuit.h"
#include "TH1.h"
#include "TF1.h"
#include "TTree.h"
#include "TCanvas.h"
#include <stdio.h>
#include "consts.h"

TGraphAsymmErrors * eh = new TGraphAsymmErrors(),
                 * egd = new TGraphAsymmErrors(),
                  * wh = new TGraphAsymmErrors(),
                 * wgd = new TGraphAsymmErrors(),
                 * gdfrac = new TGraphAsymmErrors();

double min(double a, double b)
{
  return a < b ? a : b;
}

void print()
{
  printf("H: %f %f +%f\n", gMinuit->fU[1], gMinuit->fErn[1], gMinuit->fErp[1]);
  printf("Gd: %f %f +%f\n", gMinuit->fU[4], gMinuit->fErn[4], gMinuit->fErp[4]);
}

void postmuonres_finalfit(const double emin = 60, const double emax = 110)
{
  TFile * f = new TFile(rootfilethru, "read");
  TTree * t = (TTree *)f->Get("t");

  const char * const fstring =
    " [0]/([2]*sqrt(2*TMath::Pi()))*exp(-0.5*((x-[1])/[2])**2)" // Hn
    "+[3]/([5]*sqrt(2*TMath::Pi()))*exp(-0.5*((x-[4])/[5])**2)" // Gdn
    "+[6]" // <-- flat background, non-flat bg --v
    "+[11]*([7]*exp(-0.5*((x-[1]+2.223-0.7010498)/[8])**2)"
    "      +[9]*exp(-0.5*((x-[1]+2.223-2.300148)/[10])**2))";

  TF1 * g =new TF1("g",  fstring, 0, 10);
  TF1 * bg=new TF1("bg", fstring, 0, 10);

  TCanvas * c1 = new TCanvas("c1", "", 0, 0, 500, 400);
             
  eh->SetName("eh");
  egd->SetName("egd");
  wh->SetName("wh");
  wgd->SetName("wgd");

  TH1D * ehist = new TH1D("ehist", "", 100, 1, 11);

  const float halfw = 25;
  {
    const float emean = emin + (emax - emin)/2;
    g->SetParameters(20, 2.223 + 0.001*emean, 0.055*2.223, 10, 8, 0.4, 0);

    g->SetParLimits(0, 1, 100);
    g->SetParLimits(3, 1, 100);

 
    g->FixParameter(7,0.001416994);
    g->FixParameter(8,0.5057105);
    g->FixParameter(9,0.005373248);
    g->FixParameter(10,0.855553);

    g->SetParLimits(5, 0.04, 0.4+0.002*emean);

    g->SetParLimits(11, 0, 10000);
    g->SetParameter(11, 1000);

    g->SetParLimits(6, 0, 10);

    g->SetParLimits(1, 2.2, 2.9);
    g->SetParLimits(2, 0.03*2.223, 0.12*2.223);

    t->Draw("miche >> ehist", Form("ndecay == 0 && michd < 2000 && "
      "latennear > 0 && fq/8300 < %f && fq/8300 > %f && micht < 3000",
      emax, emin), "e");
    ehist->Fit("g", "liq", "");
    gMinuit->Command("MINOS 10000 2");
    gMinuit->Command("MINOS 10000 3");
    gMinuit->Command("MINOS 10000 5");
    gMinuit->Command("MINOS 10000 6");
    gMinuit->Command("showmin");

    for(int j = 1; j <= 1; j++) bg->SetParameter(j, g->GetParameter(j));
    for(int j = 6; j <= 11; j++) bg->SetParameter(j, g->GetParameter(j));
    bg->Draw("same");

    c1->Update(); c1->Modified();
    egd->SetPoint(0, emean, g->GetParameter(4));
    eh ->SetPoint(0, emean, g->GetParameter(1));
    wgd->SetPoint(0, emean, fabs(g->GetParameter(5)/7.95));
    wh ->SetPoint(0, emean, fabs(g->GetParameter(2)/2.223));


    const double xerr = halfw*2/sqrt(12);
    eh ->SetPointError(0, xerr, xerr, -gMinuit->fErn[1], gMinuit->fErp[1]);
    wh ->SetPointError(0, xerr, xerr, min(wh->GetY()[0], fabs(gMinuit->fErn[2])/2.223), gMinuit->fErp[2]/2.223);
    egd->SetPointError(0, xerr, xerr, -gMinuit->fErn[4], gMinuit->fErp[4]);
    wgd->SetPointError(0, xerr, xerr, -gMinuit->fErn[5]/7.95, gMinuit->fErp[5]/7.95);

    print();
  }
}
