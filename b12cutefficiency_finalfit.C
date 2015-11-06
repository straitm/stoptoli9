#include "TTree.h"
#include "TROOT.h"
#include "consts.h"

#include "b12spectrum.C" // <-- C inclusion :-(

/* Spectrum from DOGS MC */
void fill_li8(TH1D * h, const double N)
{
  TH1F *li8dogs = new TH1F("li8dogs","ctEvisID  {ctX[0]**2+ctX[1]**2+ctX[2]**2 < 800**2 && ctEvisID > 0.375}",40,0,20);
   li8dogs->SetBinContent(1,12);
   li8dogs->SetBinContent(2,50);
   li8dogs->SetBinContent(3,127);
   li8dogs->SetBinContent(4,188);
   li8dogs->SetBinContent(5,282);
   li8dogs->SetBinContent(6,393);
   li8dogs->SetBinContent(7,504);
   li8dogs->SetBinContent(8,555);
   li8dogs->SetBinContent(9,624);
   li8dogs->SetBinContent(10,695);
   li8dogs->SetBinContent(11,775);
   li8dogs->SetBinContent(12,838);
   li8dogs->SetBinContent(13,846);
   li8dogs->SetBinContent(14,916);
   li8dogs->SetBinContent(15,917);
   li8dogs->SetBinContent(16,902);
   li8dogs->SetBinContent(17,896);
   li8dogs->SetBinContent(18,789);
   li8dogs->SetBinContent(19,708);
   li8dogs->SetBinContent(20,655);
   li8dogs->SetBinContent(21,547);
   li8dogs->SetBinContent(22,420);
   li8dogs->SetBinContent(23,322);
   li8dogs->SetBinContent(24,214);
   li8dogs->SetBinContent(25,185);
   li8dogs->SetBinContent(26,96);
   li8dogs->SetBinContent(27,42);
   li8dogs->SetBinContent(28,10);
   li8dogs->SetBinContent(29,1);

  li8dogs->Scale(1/li8dogs->Integral());
  h->Add(li8dogs, N);
}

void mc()
{
  TH1D *b12spec = new TH1D("b12spec","",40,0,20);
  fillit(b12spec, nb12b, b12branchp, b12branchbe, b12branchge,b12brancha, b12branch_aw, 1030e4,1.053); // fudge!
  fillit(b12spec, nb13b, b13branchp, b13branchbe, b13branchge,b13brancha, b13branch_aw,   10e4,1.000);
  fill_li8(b12spec, 4e8);

  TH1D * b12bgsubed = (TH1D*) gROOT->FindObject("b12bgsubed");
  if(b12bgsubed)
    b12spec->Scale(b12bgsubed->Integral(2, 40)/b12spec->Integral(2,40));

  b12spec->Draw("samehist");
  b12spec->SetLineColor(kRed);

  b12bgsubed->SetBinContent(1, b12spec->GetBinContent(1));
  b12bgsubed->SetBinError(1, b12spec->GetBinContent(1)/4.);
  printf("MC Efficiency %.2f %%\n", 100*b12spec->Integral(9,40)/b12spec->Integral(1,40));
}

void b12cutefficiency_finalfit(const double cutlow = 4)
{
  norm();
  TFile *_file0 = TFile::Open(rootfile0up, "Read");
  TH1D * offtime = new TH1D("offtime", "", 40, 0, 20);
  offtime->SetLineColor(kGreen+2);
  TH1D * b12e    = new TH1D("b12e",    "", 40, 0, 20);
  TTree * t = (TTree *)_file0->Get("t");
  t->Draw("e >>    b12e", "       dt > 2 && dt <102 && dist < 200 && miche < 5 && e > 0 && timeleft > 10e3", "hist");
  t->Draw("e >> offtime", "0.05*(abs(dt-8000) <1000 && dist < 200 && miche < 5 && e > 0 && timeleft > 10e3)", "histsame");
  offtime->Sumw2();
  b12e->Sumw2();

  TH1F * b12bgsubed = (TH1F *)b12e->Clone("b12bgsubed");
  b12bgsubed->Sumw2();
  b12bgsubed->Add(offtime, -1);
  mc();
  b12bgsubed->Draw("samee");

  const double cuthigh = 15;

  double bb[4] = {0, cutlow, cuthigh, 16};
  TH1D * cuth = (TH1D *)b12bgsubed->Rebin(3, "cut", bb);
  const double eff = cuth->GetBinContent(2)/cuth->Integral(1,3);
  const double error =
    sqrt(pow(cuth->GetBinContent(2)*cuth->GetBinError(2)/(
            pow(cuth->Integral(1,3), 2) ), 2)
        + pow((cuth->GetBinContent(3)+cuth->GetBinContent(1))
          *sqrt( pow(cuth->GetBinError(1), 2) +
            pow(cuth->GetBinError(3), 2) ) /(
            pow(cuth->Integral(1,3), 2)
            ), 2)
        );
  printf("TECHNOTE 4.4.1: Data Efficiency for %.1f-%.1f MeV cut "
         "for B-12+B-13: %.2f +- %.2f %%\n",
         cutlow, cuthigh, eff*100, error*100);
  printf("All the digits: %f +- %f %%\n", eff*100, error*100);

  printf("const double b12energyeff%.0fMeV = %f;\n", cutlow, eff);
  printf("const double b12energyeff%.0fMeV_e = %f;\n", cutlow, error);

  if(cutlow == 4){ // backwards compatibility
    printf("const double b12energyeff = b12energyeff4MeV;\n");
    printf("const double b12energyeff_e = b12energyeff4MeV_e;\n");
  }
}
