#include "TTree.h"
#include "TROOT.h"
#include "consts.h"

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
  TH1D *li8spec = new TH1D("li8spec","",40,0,20);
  fill_li8(li8spec, 1);

  TH1D * li8bgsubed = (TH1D*) gROOT->FindObject("li8bgsubed");
  if(li8bgsubed)
    li8spec->Scale(li8bgsubed->Integral(9, 40)/li8spec->Integral(9,40));

  li8spec->Draw("samehist");
  li8spec->SetLineColor(kRed);

  for(int i = 1; i < 8; i++){
    li8bgsubed->SetBinContent(i, li8spec->GetBinContent(i));
    li8bgsubed->SetBinError(i, li8spec->GetBinContent(i)/4.);
  }
  printf("MC Efficiency %.2f %%\n", 100*li8spec->Integral(9,40)/li8spec->Integral(1,40));
}

void li8cutefficiency_finalfit()
{
  TFile *_file0 = TFile::Open(rootfile0up, "Read");
  TH1D * offtime = new TH1D("offtime", "", 40, 0, 20);
  offtime->SetLineColor(kGreen+2);
  TH1D * li8e    = new TH1D("li8e",    "", 40, 0, 20);
  TTree * t = (TTree *)_file0->Get("t");
  t->Draw("e >>    li8e", "      dt > 300 && dt < 1100 && dist < 125 && miche < 2 && e > 0 && timeleft > 10.2e3", "hist");
  t->Draw("e >> offtime", "(1/8.)*(abs(dt-7000) < 3200 && dist < 125 && miche < 2 && e > 0 && timeleft > 10.2e3)", "histsame");
  offtime->Sumw2();
  li8e->Sumw2();

  TH1F * li8bgsubed = (TH1F *)li8e->Clone("li8bgsubed");
  li8bgsubed->Sumw2();
  li8bgsubed->Add(offtime, -1);
  mc();
  li8bgsubed->Draw("samee");

  const double cutlow = 4, cuthigh = 15;

  double bb[4] = {0, cutlow, cuthigh, 16};
  TH1D * cuth = (TH1D *)li8bgsubed->Rebin(3, "cut", bb);
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
  printf("Data Efficiency for %.1f-%.1f MeV cut for Li-8: %.2f +- %.2f %%\n",
         cutlow, cuthigh, eff*100, error*100);

};
