#include <stdio.h>
#include <string>
using std::string;
#include "TFile.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TH1.h"
#include "TF1.h"
#include "TTree.h"
#include "TMinuit.h"
#include "consts.h"
#include "sub_muon_eff.out.h"
#include "distcuteff_finalfit.out.h"
#include "totallivetime_finalfit.out.h"
#include "carbondenominators_finalfit.out.h"

void b8_finalfit(const int nn = 4)
{
  const double gf = 1-n_c12captarget/n_c12cap, // fraction in GC
               tf =   n_c12captarget/n_c12cap; // fraction in Target


  double neff =
    (nn==4?n4eff_dt_gc: nn==3?n3of4eff_dt_gc: nn==2?n2of4eff_dt_gc: 1e10)
    *gf*neff_dr_800_h
    +
    (nn==4?n4eff_dt_targ_wearly: nn==3?n3of4eff_dt_targ_wearly:
     nn==2?n2of4eff_dt_targ_wearly: 1e10)
    *tf*neff_dr_800_targ;

  const double eff = 1
    * light_noise_eff
    * mich_eff
    * sub_muon_eff05 // subsequent muons
    * wholedet_dist400eff // delta r
    * (livetime_s - num_runs*100.)/livetime_s
    * 0.9403// energy
    * neff
  ;

  printf("neutron efficiency: %.3f\n", neff);
  printf("overall efficiency: %.3f\n",  eff);

  const double captures = n_c12cap * livetime;

  const double toprob = 1./captures/eff;
  TFile * fiel = new TFile(rootfile3up, "read");
  TTree * t = (TTree *) fiel->Get("t");

  TCanvas * c = new TCanvas("c1", "c1");
  c->cd();

  const char * const ndef = "(latennear+ngdnear-latengdnear)";

  const string scut =
  Form("!earlymich && miche < 12 && dist < 400 && %s == 4 && e > 4"
    "&& e < 18 && timeleft > 100e3", ndef);

  const char * const cut = scut.c_str();

  t->Draw("dt/1000 >> hfit(10000, 0.001, 100)", cut);
  TH1 * hfit = (TH1 *)gROOT->FindObject("hfit");

  TF1 * ee = new TF1("ee", "[0]*exp(-x*log(2)/0.0202) + "
               "[1]*exp(-x*log(2)/[2]) + "
               "[3]*exp(-x*log(2)/7.13) + "
               "[4]", 0, 100);

  ee->SetParameters(1, 1, 0.7700, 1, 1);
  ee->FixParameter(3, 0);
  ee->FixParameter(2, 0.7700); // b8

  ee->SetParLimits(0, 0, 10);
  ee->SetParLimits(1, 0, 10);
  ee->SetParLimits(4, 0, 10);
  
  int p0isfixed = 0;
  p0isfixed = 1;
  ee->FixParameter(0, 0);
  ee->FixParameter(4, 0);

  if(hfit->GetEntries()){
    ee->FixParameter(1, 0);
    hfit->Fit("ee", "l");
    ee->ReleaseParameter(1);
    hfit->Fit("ee", "le");
  }

  if(ee->GetParameter(0) < 1e-6){
    ee->FixParameter(0, 0);
    p0isfixed = 1;
    if(hfit->GetEntries()) hfit->Fit("ee", "le");
  }

  if(hfit->GetEntries()){
    t->Draw("dt/1000 >> hdisp(400, 0.001, 20.001)", cut, "hist");
    TH1 * hdisp = (TH1 *)gROOT->FindObject("hdisp");
    if(hdisp->GetBinContent(2) > 5) hdisp->Draw("e");


    TF1 * eedisp = (TF1 *)ee->Clone("eedisp");
    eedisp->SetNpx(400);
    eedisp->SetLineColor(kRed);

    int tomult[4] = { 0, 1, 3, 4};
    const double mult = hdisp->GetBinWidth(1)/hfit->GetBinWidth(1);
    for(int i = 0; i < 4; i++)
      eedisp->SetParameter(tomult[i], eedisp->GetParameter(tomult[i])*mult);
    eedisp->Draw("same");

    TF1 * b12 = new TF1("b12", "[0]*exp(-x*log(2)/0.0202)", 0, 100);
    TF1 * b8  = new TF1("b8", "[0]*exp(-x*log(2)/[1])", 0, 100);
    TF1 * acc = new TF1("acc", "[0]", 0, 100);

    b12->SetNpx(400);
    b8->SetNpx(400);

    TF1 * parts[3] = { b12, b8, acc };

    b12->SetParameter(0, eedisp->GetParameter(0));
    b8 ->SetParameter(0, eedisp->GetParameter(1));
    b8 ->SetParameter(1, eedisp->GetParameter(2));
    acc->SetParameter(0, eedisp->GetParameter(4));

    for(int i = 0; i < 3; i++){
      parts[i]->SetLineStyle(7);
      parts[i]->SetLineWidth(2);
      parts[i]->Draw("Same");
    } 

    const double Nfound = b8->Integral(0, 20)/hdisp->GetBinWidth(1);

    double Nerrup, Nerrlo;
    string errtype;
    if(hfit->GetEntries() < 3){
      errtype = "HESSE";
      Nerrup = Nfound * ee->GetParError(1)/ee->GetParameter(1);
      Nerrlo = Nfound * ee->GetParError(1)/ee->GetParameter(1);
    }
    else{
      errtype = "MINOS";
      Nerrup = Nfound * gMinuit->fErp[1-p0isfixed]/ee->GetParameter(1);
      Nerrlo = Nfound * gMinuit->fErn[1-p0isfixed]/ee->GetParameter(1);
    }

    printf("%sN found: %f +%f %f %s%s\n",
        RED, Nfound, Nerrup, Nerrlo, errtype.c_str(), CLR);


    printf("%sProb: %g +%g %g%s\n", 
        RED, toprob*Nfound, toprob*Nerrup, toprob*Nerrlo, CLR);
  }
  else{
    const double N = 2.30258509299404590;
    printf("%sTECHNOTE results.tex probEightBFromTwelveC: Probability 90%% limit: %.2g%s\n",RED, toprob*N, CLR);
  }
}
