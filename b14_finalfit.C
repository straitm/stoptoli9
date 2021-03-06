#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TH1.h"
#include "TMinuit.h"
#include "consts.h"
#include "sub_muon_eff.out.h"
#include "distcuteff_b12like_finalfit.out.h"
#include "distcuteff_wholeloose_finalfit.out.h"
#include "totallivetime_finalfit.out.h"
#include "noncarbondenominators_finalfit.out.h"

void b14_finalfit()
{
  TFile fiel(rootfile3up);
  TTree * t = (TTree *) fiel.Get("t");

  const double eff = 1
    * exp(-1./b14life) // 1ms veto and B-14 half-life
    * light_noise_eff
    * mich_eff
    * sub_muon_eff05 // subsequent muons
    * wholedet_dist400eff // delta r
    * (livetime_s - num_runs*(10+1))/livetime_s // timeleft + mutime cut
    * 0.387 // energy from my toy MC, probably a bit conservative
            // since when I try to evaluate B-12 it is low
    * b12like006_dist400_eff
  ;

  printf("%sEfficiency: %.2f%%%s\n", RED, eff*100, CLR);

  const char * const cut =
    "!earlymich && miche < 12 && b12like < 0.06 && dist < 400 && latennear == 0 && e > 15"
    "&& e < 22 && timeleft > 10e3 && mutime > 1000";

  printf("%sN selected in 10s: %d%s\n", RED, (int)t->GetEntries(Form("%s && dt < 1e4", cut)), CLR);

  const int nin5 = t->GetEntries(Form("%s && dt < %f*5",
                                      cut, b14life*log(2.)));
  printf("%sN selected in 5 half-lives: %d%s\n", RED, nin5, CLR);

  const int nin10 = t->GetEntries(Form("%s && dt < %f*10",
                                       cut, b14life*log(2.)));
  printf("%sN selected in 10 half-lives: %d%s\n", RED, nin10, CLR);

  const double captures = n_o16cap_beta * livetime;

  const double toprob = 1./captures/eff;

  if(nin5 > 0){
    TCanvas c1;

    TH1D * hfit = new TH1D("hfit", "", 10000, 0.001, 10);

    t->Draw("dt/1000 >> hfit", cut);

    TF1 ee("ee", Form("[0]*exp(-x/%f) + " // b12
                 "[1]*exp(-x*log(2)/[2]) + " // b14
                 "[3]*exp(-x/%f) + " // li8
                 "[4]", b12life, li8life), 0, 100); // acc

    ee.SetParameters(1, 1, b14life/1000, 1, 1);
    ee.FixParameter(2, b14life/1000);
    ee.SetParLimits(0, 0, 100);
    ee.SetParLimits(1, 0, 100);
    ee.SetParLimits(3, 0, 100);
    ee.SetParLimits(4, 0, 100);

    hfit->Fit("ee", "le");

    TH1D * hdisp = new TH1D("hdisp", "", 100, 0.001, 10.001);

    t->Draw("dt/1000 >> hdisp", cut, "e");

    TF1 * eedisp = (TF1 *)ee.Clone("eedisp");
    eedisp->SetNpx(400);
    eedisp->SetLineColor(kRed);

    int tomult[4] = { 0, 1, 3, 4};
    const double mult = hdisp->GetBinWidth(1)/hfit->GetBinWidth(1);
    for(int i = 0; i < 4; i++)
      eedisp->SetParameter(tomult[i], eedisp->GetParameter(tomult[i])*mult);
    eedisp->Draw("same");

    TF1 b12("b12", Form("[0]*exp(-x/%f)", b12life/1e3) , 0, 100);
    TF1 b14("b14", Form("[0]*exp(-x/%f)", b14life/1e3), 0, 100);
    TF1 li8("li8", Form("[0]*exp(-x/%f)", li8life/1e3)  , 0, 100);
    TF1 acc("acc", "[0]", 0, 100);

    b12.SetNpx(400);

    TF1 * parts[4] = { &b12, &b14, &li8, &acc };

    b12.SetParameter(0, eedisp->GetParameter(0));
    b14.SetParameter(0, eedisp->GetParameter(1));
    li8.SetParameter(0, eedisp->GetParameter(3));
    acc.SetParameter(0, eedisp->GetParameter(4));

    for(int i = 0; i < 4; i++){
      parts[i]->SetLineStyle(7);
      parts[i]->SetLineWidth(2);
      parts[i]->Draw("Same");
    }

    const double Nfound = b14.Integral(0, 20)/hdisp->GetBinWidth(1);
    const double Nerrup = Nfound * gMinuit->fErp[1]/ee.GetParameter(1);
    const double Nerrlo = Nfound * gMinuit->fErn[1]/ee.GetParameter(1);

    printf("%sN found: %f +%f %f%s\n", RED, Nfound, Nerrup, Nerrlo, CLR);

    printf("%sEff: %f\nProb: %g +%g %g%s\n",
        RED, eff, toprob*Nfound, toprob*Nerrup, toprob*Nerrlo, CLR);
  }
  else{
    printf("%sTECHNOTE results.tex probFourteenBFromSixteenO: No signal, no bg, limit is <%.2f%%%s\n",
           RED, 2.30258509299404590*lim_inflation_for_obeta*toprob*100, CLR);
  }
}
