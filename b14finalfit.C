#include "consts.h"

void b14finalfit()
{
  TFile fiel(rootfile3up);
  TTree * t = (TTree *) fiel.Get("t");

  const double eff = 1
    * exp(-1.*log(2)/12.5) // 1ms veto and B-14 half-life
    * 0.981 // subsequent muons
    * 0.977 // previous muons
    * wholedet_dist400eff // delta r
    * 0.99709 // 10s from end of run
    * 0.387 // energy from my toy MC, probably a bit conservative
            // since when I try to evaluate B-12 it is low
    * 0.96 // b12like
  ;

  printf("%sEfficiency: %.2f%%%s\n", RED, eff*100, CLR);

  const char * const cut =
    "!earlymich && miche < 12 && b12like < 0.06 && dist < 400 && latennear == 0 && e > 15"
    "&& e < 22 && timeleft > 10e3";

  printf("%sN selected in 10s: %d%s\n", RED, t->GetEntries(Form("%s && dt < 1e4", cut)), CLR);

  const int nin5 = t->GetEntries(Form("%s && dt < 12.5*5", cut));
  printf("%sN selected in 5 half-lives: %d%s\n", RED, nin5, CLR);

  const int nin10 = t->GetEntries(Form("%s && dt < 12.5*10", cut));
  printf("%sN selected in 10 half-lives: %d%s\n", RED, nin10, CLR);

  const double captures = n_o16cap_beta * livetime;

  const double toprob = 1./captures/eff;

  if(nin5 > 0){
    TCanvas c1;

    t->Draw("dt/1000 >> hfit(10000, 0.001, 10)", cut);

    TF1 ee("ee", "[0]*exp(-x*log(2)/0.0202) + " // b12
                 "[1]*exp(-x*log(2)/[2]) + " // b14
                 "[3]*exp(-x*log(2)/0.8399) + " // li8
                 "[4]", 0, 100); // acc

    ee.SetParameters(1, 1, 0.0125, 1, 1);
    ee.FixParameter(2, 0.0125);
    ee.SetParLimits(0, 0, 100);
    ee.SetParLimits(1, 0, 100);
    ee.SetParLimits(3, 0, 100);
    ee.SetParLimits(4, 0, 100);

    hfit->Fit("ee", "le");

    t.Draw("dt/1000 >> hdisp(100, 0.001, 10.001)", cut, "e");

    TF1 * eedisp = ee.Clone("eedisp");
    eedisp->SetNpx(400);
    eedisp->SetLineColor(kRed);

    int tomult[4] = { 0, 1, 3, 4};
    const double mult = hdisp->GetBinWidth(1)/hfit->GetBinWidth(1);
    for(int i = 0; i < 4; i++)
      eedisp.SetParameter(tomult[i], eedisp.GetParameter(tomult[i])*mult);
    eedisp.Draw("same");

    TF1 b12("b12", "[0]*exp(-x*log(2)/0.0202)" , 0, 100);
    TF1 b14("b14", "[0]*exp(-x*log(2)/0.0125)", 0, 100);
    TF1 li8("li8", "[0]*exp(-x*log(2)/0.8399)"  , 0, 100);
    TF1 acc("acc", "[0]", 0, 100);

    b12.SetNpx(400);

    TF1 * parts[4] = { &b12, &b14, &li8, &acc };

    b12.SetParameter(0, eedisp.GetParameter(0));
    b14.SetParameter(0, eedisp.GetParameter(1));
    li8.SetParameter(0, eedisp.GetParameter(3));
    acc.SetParameter(0, eedisp.GetParameter(4));

    for(int i = 0; i < 4; i++){
      parts[i]->SetLineStyle(7);
      parts[i]->SetLineWidth(2);
      parts[i]->Draw("Same");
    }

    const double Nfound = b14->Integral(0, 20)/hdisp->GetBinWidth(1);
    const double Nerrup = Nfound * gMinuit.fErp[1]/ee->GetParameter(1);
    const double Nerrlo = Nfound * gMinuit.fErn[1]/ee->GetParameter(1);

    printf("%sN found: %f +%f %f%s\n", RED, Nfound, Nerrup, Nerrlo, CLR);

    printf("%sEff: %f\nProb: %g +%g %g%s\n",
        RED, eff, toprob*Nfound, toprob*Nerrup, toprob*Nerrlo, CLR);
  }
  else{
    printf("%sNo signal and no background, so limit is <%.2f%%%s\n",
           RED, 2.3026*lim_inflation_for_obeta*toprob*100, CLR);
  }

}
