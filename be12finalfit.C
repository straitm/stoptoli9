#include "consts.h"

void be12finalfit()
{
  const double be12eff = 
      0.894 // Energy over 3MeV
    * wholedet_dist400eff 
    * (exp(-1*log(2)/21.3) - exp(-100*log(2)/21.3)) // 1-100ms
    * 0.981 // subsequent muons
  ;


  const double b12eff = 
      0.9251 // Energy over 3 MeV
    * wholedet_dist400eff
    * exp(-1*log(2)/20.2)  // After 1ms
    * 0.981 // subsequent muons
   ;

  printf("%sBe-12 eff: %.1f%%\nB-12 eff:  %.1f%%\nTotal eff: %.1f%%%s\n",
         RED, be12eff*100, b12eff*100, b12eff*be12eff*100, CLR);

  TTree * t = new TTree("t", "t");
  t->ReadFile("/cp/s4/strait/be12-20150303.d/header");
  t->ReadFile("/cp/s4/strait/be12-20150303.d/doubles");

  const char * const cut =
    "dt < 100 && miche < 12 && latennear == 0 && "
    "dist < 400 && dist2 < 400 && "
    "e < 12 && e2 < 14 && dt < 250 && dt2 < 250 && dt2-dt > 1";

  t->Draw("dt2-dt >> hfit(249, 1, 250)", cut, "e");
  const int n = t->Scan("dt:e:dt2-dt:e2", cut);
  double upperlimcount;
  if(n == 0){
    upperlimcount = 2.30258509299404590;
  }
  else{
    hfit=hfit;

    TF1 * ee = new TF1("ee", "[0] + [1]*exp(-x*log(2)/20.20)", 0, 500);
    ee->SetNpx(400);

    ee->SetParLimits(1, 0, hfit->GetBinWidth(1));
    ee->SetParLimits(0, 0, hfit->GetBinWidth(1));

    hfit->Fit("ee", "liem", "e");
    gMinuit->fUp = 2.71/2; // 90% in 1D
    gMinuit->Command("MINOS 10000 2");

    const double upperlimpar = gMinuit->fErp[1];
    upperlimcount = upperlimpar*20.20/log(2)/hfit->GetBinWidth(1);
  }
  const double upperlimprob = upperlimcount/livetime/n_c13cap/b12eff/be12eff;

  printf("Upper limit count = %.2f\n", upperlimcount);
  printf("%sUpper limit prob = %.3f%%%s\n", RED, upperlimprob*100, CLR);
}
