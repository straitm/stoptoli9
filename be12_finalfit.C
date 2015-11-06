#include <math.h>
#include <stdio.h>
#include "TTree.h"
#include "TF1.h"
#include "TH1.h"
#include "TFile.h"
#include "TMinuit.h"
#include "consts.h"
#include "sub_muon_eff.out.h"
#include "totallivetime_finalfit.out.h"
#include "distcuteff_wholeloose_finalfit.out.h"
#include "carbondenominators_finalfit.out.h"
#include "b12cutefficiency_finalfit.out.h"

void be12_finalfit()
{
  const double be12timecut = 100.;
  const double b12timecut = 250.;

  const double be12eff = 
      0.894 // Energy over 3MeV
    * wholedet_dist400eff 
    * (exp(-1*log(2)/be12life) - exp(-be12timecut*log(2)/be12life))
    * light_noise_eff
    * mich_eff
    * sub_muon_eff05 // subsequent muons
  ;


  const double b12eff = 
      b12energyeff3MeV
    * 0.97 // dist2: from doc-5679 -- conservative given higher energies here
    * (exp(-1*log(2)/b12life) - exp(-b12timecut*log(2)/b12life))  // After 1ms
    * light_noise_eff
    * mich_eff
    * sub_muon_eff05 // subsequent muons
   ;

  printf("%sBe-12 eff: %.1f%%\nB-12 eff:  %.1f%%\nTotal eff: %.1f%%%s\n",
         RED, be12eff*100, b12eff*100, b12eff*be12eff*100, CLR);

  TFile * feil = new TFile(rootfile_be12, "read");
  if(!feil || feil->IsZombie()){
    fprintf(stderr, "Couldn't read %s\n", rootfile_be12);
    exit(1);
  }
  TTree * t = dynamic_cast<TTree *>(feil->Get("t"));
  if(!t){
    fprintf(stderr, "Couldn't get TTree 't' out of %s\n", rootfile_be12);
    exit(1);
  }

  char cut[1000];
  sprintf(cut, "miche < 12 && latennear == 0 && "
    "dist < 400 && dist2 < 400 && "
    "e < 12 && e2 < 14 && dt < %f && dt2 < %f && dt2-dt > 1",
    be12timecut, b12timecut);

  TH1D * hfit = new TH1D("hfit", "", 249, 1, 250);
  t->Draw("dt2-dt >> hfit", cut, "e");
  const int n = t->Scan("dt:e:dt2-dt:e2", cut);
  double upperlimcount;
  if(n == 0){
    upperlimcount = 2.30258509299404590;
  }
  else{
    hfit=hfit;

    TF1 * ee = new TF1("ee",
      Form("[0] + [1]*exp(-x*log(2)/%f)", b12life), 0, 500);
    ee->SetNpx(400);

    ee->SetParLimits(1, 0, hfit->GetBinWidth(1));
    ee->SetParLimits(0, 0, hfit->GetBinWidth(1));

    hfit->Fit("ee", "liem", "e");
    gMinuit->fUp = 2.71/2; // 90% in 1D
    gMinuit->Command("MINOS 10000 2");

    const double upperlimpar = gMinuit->fErp[1];
    upperlimcount = upperlimpar*b12life/log(2)/hfit->GetBinWidth(1);
  }
  const double upperlimprob = upperlimcount/livetime/n_c13cap/b12eff/be12eff;

  printf("TECHNOTE 9: Be-12 count = %d\n", n);
  printf("TECHNOTE \\newcommand{probTwelveBeFromThirteenC}"
         "{$<%.2f\\%%$}\n", upperlimprob*100);
}
