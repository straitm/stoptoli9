#include <math.h>
#include <string>
using std::string;
#include "TTree.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TFile.h"
#include "TMinuit.h"
#include "TF1.h"
#include "TH1.h"
#include "consts.h"
#include "sub_muon_eff.out.h"
#include "neff_dt_finalfit.out.h"
#include "distcuteff_targetgc_finalfit.out.h"
#include "totallivetime_finalfit.out.h"
#include "noncarbondenominators_finalfit.out.h"

void c9_finalfit(const char elem = 'o')
{
  const int nncut = 4, nncuthigh = 6;

  // oxygen
  double targrate = mass_o16targ,
         targvesrate =  mass_o16targves,
         targbitsrate = mass_o16targbits,
         gcrate = mass_o16gc,
         gcvesrate = mass_o16gcves_effective; 

  // nitrogen
  if(elem != 'o'){
    targrate = mass_n14targ,
    targvesrate =  0,
    targbitsrate = 0,
    gcrate = mass_n14gc,
    gcvesrate = 0;
  }

  const double totalrate =
    gcvesrate+gcrate+targvesrate+targbitsrate+targrate;

  const double targf    =              targrate/totalrate,
               targedgef=           targvesrate/totalrate,
               gcf      = (targbitsrate+gcrate)/totalrate,
               gcedgef  =             gcvesrate/totalrate;

  const double spillin = 0.45; // could be 0.41 -- 0.49

  const double tpneff = n5of5eff_dt_dr_800_targ + n4of5eff_dt_dr_800_targ,
               gpneff = n5of5eff_dt_dr_800_gc + n4of5eff_dt_dr_800_gc,
               tpedgeneff = (tpneff+gpneff)/2, // not really right, but certainly it is between the above two
               gpedgeneff =  n5of5eff_dt_dr_800_gc*pow(spillin, 5)
                           + n4of5eff_dt_dr_800_gc*pow(spillin, 4)*5;

  const double neff = tpneff*targf + tpedgeneff*targedgef 
                    + gpneff*gcf   + gpedgeneff*gcedgef;

  const double deltar_eff = elem == 'o'?
  (mass_o16targ * targ_dist400eff +
   mass_o16targves * targves_dist400eff +
   (mass_o16targbits + mass_o16gc) * gc_dist400eff +
   mass_o16gcves_effective * gcves_dist400eff)/
  (mass_o16targ + mass_o16targves + mass_o16targbits +
   mass_o16gc + mass_o16gcves_effective)
  :
  (mass_n14targ * targ_dist400eff + mass_n14gc * gc_dist400eff)/
  (mass_n14targ + mass_n14gc);

  printf("TECHNOTE 10.2.3: Delta r efficiency for C-9 from %s is %.2f%%\n",
         elem=='o'?"O-16":"N-14", deltar_eff*100);

  const double eff = 1
    * exp(-1/c9life) // half-life and 1ms veto
    * (1-exp(-1000./c9life)) // up to 1s
    * light_noise_eff
    * mich_eff
    * sub_muon_eff05 // subsequent muons
    * deltar_eff
    * (livetime_s - num_runs*10.)/livetime_s
    * 0.969 // energy -- maybe a little optimistic: 60%
            // goes to the ground state, rest to excited
            // states between 2-3 MeV, which isn't too bad,
            // but will lower (raise?) the efficiency a little.
    * neff
  ;

  const double captures = (elem == 'o'?n_o16cap_beta:n_n14cap) * livetime;

  const double toprob = 1./captures/eff;

  printf("Efficiencies: neutron, total: %.3f (%.3f, %.3f, %.3f, %.3f), %.3f\n",
         neff, tpneff, tpedgeneff, gpneff, gpedgeneff, eff);


  TFile * fiel = new TFile(rootfile3up, "read");
  TTree * t = (TTree *) fiel->Get("t");

  TCanvas * c = new TCanvas(Form("c%d", nncut), Form("c%d", nncut));
  c->cd();

  const char * const ndef = "(latennear+ngdnear-latengdnear)";

  const string scut =
  Form("!earlymich && miche < 12 && dist < 400 && %s >= %d && %s <= %d && e > 4"
    "&& e < 18 && timeleft > 10e3", ndef, nncut, ndef, nncuthigh);

  const char * const cut = scut.c_str();

  const int nsel = t->GetEntries(Form("%s && dt < 1000", cut));
  printf("N selected: %d\n", nsel);

  if(nsel > 0){

    t->Draw(Form("dt/1000 >> hfit%d(1000, 0.001, 10)", nncut), cut);
    TH1 * hfit = dynamic_cast<TH1*>(gROOT->FindObject(Form("hfit%d", nncut)));

    TF1 * ee = new TF1(Form("ee%d", nncut), "[0]*exp(-x*log(2)/0.0202) + "
                 "[1]*exp(-x*log(2)/[2]) + "
                 "[3]*exp(-x*log(2)/7.13) + "
                 "[4]", 0, 100);

    ee->SetParameters(1, 1, 0.7700, 1, 1);
    ee->FixParameter(3, 0);
    ee->FixParameter(2, 0.1270);

    ee->SetParLimits(1, 0, 10);
    if(nncut >= 3){
      ee->SetParLimits(0, 0, 10);
      ee->SetParLimits(4, 0, 10);
    }
    int p0isfixed = 0;
    if(nncut >= 4){
      p0isfixed = 1;
      ee->FixParameter(0, 0);
      ee->FixParameter(4, 0);
    }

    hfit->Fit(Form("ee%d", nncut), "le");

    if(ee->GetParameter(0) < 1e-6){
      p0isfixed = 1;
      ee->FixParameter(0, 0);
      hfit->Fit(Form("ee%d", nncut), "le");
    }

    t->Draw(Form("dt/1000 >> hdisp%d(20, 0.001, 10.01)", nncut), cut, "hist");
    TH1 * hdisp = dynamic_cast<TH1*>(gROOT->FindObject(Form("hdisp%d", nncut)));
    if(hdisp->GetBinContent(2) > 5) hdisp->Draw("e");

    TF1 * eedisp = dynamic_cast<TF1*>(ee->Clone(Form("eedisp%d", nncut)));
    eedisp->SetNpx(400);
    eedisp->SetLineColor(kRed);

    int tomult[4] = { 0, 1, 3, 4};
    const double mult = hdisp->GetBinWidth(1)/hfit->GetBinWidth(1);
    for(int i = 0; i < 4; i++)
      eedisp->SetParameter(tomult[i], eedisp->GetParameter(tomult[i])*mult);
    eedisp->Draw("same");

    TF1 * b12 = new TF1("b12", "[0]*exp(-x*log(2)/0.0202)" , 0, 100);
    TF1 * b8  = new TF1("b8",  "[0]*exp(-x*log(2)/[1])", 0, 100);
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

    printf("TECHNOTE: c9 error.  Some events selected, so you may need"
      " to rewrite c9_finalfit.C\n");
  }
  else{
    if(elem == 'o')
      printf("%sTECHNOTE results.tex %s: "
           "If no events and no background: <%.2f%%%s\n", 
        RED, "probNineCFromSixteenO",
        2.30258509299404590*toprob*lim_inflation_for_obeta*100, CLR);

    else 
      printf("%sTECHNOTE results.tex %s: "
           "If no events and no background: < %.1f%%%s\n", 
        RED, "probNineCFromFourteenN",
        2.30258509299404590*toprob*lim_inflation_for_obeta*100, CLR);
  }


/*  TF1 gaus("gaus", "gaus(0)", 0, 20);
  gaus.SetParameters(1, toprob*Nfound, toprob*Nerrup);

  for(int i = 1; i < 400; i++){
    const double up = 0.01*i;
    const double frac = gaus.Integral(0, up)/gaus.Integral(0, 20);
    printf("%f %f\n", up, frac);
    if(frac > 0.9){
      printf("90%% limit = %f\n", up);
      printf("90%% limit prob = %f\n", up*toprob);
      printf("90%% limit prob *0.1/1.22 = %f\n", up*toprob *(1+0.1/1.22));
      break;
    }
  } */

/*
  c->cd(2);

  const string escut =
    Form("!earlymich && miche < 12 && dist < 400 && %s >= %d && %s <= %d "
    "&& timeleft > 10e3", ndef, nncut, ndef, nncuthigh);
  const char * const ecut = escut.c_str();

  t->Draw(Form("e >> ehist%d(250, 0, 25)", nncut), ecut, "e");
*/
}
