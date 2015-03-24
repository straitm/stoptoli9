#include "consts.h"

void b8finalfit(const int nncut = 4, const int nncuthigh = 4)
{
  const double
    tp = (neff_dt_targ+0.0726)*neff_dr_800_targ, // accepting early Gd-n
    gp = neff_dt_gc  * neff_dr_800_h, // since not accepting early H-n
    gf = 1-n_c12captarget/n_c12cap,
    tf =   n_c12captarget/n_c12cap;

  double neffs[5] = {
       gf*          pow(1-gp,4)+tf*          pow(1-tp,4),
    4*(gf*    gp   *pow(1-gp,3)+tf*    tp   *pow(1-tp,3)),
    6*(gf*pow(gp,2)*pow(1-gp,2)+tf*pow(tp,2)*pow(1-tp,2)),
    4*(gf*pow(gp,3)*   (1-gp)  +tf*pow(tp,3)*   (1-tp)),
       gf*pow(gp,4)            +tf*pow(tp,4)
  };

  double neff = 0;
  for(int i = nncut; i <= nncuthigh && i <= 4; i++)
    neff += neffs[i];

  const double eff = 1
    * 0.981 // subsequent muons
    * 0.977 // previous muons
    * wholedet_dist400eff // delta r
    * 0.9709 // 100s from end of run
    * 0.969 // energy
    * neff
  ;

  printf("neutron efficiency: %.3f\n", neff);
  printf("overall efficiency: %.3f\n",  eff);

  const double captures = n_c12cap * livetime;

  const double toprob = 1./captures/eff;
  TFile * fiel = new TFile(rootfile3up, "read");
  TTree * t = (TTree *) fiel->Get("t");

  TCanvas * c = new TCanvas(Form("c%d", nncut), Form("c%d", nncut));

  const char * const ndef = "(latennear+ngdnear-latengdnear)";

  const string scut =
  Form("!earlymich && miche < 12 && dist < 400 && %s >= %d && %s <= %d && e > 4"
    "&& e < 18 && timeleft > 100e3", ndef, nncut, ndef, nncuthigh);

  const char * const cut = scut.c_str();

  t->Draw(Form("dt/1000 >> hfit%d(10000, 0.001, 100)", nncut), cut);
  TH1 * hfit = gROOT->FindObject(Form("hfit%d", nncut));

  TF1 * ee = new TF1(Form("ee%d", nncut), "[0]*exp(-x*log(2)/0.0202) + "
               "[1]*exp(-x*log(2)/[2]) + "
               "[3]*exp(-x*log(2)/7.13) + "
               "[4]", 0, 100);

  ee->SetParameters(1, 1, 0.7700, 1, 1);
  ee->FixParameter(3, 0);
  ee->FixParameter(2, 0.7700); // b8
  //ee->FixParameter(2, 0.0110); // n12

  if(nncut >= 3){
    ee->SetParLimits(0, 0, 10);
    ee->SetParLimits(1, 0, 10);
    ee->SetParLimits(4, 0, 10);
  }
  int p0isfixed = 0;
  if(nncut >= 4){
    p0isfixed = 1;
    ee->FixParameter(0, 0);
    ee->FixParameter(4, 0);
  }

  if(hfit->GetEntries()){
    ee->FixParameter(1, 0);
    hfit->Fit(Form("ee%d", nncut), "l");
    ee->ReleaseParameter(1);
    hfit->Fit(Form("ee%d", nncut), "le");
  }

  if(ee->GetParameter(0) < 1e-6){
    ee->FixParameter(0, 0);
    p0isfixed = 1;
    if(hfit->GetEntries()) hfit->Fit(Form("ee%d", nncut), "le");
  }

  if(hfit->GetEntries()){
    t->Draw(Form("dt/1000 >> hdisp%d(400, 0.001, 20.001)", nncut), cut, "hist");
    TH1 * hdisp = gROOT->FindObject(Form("hdisp%d", nncut));
    if(hdisp->GetBinContent(2) > 5) hdisp->Draw("e");


    TF1 * eedisp = ee->Clone(Form("eedisp%d", nncut));
    eedisp->SetNpx(400);
    eedisp->SetLineColor(kRed);

    int tomult[4] = { 0, 1, 3, 4};
    const double mult = hdisp->GetBinWidth(1)/hfit->GetBinWidth(1);
    for(int i = 0; i < 4; i++)
      eedisp->SetParameter(tomult[i], eedisp->GetParameter(tomult[i])*mult);
    eedisp->Draw("same");

    TF1 * b12 = new TF1(Form("b12", nncut), "[0]*exp(-x*log(2)/0.0202)", 0, 100);
    TF1 * b8  = new TF1(Form("b8", nncut), "[0]*exp(-x*log(2)/[1])", 0, 100);
    TF1 * acc = new TF1(Form("acc", nncut), "[0]", 0, 100);

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
    char * errtype = NULL;
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
        RED, Nfound, Nerrup, Nerrlo, errtype, CLR);


    printf("%sProb: %g +%g %g%s\n", 
        RED, toprob*Nfound, toprob*Nerrup, toprob*Nerrlo, CLR);
  }
  else{
    const double N = 2.3026;
    printf("%sProbability 90%% limit: %.2g%s\n",RED, toprob*N, CLR);
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
    "&& timeleft > 100e3", ndef, nncut, ndef, nncuthigh);
  const char * const ecut = escut.c_str();

  t->Draw(Form("e >> ehist%d(250, 0, 25)", nncut), ecut, "e");
*/
}
