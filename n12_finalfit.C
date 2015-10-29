#include "consts.h"
#include "noncarbondenominators_finalfit.out.h"

void n12_finalfit()
{
  const int nncut = 3, nncuthigh = 4;
  
  // Relative amounts of effective oxygen in these regions
  const double o_targ = mass_o16targ,
               o_targacrlyic = mass_o16targves+mass_o16targbits,
               o_gc = mass_o16gc,

               // use the beta number here, since I'm going to apply
               // different n efficiencies
               o_gcacrlyic = mass_o16gcves_effective;

  const double o_sum = o_targ+o_targacrlyic+o_gc+o_gcacrlyic;

  const double
    targf    =                             o_targ/o_sum,
    targedgef= (o_targacrlyic*85./(85.+58.)     )/o_sum,
    gcf      = (o_targacrlyic*58./(85.+58.)+o_gc)/o_sum,
    gcedgef  =                        o_gcacrlyic/o_sum;


  const double tpneff = 0.55, tpedgeneff = 0.748,
               gpneff = 0.9323, gpedgeneff = 0.1843 /* vary this */;

  const double neff = tpneff*targf
                    + tpedgeneff*targedgef
                    + gpneff*gcf
                    + gpedgeneff*gcedgef;


  const double eff = 1
    * exp(-1.*log(2)/11.00) // n12 half-life and 1ms veto
    * (1 - exp(-100.*log(2)/11.00)) // 100ms window
    * light_noise_eff
    * mich_eff
    * sub_muon_eff05 // subsequent muons
    * 0.897 // delta r
    * 0.99709 // 10s from end of run
    * 0.9481 // energy
    * neff
  ;

  const double captures = n_o16cap_beta * livetime;

  const double toprob = 1./captures/eff;

  printf("Efficiencies: neutron, total: %.3f (%.3f, %.3f, %.3f, %.3f), %.3f\n",
         neff, tpneff, tpedgeneff, gpneff, gpedgeneff, eff);


  TFile * fiel = new TFile(rootfile3up, "read");
  TTree * t = (TTree *) fiel->Get("t");

  TCanvas * c = new TCanvas(Form("c%d", nncut), Form("c%d", nncut));

  const char * const ndef = "(latennear+ngdnear-latengdnear)";

  const string scut =
  Form("!earlymich && miche < 12 && dist < 400 && %s >= %d && %s <= %d && e > 4"
    "&& e < 18 && timeleft > 10e3", ndef, nncut, ndef, nncuthigh);

  const char * const cut = scut.c_str();

  const int nsel = t->GetEntries(Form("%s && dt < 100", cut));
  printf("Number selected in 100ms: %d\n", nsel);

  if(nsel > 0){
    t->Draw(Form("dt/1000 >> hfit%d(10000, 0.001, 100)", nncut), cut);
    TH1 * hfit = gROOT->FindObject(Form("hfit%d", nncut));

    TF1 * ee = new TF1(Form("ee%d", nncut), "[0]*exp(-x*log(2)/0.0202) + "
                 "[1]*exp(-x*log(2)/[2]) + "
                 "[3]*exp(-x*log(2)/7.13) + "
                 "[4]", 0, 100);

    ee->SetParameters(1, 1, 0.7700, 1, 1);
    ee->FixParameter(3, 0);
    ee->FixParameter(2, 0.0110);

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

    t->Draw(Form("dt/1000 >> hdisp%d(200, 0.001, 2.001)", nncut), cut, "hist");
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

    TF1 * b12 = new TF1(Form("b12", nncut), "[0]*exp(-x*log(2)/0.0202)" , 0, 100);
    TF1 * b8 = new TF1(Form("b8", nncut), "[0]*exp(-x*log(2)/[1])", 0, 100);
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

    printf("TECHNOTE: Selected >0 events for N-12, look at the code\n");
  }
  else{
    printf("%sTECHNOTE results.tex probTwelveNFromSixteenO: N-12 limit < %.0e\n%s",
          RED, 2.30258509299404590/eff/captures/lim_inflation_for_obeta, CLR);
  }
}
