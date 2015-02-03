{
  const char * const RED     = "\033[31;1m"; // bold red
  const char * const CLR      = "\033[m"    ; // clear

  TFile fiel("/cp/s4/strait/fullfido-300s-3-25MeV-20141117.root");
  TTree * t = (TTree *) fiel->Get("t");

  const char * const cut =
    "!earlymich && miche < 12 && b12like < 0.02 && dist < 400 && latennear == 0 && e > 15"
    "&& e < 22 && timeleft > 10e3";

  TCanvas c1;

  t.Draw("dt/1000 >> hfit(10000, 0.001, 10)", cut);

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

  const double eff = 1
    * 0.981 // subsequent muons
    * 0.977 // previous muons
    * 0.944 // delta r
    * 0.99709 // 10s from end of run
    * 0.146 // energy
    * 0.906 // b12like
  ;

  const double captures = (1.1+6.2+0.1+2.1) * 489.509;

  const double toprob = 1./captures/eff;

  printf("%sEff: %f\nProb: %g +%g %g%s\n", 
      RED, eff, toprob*Nfound, toprob*Nerrup, toprob*Nerrlo, CLR);

  return;
  new TCanvas;

  gMinuit->Command("set print 0");

  gMinuit->fUp = 2.3/2; // 90% in 1D
  gMinuit->Command("mncont 1 2 200");
  TGraph * ninty_1d = ((TGraph*)gMinuit->GetPlot())->Clone();

  gMinuit->fUp = 4.61/2; // 90% CL contour in 2D
  gMinuit->Command("mncont 1 2 200");
  TGraph * ninty_2d = ((TGraph*)gMinuit->GetPlot())->Clone();

  gMinuit->fUp = 11.83/2; // 99.73% contour in 2D
  gMinuit->Command("MNC 1 2 200");
  TGraph * ninty983_2d = ((TGraph*)gMinuit->GetPlot())->Clone();

  ninty983_2d->Draw("al");
  ninty_2d->Draw("l");
  ninty_1d->Draw("l");
}
