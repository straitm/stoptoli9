{
  const char * const RED     = "\033[31;1m"; // bold red
  const char * const CLR      = "\033[m"    ; // clear

  TFile fiel("/cp/s4/strait/fullfido-100s-3-25MeV-20140925-earlymich.root");
  TTree * t = (TTree *) fiel->Get("t");

  const char * const cut =
    "!earlymich && miche < 12 && dist < 400 && nnear == 0 && e > 5"
    "&& e < 14 && timeleft > 100e3";

  TCanvas c1;
  c1.SetLogy();

  t.Draw("dt/1000 >> hfit(10000, 0.001, 100)", cut);

  TF1 ee("ee", "[0]*exp(-x*log(2)/0.0202) + "
               "[1]*exp(-x*log(2)/[2]) + "
               "[3]*exp(-x*log(2)/7.13) + "
               "[4]", 0, 100);

  ee.SetParameters(1, 1, 0.8399, 1, 1);
  ee.FixParameter(2, 0.8399);
  hfit->Fit("ee", "l");
  ee.ReleaseParameter(2);
  hfit->Fit("ee", "le");

  printf("%sli8 lifetime: %f +%f %f%s\n", RED,
         ee->GetParameter(2), gMinuit.fErp[2], gMinuit.fErn[2], CLR);

  ee.FixParameter(2, 0.8399);

  hfit->Fit("ee", "le");

  t.Draw("dt/1000 >> hdisp(80, 0.1, 20.1)", cut, "e");

  TF1 * eedisp = ee.Clone("eedisp");
  eedisp->SetNpx(400);
  eedisp->SetLineColor(kRed);

  int tomult[4] = { 0, 1, 3, 4};
  const double mult = hdisp->GetBinWidth(1)/hfit->GetBinWidth(1);
  for(int i = 0; i < 4; i++)
    eedisp.SetParameter(tomult[i], eedisp.GetParameter(tomult[i])*mult);
  eedisp.Draw("same");

  TF1 b12("b12", "[0]*exp(-x*log(2)/0.0202)" , 0, 100);
  TF1 li8("li8", "[0]*exp(-x*log(2)/0.8399)", 0, 100);
  TF1 n16("n16", "[0]*exp(-x*log(2)/7.13)"  , 0, 100);
  TF1 acc("acc", "[0]", 0, 100);

  b12.SetNpx(400);

  TF1 * parts[4] = { &b12, &li8, &n16, &acc };


  b12.SetParameter(0, eedisp.GetParameter(0));
  li8.SetParameter(0, eedisp.GetParameter(1));
  n16.SetParameter(0, eedisp.GetParameter(3));
  acc.SetParameter(0, eedisp.GetParameter(4));

  hdisp->GetYaxis()->SetRangeUser(acc->GetParameter(0)/50,
                                  acc->GetParameter(0)*10);

  for(int i = 0; i < 4; i++){
    parts[i]->SetLineStyle(7);
    parts[i]->SetLineWidth(2);
    parts[i]->Draw("Same");
  } 

  const double Nfound = li8->Integral(0, 20)/hdisp->GetBinWidth(1);
  const double Nerrup = Nfound * gMinuit.fErp[1]/ee->GetParameter(1);
  const double Nerrlo = Nfound * gMinuit.fErn[1]/ee->GetParameter(1);

  printf("%sN found: %f +%f %f%s\n", RED, Nfound, Nerrup, Nerrlo, CLR);

  const double eff = 1
    * 0.981 // subsequent muons
    * 0.977 // previous muons
    * 0.948 // delta r
    * 0.9709 // 100s from end of run
    * 0.71 // energy
  ;

  const double captures = 358. * 489.509;

  const double toprob = 1./captures/eff;

  printf("%sProb: %g +%g %g%s\n", 
      RED, toprob*Nfound, toprob*Nerrup, toprob*Nerrlo, CLR);

return;

  TCanvas c2;

  const char * const cut =
    "!earlymich && miche < 12 && dist < 200 && nnear == 0"
    "&& timeleft > 100e3  && dt/1000 > 0.3 && dt/1000 < 3";

  t.Draw("e >> ehist(50, 0, 25)", cut, "e");

  TCanvas c3;

  const char * const cut =
    "!earlymich && dist < 300 && e > 4 && e < 14 && nnear == 0"
    "&& timeleft > 100e3  && dt/1000 > 0.3 && dt/1000 < 3";

  t.Draw("miche >> emichhist(80, 0, 20)", cut, "e");
  
}
