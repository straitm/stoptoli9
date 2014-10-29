{
  TFile fiel("/cp/s4/strait/fullfido-100s-3-25MeV-20141022.root");
  TTree * t = (TTree *) fiel.Get("t");

  const int nncut = 1;
  TCanvas * c = new TCanvas(Form("c%d", nncut), Form("c%d", nncut));
  c->Divide(2, 1);
  c->cd(1)->SetLogy();

  const char * const RED = "\033[31;1m"; // bold red
  const char * const CLR = "\033[m"    ; // clear

  const char * const cut =
       "!earlymich && miche < 12 && dist < 400 &&"
       "(latennear+ngdnear-latengdnear) == 1 && e > 5"
       "&& e < 12 && timeleft > 10e3 && ttlastvalid > 1 && ttlastmuon > 1";

  t->Draw(Form("dt/1000 >> hfit%d(10000, 0.001, 10)", 1), cut);
  TH1 * hfit = gROOT->FindObject(Form("hfit%d", 1));

  TF1 ee(Form("ee%d", 1), "[0]*exp(-x*log(2)/0.0202) + "
               "[1]*exp(-x*log(2)/[2]) + "
               "[3]*exp(-x*log(2)/0.8399) + "
               "[4]", 0, 10);

  ee.SetParameters(1, 1, 0.1783, 1, 1);
  ee.FixParameter(2, 0.1783);
  ee.FixParameter(3, 0); //XXX

  hfit->Fit(Form("ee%d", nncut), "l");
  ee.ReleaseParameter(2);
  hfit->Fit(Form("ee%d", nncut), "le");
  printf("%sli9 lifetime: %f +%f %f%s\n", RED,
     ee->GetParameter(2), gMinuit.fErp[2], gMinuit.fErn[2], CLR);

  ee.FixParameter(2, 0.1783);


  hfit->Fit(Form("ee%d", nncut), "le");

  t->Draw(Form("dt/1000 >> hdisp%d(100, 0.001, 10.001)", nncut), cut, "e");
  TH1 * hdisp = gROOT->FindObject(Form("hdisp%d", nncut));

  TF1 * eedisp = ee.Clone(Form("eedisp%d", nncut));
  eedisp->SetNpx(400);
  eedisp->SetLineColor(kRed);

  int tomult[4] = { 0, 1, 3, 4};
  const double mult = hdisp->GetBinWidth(1)/hfit->GetBinWidth(1);
  for(int i = 0; i < 4; i++)
    eedisp.SetParameter(tomult[i], eedisp.GetParameter(tomult[i])*mult);
  eedisp.Draw("same");

  TF1 b12(Form("b12", nncut), "[0]*exp(-x*log(2)/0.0202)" , 0, 10);
  TF1 li9(Form("b8", nncut), "[0]*exp(-x*log(2)/[1])", 0, 10);
  TF1 li8(Form("li8", nncut), "[0]*exp(-x*log(2)/0.8399)" , 0, 10);
  TF1 acc(Form("acc", nncut), "[0]", 0, 10);

  b12.SetNpx(400);
  b8.SetNpx(400);

  TF2 * parts[4] = { &b12, &li9, &li8, &acc };

  b12.SetParameter(0, eedisp.GetParameter(0));
  b8 .SetParameter(0, eedisp.GetParameter(1));
  b8 .SetParameter(1, eedisp.GetParameter(2));
  li8.SetParameter(0, eedisp.GetParameter(3));
  acc.SetParameter(0, eedisp.GetParameter(4));

  for(int i = 0; i < 4; i++){
    parts[i]->SetLineStyle(7);
    parts[i]->SetLineWidth(2);
    parts[i]->Draw("Same");
  } 

  const double Nfound = b8->Integral(0, 20)/hdisp->GetBinWidth(1);
  const double Nerrup = Nfound * gMinuit.fErp[1]/ee->GetParameter(1);
  const double Nerrlo = Nfound * gMinuit.fErn[1]/ee->GetParameter(1);

  printf("%sN found: %f +%f %f%s\n", RED, Nfound, Nerrup, Nerrlo, CLR);

  const double tp = 0.64,
               gp = 0.93-0.0303, // since not accepting early nH
               gf = 0.688, tf = 1-0.688;

  neff = gf*gp + tf*tp;

  const double eff = 1
    * 0.981 // subsequent muons
    * 0.977 // previous muons
    * 0.948 // delta r
    * 0.99709 // 10s from end of run
    * 0.8 // energy XXX
    * neff
  ;

  const double captures = 358. * 489.509;

  const double toprob = 1./captures/eff;

  printf("%sProb: %g +%g %g%s\n", 
      RED, toprob*Nfound, toprob*Nerrup, toprob*Nerrlo, CLR);

  TF1 gaus("gaus", "gaus(0)", 0, 20);
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
  }

  c->cd(2)->SetLogy();

  const string escut =
    Form("!earlymich && miche < 12 && dist < 400 && (latennear+ngdnear-latengdnear) == 1"
    "&& timeleft > 10e3 && ttlastvalid > 1 && ttlastmuon > 1");
  const char * const ecut = escut.c_str();

  t->Draw(Form("e >> ehist%d(250, 0, 25)", nncut), ecut, "e");
}
