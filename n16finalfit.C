{
  TFile * fiel = new TFile("/cp/s4/strait/fullfido-100s-3-25MeV-20140925-earlymich.root", "read");
  TTree * t = (TTree *) fiel->Get("t");

  TCanvas * c = new TCanvas("c1", "c1");
  c->Divide(2, 1);
  c->cd(1)->Divide(1, 2);
  c->cd(1)->cd(1)->SetLogy();

  const char * const RED = "\033[31;1m"; // bold red
  const char * const CLR = "\033[m"    ; // clear

  const char * const scutwoe = "!earlymich && miche < 12 && dist < 200 "
    "&& ttlastmuon > 1 && ttlastvalid > 1 && "
    "nnear == 0 && timeleft > 100e3";

  const char * const cutx =
    "(abs(dz) < 1200 || dx**2+dy**2 > 500**2) && "
    "(dz > -1850 || dy < 1600) &&"
    "(abs(dz + 900) > 150 || abs(dy + 1150) > 150 || abs(dx - 50) > 150)";

  const string scut = Form("%s && %s && e > 5 && e < 10 && (e < 7.5 || e > 8.5)", scutwoe, cutx);

  t->Draw("dt/1000 >> hfit(10000, 0.001, 100)", scut.c_str());

  TF1 * ee = new TF1("ee",
               "[0]*exp(-x*log(2)/[1]) + "
               "[2]*exp(-x*log(2)/0.0202) + "
               "[3]*exp(-x*log(2)/0.8399) + "
               "[4]", 0, 100);

  ee->SetParameters(1, 1, 1, 1, 1);
  ee->FixParameter(1, 7.13);

  hfit->Fit("ee", "le");

  t->Draw("dt/1000 >> hdisp(200, 0.001, 100)", scut.c_str(), "hist");
  if(hdisp->GetBinContent(2) > 5) hdisp->Draw("e");
  c->Modified();
  c->Update();

  TF1 * eedisp = (TF1 *)ee->Clone("eedisp");
  eedisp->SetNpx(400);
  eedisp->SetLineColor(kRed);

  int tomult[4] = { 0, 2, 3, 4};
  const double mult = hdisp->GetBinWidth(1)/hfit->GetBinWidth(1);
  for(int i = 0; i < 4; i++)
    eedisp->SetParameter(tomult[i], eedisp->GetParameter(tomult[i])*mult);
  eedisp->Draw("same");
  c->Update();

  TF1 * b12 = new TF1("b12", "[0]*exp(-x*log(2)/0.0202)" , 0, 100);
  TF1 * li8 = new TF1("n16", "[0]*exp(-x*log(2)/0.8399)", 0, 100);
  TF1 * n16 = new TF1("n16", "[0]*exp(-x*log(2)/[1])", 0, 100);
  TF1 * acc = new TF1("acc", "[0]", 0, 100);

  b12->SetNpx(400);
  n16->SetNpx(400);

  TF1 * parts[4] = { b12, li8, n16, acc };

  n16->SetParameter(0, eedisp->GetParameter(0));
  n16->SetParameter(1, eedisp->GetParameter(1));
  b12->SetParameter(0, eedisp->GetParameter(2));
  li8->SetParameter(0, eedisp->GetParameter(3));
  acc->SetParameter(0, eedisp->GetParameter(4));

  for(int i = 0; i < 4; i++){
    parts[i]->SetLineStyle(7);
    parts[i]->SetLineWidth(2);
    parts[i]->Draw("Same");
  } 

  const double Nfound = n16->Integral(0, 20)/hdisp->GetBinWidth(1);
  double Nerrup, Nerrlo;
  char * errtype = NULL;
  if(hfit->GetEntries() < 3){
    errtype = "HESSE";
    Nerrup = Nfound * ee->GetParError(0)/ee->GetParameter(0);
    Nerrlo = Nfound * ee->GetParError(0)/ee->GetParameter(0);
  }
  else{
    errtype = "MINOS";
    Nerrup = Nfound * gMinuit.fErp[0]/ee->GetParameter(0);
    Nerrlo = Nfound * gMinuit.fErn[0]/ee->GetParameter(0);
  }

  printf("%sN found: %f +%f %f %s%s\n",
         RED, Nfound, Nerrup, Nerrlo, errtype, CLR);

  const double eff = 1
    * 0.981 // subsequent muons
    * 0.977 // previous muons
    * 0.565 // delta r
    * 0.9709 // 100s from end of run
    * 0.8 // energy XXX
    * 0.9 // spatial cuts XXX
  ;

  const double captures = (7.5 + 7.0 + 0.2 + 0.3) * 489.509;

  const double toprob = 1./captures/eff;

  printf("%sProb: %g +%g %g%s\n", 
      RED, toprob*Nfound, toprob*Nerrup, toprob*Nerrlo, CLR);

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

  c->cd(1)->cd(2);

  const string escut = Form("%s && %s && dt > 20e3", scutwoe, cutx);

  t->Draw("e >> ehist(250, 0, 25)", escut.c_str(), "e");

  c->Update();

  c->cd(2)->Divide(2,2);

  const string xscut = Form("%s && dt > 20e3", scutwoe);

  c->cd(2)->cd(1);
  t->Draw("dz:dx**2+dy**2", xscut.c_str(), "");

  c->cd(2)->cd(2);
  t->Draw("dz:dy",
    Form("%s && dt>20e3 && (abs(dz)<1200 || dx**2+dy**2>500**2)", scutwoe), "");

  c->cd(2)->cd(3);
  t->Draw("dy:dx",
    Form("%s && dt>20e3 && (abs(dz)<1200 || dx**2+dy**2>500**2) && (dz > -1850 || dy < 1600)", scutwoe), "");

  c->cd(2)->cd(4);
  t->Draw("dz:dx",
    Form("%s && dt>20e3 && %s", scutwoe, cutx), "");
}
