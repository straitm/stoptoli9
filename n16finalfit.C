{
  gROOT->ProcessLine(".L ~/root_code_fragments/distfromtarg.C");

  const bool fit = true;

  TFile * fiel = new TFile("/cp/s4/strait/fullfido-100s-3-25MeV-20141006-earlymich.root", "read");
  TTree * t = (TTree *) fiel->Get("t");

  TCanvas * c = new TCanvas("c1", "c1", 1000, 1000);
  c->Divide(2, 1);
  c->cd(1)->Divide(1, 2);
  c->cd(1)->cd(2)->Divide(2, 2);
  c->cd(2)->Divide(2,4);

  const char * const RED = "\033[31;1m"; // bold red
  const char * const CLR = "\033[m"    ; // clear

  const char * const scutwoe = "!earlymich && nnear == 0 && miche < 12 "
    "&& dist < 400 &&  timeleft > 100e3";

  const char * const acrcut = "distfromacrylic(dx, dy, dz) < 100";

  const char * const cute = "e > 5 && e < 12 && (e < 7.4 || e > 8.4)";

  const char * const cutb12 = "b12like < 0.4";

  const char * const cutx[3] = {
   "dz > -(1229+200)",
   "(dx**2 + dy**2 < (1150+200)**2 || "
   "(atan2(dy,dx) > -1.7+3.14159/2 && atan2(dy,dx) < -1.7+3*3.14159/2))",
   "1"
  };

  const char * const cutxx = (string(cutx[1]) + "&&" + string(cutx[0])).c_str(); 

  const string scut = Form("%s && %s && %s && %s", scutwoe, cutb12, cute, acrcut);

  const char * xscut = string(Form("%s && dt > 10e3 && %s && %s && %s", scutwoe, cute, cutb12, acrcut)).c_str();

  c->cd(2)->cd(1);t->Draw("dz:dx**2+dy**2;r",xscut,".");
  c->cd(2)->cd(2);t->Draw("dz:dx**2+dy**2",Form("%s&&%s",xscut,cutxx),".");
  c->cd(2)->cd(3);t->Draw("dz:dy*sqrt(1708^2-dy^2)+1708^2*atan(dy/sqrt(1708^2-dy^2))",xscut,".");
  c->cd(2)->cd(4);t->Draw("dz:dy*sqrt(1708^2-dy^2)+1708^2*atan(dy/sqrt(1708^2-dy^2))",Form("%s&&%s",xscut,cutxx),".");
  c->cd(2)->cd(5);t->Draw("dz:dx*sqrt(1708^2-dx^2)+1708^2*atan(dx/sqrt(1708^2-dx^2))",xscut,".");
  c->cd(2)->cd(6);t->Draw("dz:dx*sqrt(1708^2-dx^2)+1708^2*atan(dx/sqrt(1708^2-dx^2))",Form("%s&&%s",xscut,cutxx),".");
  c->cd(2)->cd(7);t->Draw("dy:dx",xscut,".");
  c->cd(2)->cd(8);t->Draw("dy:dx",Form("%s&&%s",xscut,cutxx),".");

  c->cd(1)->cd(1)->SetLogy();

  t->Draw("dt/1000 >> hfit(10000, 0.001, 100)", scut.c_str());

  TF1 * ee = new TF1("ee",
               "[0]*exp(-x*log(2)/[1]) + "
               "[2]*exp(-x*log(2)/0.0202) + "
               "[3]*exp(-x*log(2)/0.8399) + "
               "[4] +"
               "[5]*exp(-x*log(2)/13.81)", 0, 100);

  ee->SetParameters(1, 7.13, 1, 1, 1, 1);
  ee->FixParameter(0, 0);
  ee->FixParameter(1, 7.13);
  ee->FixParameter(5, 0);
  if(fit) hfit->Fit("ee", "l");  // fit with no n16 or be11

  double likenon16 = 0; if(fit) likenon16 = gMinuit->fAmin;
  ee->ReleaseParameter(0);
  if(fit) hfit->Fit("ee", "l"); // fit with n16

  double likewn16 = 0; if(fit) likewn16 = gMinuit->fAmin;

  if(likenon16 > likewn16)
    printf("%sSignificance %f%s\n", RED,
           sqrt(2)*sqrt(likenon16-likewn16), CLR);

  ee->ReleaseParameter(1);
  if(fit) hfit->Fit("ee", "le");  // fit with free lifetime n16

  if(fit)
    printf("%sHalf life %f +%f %f%s\n", RED, ee->GetParameter(1),
           gMinuit->fErp[1], gMinuit->fErn[1], CLR);

  ee->FixParameter(1, 7.13);
  ee->ReleaseParameter(5);
  ee->SetParLimits(5, 0, 10);

  if(fit) hfit->Fit("ee", "le");  // fit with b16 and be11

  t->Draw("dt/1000 >> hdisp(200, 0.001, 100)", scut.c_str(), "hist");
  if(hdisp->GetBinContent(2) > 5) hdisp->Draw("e");
  c->Modified();
  c->Update();

  TF1 * eedisp = (TF1 *)ee->Clone("eedisp");
  eedisp->SetNpx(400);
  eedisp->SetLineColor(kRed);

  bool atlimit[5] = {false, false, false, false, false};
  int tomult[5] = { 0, 2, 3, 4, 5};
  const double mult = hdisp->GetBinWidth(1)/hfit->GetBinWidth(1);
  for(int i = 0; i < 5; i++){
    double setit = eedisp->GetParameter(tomult[i]);
    if(eedisp->GetParameter(tomult[i]) < eedisp->GetParError(tomult[i])){
      setit += eedisp->GetParError(tomult[i]);
      atlimit[i] = true;
    }

    eedisp->SetParameter(tomult[i], eedisp->GetParameter(tomult[i])*mult);
  }
  eedisp->Draw("same");
  c->Update();

  TF1 * b12 = new TF1("b12", "[0]*exp(-x*log(2)/0.0202)" , 0, 100);
  TF1 * li8 = new TF1("n16", "[0]*exp(-x*log(2)/0.8399)", 0, 100);
  TF1 * n16 = new TF1("n16", "[0]*exp(-x*log(2)/[1])", 0, 100);
  TF1 * acc = new TF1("acc", "[0]", 0, 100);
  TF1 * be11= new TF1("be11","[0]*exp(-x*log(2)/13.81)", 0, 100);

  b12->SetNpx(400);
  n16->SetNpx(400);
  be11->SetNpx(400);

  TF1 * parts[5] = { b12, li8, n16, acc, be11 };

  n16->SetParameter(0, eedisp->GetParameter(0));
  n16->SetParameter(1, eedisp->GetParameter(1));
  b12->SetParameter(0, eedisp->GetParameter(2));
  li8->SetParameter(0, eedisp->GetParameter(3));
  acc->SetParameter(0, eedisp->GetParameter(4));
  be11->SetParameter(0, eedisp->GetParameter(5));
  
  for(int i = 0; i < 5; i++){
    if(atlimit[i]) parts[i]->SetLineStyle(kDotted);
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
    if(fit) Nerrup = Nfound * gMinuit.fErp[0]/ee->GetParameter(0);
    if(fit) Nerrlo = Nfound * gMinuit.fErn[0]/ee->GetParameter(0);
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

  c->cd(1)->cd(2)->cd(1);

  const string escut = Form("%s && %s && %s && %s && dt > 10e3",
                            scutwoe, cutxx, cutb12, acrcut);

  t->Draw("e >> ehist(250, 0, 25)", escut.c_str(), "e");

  c->cd(1)->cd(2)->cd(2);


  const string b12scut = Form("%s && %s && %s && %s && dt > 10e3",
                            scutwoe, cutxx, cute, acrcut);

  t->Draw("b12like >> b12hist(100, 0, 1)", b12scut.c_str(), "e");

  c->cd(1)->cd(2)->cd(3);


  const string acrscut = Form("%s && %s && %s && %s && dt > 0e3",
                            scutwoe, cutxx, cutb12, cute );

  t->Draw("distfromacrylic(dx, dy, dz) >> acrhist(42, 0, 1300)", acrscut.c_str(), "e");
}
