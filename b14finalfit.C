{
  gROOT->ProcessLine(".L ~/root_code_fragments/distfromtarg.C");

  const bool fit = true;

  TFile * fiel = new TFile("/cp/s4/strait/fullfido-100s-3-25MeV-20141006-earlymich.root", "read");
  TTree * t = (TTree *) fiel->Get("t");

  TCanvas * c = new TCanvas("c1", "c1", 1000, 1000);
/*  c->Divide(2, 1);
  c->cd(1)->Divide(1, 2);
  c->cd(1)->cd(2)->Divide(2, 2);
  c->cd(2)->Divide(2,4); */

  const char * const RED = "\033[31;1m"; // bold red
  const char * const CLR = "\033[m"    ; // clear

  const char * const scutwoe = "!earlymich && nnear == 0 && miche > 0.7 && miche < 12 "
    "&& dist < 200 && ttlastvalid > 0.1 && ttlastmuon > 1 && timeleft > 100e3 ";

  const char * const acrcut = "1"; 

  const char * const cute = "e > 5 && e < 12";

  const char * const cutb12 = "b12like < 0.4";

  const char * const cutx[3] = {
   "1", //"mz > -1500",
   "1",
   //"(mx**2 + my**2 < 1400**2 || "
   //"(atan2(my,mx) > -1.7+3.14159/2 && atan2(my,mx) < -1.7+3*3.14159/2))",
   "1"
  };

  const char * const cutxx = string(cutx[0]).c_str(); 

  const string scut = Form("%s && %s && %s && %s", scutwoe, cutb12, cute, acrcut);

  const char * xscut = string(Form("%s && dt > 4e3 && dt < 20e3 && %s && %s && %s",
                                   scutwoe, cute, cutb12, acrcut)).c_str();

/*const char* drawx ="dx";
  const char* drawy ="dy";
  const char* drawz ="dz";
  const char* drawsx="dx*sqrt(1708^2-dx^2)+1708^2*atan(dx/sqrt(1708^2-dx^2))";
  const char* drawsy="dy*sqrt(1708^2-dy^2)+1708^2*atan(dy/sqrt(1708^2-dy^2))";
  const char* drawr ="dx**2+dy**2"; */
  const char* drawx ="mx";
  const char* drawy ="my";
  const char* drawz ="mz";
  const char* drawsx="mx*sqrt(1708^2-mx^2)+1708^2*atan(mx/sqrt(1708^2-mx^2))";
  const char* drawsy="my*sqrt(1708^2-my^2)+1708^2*atan(my/sqrt(1708^2-my^2))";
  const char* drawr ="mx**2+my**2";

/*
  c->cd(2)->cd(1);t->Draw(Form("%s:%s",drawz,drawr),xscut,".");
  c->cd(2)->cd(2);t->Draw(Form("%s:%s",drawz,drawr),Form("%s&&%s",xscut,cutxx),".");
  c->cd(2)->cd(3);t->Draw(Form("%s:%s",drawz,drawsy),xscut,".");
  c->cd(2)->cd(4);t->Draw(Form("%s:%s",drawz,drawsy),Form("%s&&%s",xscut,cutxx),".");
  c->cd(2)->cd(5);t->Draw(Form("%s:%s",drawz,drawsx),xscut,".");
  c->cd(2)->cd(6);t->Draw(Form("%s:%s",drawz,drawsx),Form("%s&&%s",xscut,cutxx),".");
  c->cd(2)->cd(7);t->Draw(Form("%s:%s",drawy,drawx),xscut,".");
  c->cd(2)->cd(8);t->Draw(Form("%s:%s",drawy,drawx),Form("%s&&%s",xscut,cutxx),".");
*/

//  c->cd(1)->cd(1);
  printf("Drawing...\n");
  t->Draw("dt/1000 >> hfit(1000, 0.001, 100)", scut.c_str());

 // This is the n16 efficiency.  The energy efficiency is made up.
 // 67% of the events should be nearly 100% efficient because they give a 6.1MeV gamma.
 // The other 33% is maybe 50% efficient given the 10.5 MeV endpoint and 5 MeV cut.
 // What is the be11 efficiency?
  const double n16eff = 1
    * 0.981 // subsequent muons
    * 0.977 // previous muons
    * 0.565 // delta r for 200mm, maybe not valid for the 6.1MeV gammas, so perhaps should be lower
    * 0.9709 // 100s from end of run
    * 0.84 // energy, see above
    * 0.986 // from ttlastvalid cut, very naive
    * 0.96 // from ttlastmuon cut, vary naive
    * (1-0.00504 * 20.20/178.3) // from b12like cut, using li-9 as a guide
    * 1.   // from mz cut: handled in the denominator
  ;

  const double be11eff = 1
    * 0.981 // subsequent muons
    * 0.977 // previous muons
    * 0.565 // delta r for 200mm
    * 0.9709 // 100s from end of run
    * 0.59 // energy, estimated from scaled b12 mc
    * 0.986 // from ttlastvalid cut, very naive
    * 0.96 // from ttlastmuon cut, vary naive
    * (1-0.00504 * 20.20/178.3) // from b12like cut, using li-9 as a guide
    * 1 //(1600.+1500.)/(1600.+1786.) // from mz cut: 0.916
  ;

  printf("%sN-16 eff: %f; Be-11 eff: %f%s\n", RED, n16eff, be11eff, CLR);

  TF1 * ee = new TF1("ee",
               Form("%f*(%f*[0]/[1]*log(2)*exp(-x*log(2)/[1]) + " // n16
               "[2]*exp(-x*log(2)/0.0202) + " // b12
               "[3]*exp(-x*log(2)/0.8399) + " // li8
               "[4] +" // accidentals
               "%f*[5]/[6]*log(2)*exp(-x*log(2)/[6]))", hfit->GetBinWidth(1), n16eff, be11eff), 0, 100); // be11

  ee->SetParameters(0.1, 7.13, 3000, 2, 0.1, 1);
  ee->FixParameter(0, 0);
  ee->FixParameter(1, 7.13);
  ee->FixParameter(6, 13.81);
  ee->FixParameter(5, 0);
  if(fit) printf("fitting...\n"), hfit->Fit("ee", "l");  // fit with no n16 or be11

  double likenon16 = 0; if(fit) likenon16 = gMinuit->fAmin;
  ee->ReleaseParameter(0);
  ee->ReleaseParameter(5);
  ee->SetParLimits(0, 0, 2000);
  ee->SetParLimits(5, 0, 2000);
  if(fit) printf("fitting...\n"), hfit->Fit("ee", "l"); // fit with n16 and be11

  double likewn16 = 0; if(fit) likewn16 = gMinuit->fAmin;

  if(likenon16 > likewn16)
    printf("%sSignificance for n16 and/or be11 rather than neither %f%s\n", RED,
           sqrt(2)*sqrt(likenon16-likewn16), CLR);

  ee->ReleaseParameter(1);
  ee->FixParameter(5, 0);
  if(fit) printf("fitting...\n"), hfit->Fit("ee", "le");  // fit with free lifetime n16

  if(fit)
    printf("%sHalf life %f +%f %f%s\n", RED, ee->GetParameter(1),
           gMinuit->fErp[1], gMinuit->fErn[1], CLR);

  ee->ReleaseParameter(5);
  ee->FixParameter(1, 7.13);
/*
  ee->ReleaseParameter(6);
  ee->SetParLimits(6, 8, 2000);
  if(fit) printf("fitting...\n"), hfit->Fit("ee", "le");  // fit with free lifetime be11

  if(fit)
    printf("%sHalf life %f +%f %f%s\n", RED, ee->GetParameter(6),
           gMinuit->fErp[5], gMinuit->fErn[5], CLR);

  ee->FixParameter(6, 13.81);
*/

 
  ee->FixParameter(0, 0);
  if(fit) printf("fitting...\n"), hfit->Fit("ee", "le");  // fit with not n16 and be11

  const double be11norm = ee->GetParameter(5);

  printf("Drawing...\n");
  t->Draw("dt/1000 >> hdisp(49, 2, 100)", scut.c_str(), "hist");
  if(hdisp->GetBinContent(2) > 5) hdisp->Draw("e");

  ee->ReleaseParameter(0);
  ee->FixParameter(5, 0);
  if(fit) printf("fitting...\n"), hfit->Fit("ee", "le");  // fit with n16 and no be11

  TF1 * eedisp = (TF1 *)ee->Clone("eedisp");
  eedisp->SetNpx(400);
  eedisp->SetLineColor(kRed);

  bool nearlimit[5] = {false, false, false, false, false};
  int tomult[5] = { 0, 2, 3, 4, 5 };
  const double mult = hdisp->GetBinWidth(1)/hfit->GetBinWidth(1);
  for(int i = 0; i < 5; i++){
    //if(eedisp->GetParameter(tomult[i]) < eedisp->GetParError(tomult[i]))
    if(tomult[i] == 5)
      nearlimit[i] = true;

    eedisp->SetParameter(tomult[i], ee->GetParameter(tomult[i])*mult);
    eedisp->SetParError( tomult[i], ee->GetParError (tomult[i])*mult);
  }

  hdisp->Draw("e");
  eedisp->Draw("same");
  c->Update();

  TF1* acc= new TF1("acc", Form("%f*[0]",                                             hfit->GetBinWidth(1)), 0, 100);
  TF1* b12= new TF1("b12", Form("%f*([0]*exp(-x*log(2)/0.0202) + [1])",               hfit->GetBinWidth(1)), 0, 100);
  TF1* li8= new TF1("li8", Form("%f*([0]*exp(-x*log(2)/0.8399) + [1])",               hfit->GetBinWidth(1)), 0, 100);
  TF1* n16= new TF1("n16", Form("%f*(%f*[0]/ 7.13*log(2)*exp(-x*log(2)/ 7.13) + [1])",hfit->GetBinWidth(1), n16eff), 0,100);
  TF1* be11=new TF1("be11",Form("%f*(%f*[0]/13.81*log(2)*exp(-x*log(2)/13.81) + [1])",hfit->GetBinWidth(1), be11eff),0,100);

  b12->SetNpx(400);
  n16->SetNpx(400);
  be11->SetNpx(400);
  li8->SetNpx(400);

  TF1 * parts[5] = { b12, li8, n16, acc, be11 };

  n16->SetParameter(0, eedisp->GetParameter(0));
  n16->SetParError(0, eedisp->GetParError(0));
  n16->SetParameter(1, eedisp->GetParameter(4));
  n16->SetParError(1, eedisp->GetParError(4));

  b12->SetParameter(0, eedisp->GetParameter(2));
  b12->SetParError(0, eedisp->GetParError(2));
  b12->SetParameter(1, eedisp->GetParameter(4));
  b12->SetParError(1, eedisp->GetParError(4));

  li8->SetParameter(0, eedisp->GetParameter(3));
  li8->SetParError(0, eedisp->GetParError(3));
  li8->SetParameter(1, eedisp->GetParameter(4));
  li8->SetParError(1, eedisp->GetParError(4));

  acc->SetParameter(0, eedisp->GetParameter(4));
  acc->SetParError(0, eedisp->GetParError(4));

  be11->SetParameter(0, eedisp->GetParameter(5));
  be11->SetParError(0, eedisp->GetParError(5));
  be11->SetParameter(1, eedisp->GetParameter(4));
  be11->SetParError(1, eedisp->GetParError(4));

  
  for(int i = 0; i < 5; i++){
    parts[i]->SetLineStyle(7);
    parts[i]->SetLineWidth(2);
    if(nearlimit[i]){
      printf("plot %d (%d) is near limit\n", i, tomult[i]);
      parts[i]->SetParameter(0, be11norm*mult); // XXX parts[i]->GetParameter(0)+
                                //parts[i]->GetParError(0));
      parts[i]->SetLineStyle(kDotted);
      parts[i]->SetLineWidth(1);
    }
    parts[i]->Draw("Same");
  } 

  const double Nfound = ee->GetParameter(0);
  double Nerrup, Nerrlo;
  char * errtype = NULL;
  if(hfit->GetEntries() < 3){
    errtype = "HESSE";
    Nerrup = Nfound * ee->GetParError(0)/ee->GetParameter(0);
    Nerrlo = Nfound * ee->GetParError(0)/ee->GetParameter(0);
  }
  else{
    errtype = "MINOS";
    if(fit) Nerrup = Nfound * gMinuit->fErp[0]/ee->GetParameter(0);
    if(fit) Nerrlo = Nfound * gMinuit->fErn[0]/ee->GetParameter(0);
  }

  printf("%sN N-16 found assuming no Be-11: %f +%f %f %s%s\n",
         RED, Nfound, Nerrup, Nerrlo, errtype, CLR);

  const double Ocaptures = (0.2 + 5.67 + 0.3 + 1.248) * 489.509;
  const double Ccaptures = 358 * 489.509;

  const double toprob = 1./Ocaptures;

  printf("%sProb: %g +%g %g%s\n", 
      RED, toprob*Nfound, toprob*Nerrup, toprob*Nerrlo, CLR);


/*  c->cd(1)->cd(2)->cd(1)->SetLogy();

  const string escut = Form("%s && %s && %s && %s && dt > 4e3 && dt < 20e3",
                            scutwoe, cutxx, cutb12, acrcut);
  printf("Drawing...\n");
  t->Draw("e >> ehist(50, 4, 14)", escut.c_str(), "e");

  c->cd(1)->cd(2)->cd(2);

  const string b12scut = Form("%s && %s && %s && %s && dt > 4e3 && dt < 20e3",
                              scutwoe, cutxx, cute, acrcut);

  printf("Drawing...\n");
  t->Draw("b12like >> b12hist(500, 0, 1)", b12scut.c_str(), "e");
*/
  ee->ReleaseParameter(5);
  ee->ReleaseParameter(0);
 // ee->SetParLimits(5, 0, 2000);
 // ee->SetParLimits(0, 0, 2000);

  new TCanvas;

  hfit->Fit("ee", "le");
  
  gMinuit->Command("set print 0");

  gMinuit->fUp = 2.3/2; // 90% in 1D
  gMinuit->Command("mncont 1 6 200");
  TGraph * ninty_1d = ((TGraph*)gMinuit->GetPlot())->Clone();

  double be11lim = 0;
  for(int i = 0; i < ninty_1d->GetN(); i++)
    if(ninty_1d->GetY()[i] > be11lim)
      be11lim = ninty_1d->GetY()[i];

  printf("%sBe-11 90%% CL upper limit events: %.3lf %s\n", RED, be11lim, CLR);
  printf("%sBe-11 90%% CL upper limit prob: %.3lg %s\n", RED, be11lim/Ccaptures, CLR);

  gMinuit->fUp = 4.61/2; // 90% CL contour in 2D
  gMinuit->Command("mncont 1 6 200");
  TGraph * ninty_2d = ((TGraph*)gMinuit->GetPlot())->Clone();

  gMinuit->fUp = 11.83/2; // 99.73% contour in 2D
  gMinuit->Command("MNC 1 6 200");
  TGraph * ninty983_2d = ((TGraph*)gMinuit->GetPlot())->Clone();

  ninty983_2d->Draw("al");
  ninty_2d->Draw("l");
  ninty_1d->Draw("l");

}
