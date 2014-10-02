{
  const char * const RED = "\033[31;1m"; // bold red
  const char * const CLR = "\033[m"    ; // clear

  // li11
  //const double shortlife = 0.00875;
  //const double betanprob = 0.83;

  // b13bn
  const double shortlife = 0.01733;
  const double betanprob = 0.0029;

  TTree tg("t", "t");
  TTree th("t", "t");
  tg.ReadFile("li9-20140925-earlymich.Gd.ntuple");
  th.ReadFile("li9-20140925-earlymich.H.ntuple");

  TF1 ee("ee", "[0]*exp(-x*log(2)/[1]) + [2]*exp(-x*log(2)/[3]) + [4]", 0, 10);
  ee.SetNpx(400);
  ee.SetLineColor(kRed);
  ee.FixParameter(1, 0.1783);
  ee.FixParameter(3, shortlife);
  ee.SetParLimits(0, 0, 10);
  ee.SetParLimits(2, 0, 10);

  TF1 e("e", "[0]*exp(-x*log(2)/[1])", 0, 10);
  
  const char * const cut = "n==0 && dist < 300  && miche < 12 && !earlymich";

  ////////////////////////////////////////////////////
  tg.Draw("dt/1000 >> hfitg(10000, 0.001, 10)", cut);

  hfitg->Fit("ee", "le", "");

  for(int i = 0; i < 2; i++) e->SetParameter(i, ee->GetParameter(i+2));

  const double Nrawg = e->Integral(0, 10)/hfitg.GetBinWidth(1);
  const double Nrawgup = Nrawg / ee.GetParameter(2) * gMinuit.fErp[1];
  const double Nrawglo = Nrawg / ee.GetParameter(2) * gMinuit.fErn[1];

  const double Nrawguph = Nrawg / ee.GetParameter(2) * ee.GetParError(2);
  const double Nrawgloh = Nrawg / ee.GetParameter(2) * ee.GetParError(2);

  printf("%sgd Nraw = %.2f +%f %f%s\n", RED, Nrawg, Nrawgup, Nrawglo, CLR);

  ////////////////////////////////////////////////////
  th.Draw("dt/1000 >> hfitg(10000, 0.001, 10)", cut);

  hfitg->Fit("ee", "le", "");

  for(int i = 0; i < 2; i++) e->SetParameter(i, ee->GetParameter(i+2));
 
  const double Nrawh = e->Integral(0, 10)/hfitg.GetBinWidth(1);
  const double Nrawhup = Nrawh/ ee.GetParameter(2) * gMinuit.fErp[1];
  const double Nrawhlo = Nrawh/ ee.GetParameter(2) * gMinuit.fErn[1];

  const double Nrawhuph = Nrawh / ee.GetParameter(2) * ee.GetParError(2);
  const double Nrawhloh = Nrawh / ee.GetParameter(2) * ee.GetParError(2);

  printf("%sH Nraw = %.2f +%f %f%s\n", RED, Nrawh, Nrawhup, Nrawhlo, CLR);

  ////////////////////////////////////////////////////
  th.Draw("dt/1000 >> hfitb(10000, 0.001, 10)", cut);
  tg.Draw("dt/1000 >> +hfitb", cut);

  hfitb->Fit("ee", "le", "");


  for(int i = 0; i < 2; i++) e->SetParameter(i, ee->GetParameter(i+2));

  const double Nrawb = e->Integral(0, 10)/hfitb.GetBinWidth(1);
  const double Nrawbup = Nrawb / ee.GetParameter(2) * gMinuit.fErp[1];
  const double Nrawblo = Nrawb / ee.GetParameter(2) * gMinuit.fErn[1];

  printf("%sBoth Nraw = %.2f +%f %f%s\n", RED, Nrawb, Nrawbup, Nrawblo, CLR);

  TF1 * eewith = (TF1 *)ee->Clone("eewith");

  ee->FixParameter(2, 0);
  hfitb->Fit("ee", "le", "");

  
  /////////////////////////////////////////////////////////
  th.Draw("dt/1000 >> hdispb(50, 0.001, 0.501)", cut, "e");
  tg.Draw("dt/1000 >> +hdispb", cut, "e");

  for(int i = 0; i <= 4; i+=2){
    ee    ->SetParameter(i, ee    ->GetParameter(i) * hdispb->GetBinWidth(1)/
                                                      hfitb->GetBinWidth(1));
    eewith->SetParameter(i, eewith->GetParameter(i) * hdispb->GetBinWidth(1)/
                                                      hfitb->GetBinWidth(1));
  }

  ee->Draw("same");
  eewith->Draw("same");


  ///////////////////////////////////////////////////////////
  const double Geff = 0.9 /* Made up! */,
               Heff = 0.77 /* doc-5787, but maybe misinterpreted */;

  const double Ncorr = Nrawg/Geff + Nrawh/Heff;
  const double Ncorrup = sqrt((Nrawgup/Geff)**2 + (Nrawhup/Heff)**2);
  const double Ncorrlo = sqrt((Nrawglo/Geff)**2 + (Nrawhlo/Heff)**2);
  const double Ncorruph = sqrt((Nrawguph/Geff)**2 + (Nrawhuph/Heff)**2);
  const double Ncorrloh = sqrt((Nrawgloh/Geff)**2 + (Nrawhloh/Heff)**2);

  printf("%s\n--> Eff corr'd sum using MINOS summed quadrature errors:\n", RED);
  printf("Both with eff = %.2f +%f %f%s\n", Ncorr, Ncorrup, Ncorrlo, CLR);


  printf("%s\n--> Eff corr'd sum using HESSE summed quadrature errors:\n", RED);
  printf("Both with eff = %.2f +%f %f%s\n", Ncorr, Ncorruph, Ncorrloh, CLR);

  const double toprob = 1.0/ (3.6*489.509)/ betanprob / 0.852 / 0.981; 

  printf("Both with eff prob = %.2g +%f %f%s\n", Ncorr*toprob, Ncorruph*toprob, Ncorrloh*toprob, CLR);

  TF1 gaus("gaus", "gaus(0)", 0, 20);
  gaus.SetParameters(1, Ncorr, Ncorruph);

  for(int i = 1; i < 400; i++){
    const double up = 0.1*i;
    const double frac = gaus.Integral(0, up)/gaus.Integral(0, 20);
    if(frac > 0.9){
      printf("90%% limit = %f\n", up);
      printf("90%% limit prob = %f\n", up*toprob);
      printf("90%% limit prob *0.1/1.22 = %f\n", up*toprob *(1+0.1/1.22));
      break;
    }
  }

}
