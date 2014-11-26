static void scalemarker(TMarker * m)
{
  const double newx = m->GetX() * c15life/c15eff;
  const double newy = m->GetY() * be11life/be11eff;
  m->SetX(newx);
  m->SetY(newy);
}

static void scalegraph(TGraph * g)
{
  for(int i = 0; i < g->GetN(); i++){
    const double newx = g->GetX()[i] * c15life/c15eff;
    const double newy = g->GetY()[i] * be11life/be11eff;
    g->SetPoint(i, newx, newy);
  }
}

void li9finalfit(bool neutron = false)
{
  const double li9ebn = 0.5080,
               he8ebn = 0.1600;

  TCanvas * c1 = new TCanvas;
  const char * const RED     = "\033[31;1m"; // bold red
  const char * const CLR      = "\033[m"    ; // clear

  TTree tg("t", "t");
  TTree th("t", "t");
  tg.ReadFile("li9-20140925-earlymich.Gd.ntuple");
  th.ReadFile("li9-20141119.H.ntuple");

  tg.Draw("dt/1000 >> hfitg(10000, 0.001, 10)", Form("%sdist < 300  && miche < 12 && !earlymich", neutron?"n==1&&":" "));
  tg.Draw("dt/1000 >> hdispg(50, 0.001, 10.001)", Form("%sdist < 300  && miche < 12 && !earlymich", neutron?"n==1&&":" "));

  TF1 ee("ee", "[0]*exp(-x*log(2)/[1]) + [2]*exp(-x*log(2)/[3]) + [4]", 0, 10);
  ee.FixParameter(1, 0.1783);
  ee.FixParameter(3, 0.1191);
  ee.FixParameter(2, 0);

  hfitg->Fit("ee", "le", "");

  TF1 e("e", "[0]*exp(-x*log(2)/[1])", 0, 10);
  for(int i = 0; i < 2; i++)
    e.SetParameter(i, ee.GetParameter(i));

  const double Nrawg = e.Integral(0, 10)/hfitg->GetBinWidth(1);
  const double Nrawgup = e.Integral(0, 10)/hfitg->GetBinWidth(1) / ee.GetParameter(0) * gMinuit->fErp[0];
  const double Nrawglo = e.Integral(0, 10)/hfitg->GetBinWidth(1) / ee.GetParameter(0) * gMinuit->fErn[0];

  printf("%sgd Nraw = %.2f +%f %f%s\n", RED, Nrawg, Nrawgup, Nrawglo, CLR);

  for(int i = 0; i <= 4; i+=2)
    ee.SetParameter(i, ee.GetParameter(i)*hdispg->GetBinWidth(1)/
                                           hfitg->GetBinWidth(1));

  hdispg->Draw("e");

  th.Draw("dt/1000 >> hfitg(10000, 0.001, 10)", Form("%sdist < 300  && miche < 12 && !earlymich", neutron?"n==1&&":" "));
  th.Draw("dt/1000 >> hdispg(50, 0.001, 10.001)", Form("%sdist < 300  && miche < 12 && !earlymich", neutron?"n==1&&":" "));

  TF1 ee("ee", "[0]*exp(-x*log(2)/[1]) + [2]*exp(-x*log(2)/[3]) + [4]", 0, 10);
  ee.FixParameter(1, 0.1783);
  ee.FixParameter(3, 0.1191);
  ee.FixParameter(2, 0);
  ee.SetParLimits(0, 0, 10);

  hfitg->Fit("ee", "le", "");

  TF1 e("e", "[0]*exp(-x*log(2)/[1])", 0, 10);
  for(int i = 0; i < 2; i++)
    e.SetParameter(i, ee.GetParameter(i));
 
  const double Nrawh = e.Integral(0, 10)/hfitg->GetBinWidth(1);
  const double Nrawhup = e.Integral(0, 10)/hfitg->GetBinWidth(1) / ee.GetParameter(0) * gMinuit->fErp[0];
  const double Nrawhlo = e.Integral(0, 10)/hfitg->GetBinWidth(1) / ee.GetParameter(0) * gMinuit->fErn[0];

  printf("%sH Nraw = %.2f +%f %f%s\n", RED, Nrawh , Nrawhup , Nrawhlo, CLR);

  for(int i = 0; i <= 4; i+=2)
    ee.SetParameter(i, ee.GetParameter(i)*hdispg->GetBinWidth(1)/
                                           hfitg->GetBinWidth(1));

  th.Draw("dt/1000 >> hfitb(10000, 0.001, 10)", Form("%sdist < 300  && miche < 12 && !earlymich", neutron?"n==1&&":" "));
  th.Draw("dt/1000 >> hdispb(40, 0, 10)", Form("%sdist < 300  && miche < 12 && !earlymich", neutron?"n==1&&":" "));
  tg.Draw("dt/1000 >> +hfitb", Form("%sdist < 300  && miche < 12 && !earlymich", neutron?"n==1&&":" "));
  tg.Draw("dt/1000 >> +hdispb", Form("%sdist < 300  && miche < 12 && !earlymich", neutron?"n==1&&":" "));

  TH1D * hdispb = (TH1D*) gROOT->FindObject("hdispb");
  TH1D * hfitb = (TH1D*) gROOT->FindObject("hfitb");

  TF1 ee("ee", "[0]*exp(-x*log(2)/[1]) + [2]*exp(-x*log(2)/[3]) + [4]", 0, 10);
  ee.FixParameter(1, 0.1783);
  ee.FixParameter(3, 0.1191);
  ee.FixParameter(2, 0);

  hfitb->Fit("ee", "le", "");

  TF1 e("e", "[0]*exp(-x*log(2)/[1])", 0, 10);
  for(int i = 0; i < 2; i++)
    e.SetParameter(i, ee.GetParameter(i));

  // First factor takes into account the prompt energy cut, second
  // as documented 
  const double Geff = (neutron?0.97*0.64:1)* // neutron eff
               0.996*0.9217, /* DC3rdPub product of muon, light noise, OV, multiplicity, neutron (E, t, R), FV and IV efficiencies */
               Heff = (neutron?0.97*0.93:1)* // neutron eff
               0.993*0.870 /* doc-5787 and others, but maybe misinterpreted */;

  const double Ncorr = Nrawg/Geff + Nrawh/Heff;
  const double Ncorrup = sqrt((Nrawgup/Geff)**2 + (Nrawhup/Heff)**2);
  const double Ncorrlo = sqrt((Nrawglo/Geff)**2 + (Nrawhlo/Heff)**2);

  printf("%s\n--> Sum using summed quadrature errors:\n", RED);
  printf("Both with eff = %.2f +%f %f%s\n", Ncorr , Ncorrup , Ncorrlo, CLR);

  TF1 * eedisp = new TF1("eedisp", "[0]*exp(-x*log(2)/[1]) + [2]*exp(-x*log(2)/[3]) + [4]", 0, 10);
  eedisp->FixParameter(1, 0.1783);
  eedisp->FixParameter(3, 0.1191);
  eedisp->FixParameter(2, 0);


  TF1 ee2("ee2",
    Form("(([0]-[2])*%f*%f*exp(-    x    *log(2)/0.1783)+[1])*(x < 10) +"
         "(   [2]   *%f*%f*exp(-(x-9.999)*log(2)/0.1783)+[3])*(x >=10)",
         Heff, li9ebn, Geff, li9ebn), 0, 20);
  ee2.SetParNames("totalli9", "hbg", "gdli9", "gdbg");
  ee2.SetParameters(1, 1, 1, 1);
  ee2.SetNpx(400);
  ee2.SetLineColor(kRed);
  ee2.SetParLimits(0, 0, 20);
  ee2.SetParLimits(2, 0, 50);
  th.Draw("dt/1000 >> h2fit(1000, 0.001, 19.998)", Form("%sdist < 300  && miche < 12 && !earlymich", neutron?"n==1&&":" "));
  tg.Draw("dt/1000+9.999 >> +h2fit", Form("%sdist < 300  && miche < 12 && !earlymich", neutron?"n==1&&":" "));
  h2fit->Fit("ee2", "l", "", 0, 20);

  const double denominator = 489.509*367*0.852;

  //////////////////
  new TCanvas;

  TF1 ee2str("ee2str",
    Form("%f*%f*(%f*(%f*([0]-[2])/0.1783*log(2)*exp(-    x    *log(2)/0.1783) + %f*([5]-[4])/0.1191*log(2)*exp(-    x    *log(2)/0.1191) + [1])*(x < 10) +"
             "%f*(%f*   [2]   /0.1783*log(2)*exp(-(x-9.999)*log(2)/0.1783) + %f*   [4]   /0.1191*log(2)*exp(-(x-9.999)*log(2)/0.1191) + [3])*(x >=10))",
         h2fit->GetBinWidth(1), denominator, Heff, li9ebn, he8ebn, Geff, li9ebn, he8ebn), 0, 20);
  ee2str.SetParameters(ee2.GetParameter(0)/denominator,
                       ee2.GetParameter(1)/denominator, 
                       ee2.GetParameter(2)/denominator,
                       ee2.GetParameter(3)/denominator, 1, 1);
  ee2str.SetParNames("totalli9" /* c0, f1 */, "hbg", "gdli9", "gdbg", "gdhe8", "totalhe8");
  ee2str.SetNpx(400);
  ee2str.SetLineColor(kRed);
  ee2str.SetParLimits(0, 0, 120/denominator);
  ee2str.SetParLimits(2, 0, 120/denominator);
  ee2str.SetParLimits(4, 0, 300/denominator);
  ee2str.SetParLimits(5, 0, 300/denominator);

  h2fit->Fit("ee2str", "le", "", 0, 20);

  gMinuit->Command("set print 0");
  gMinuit->Command("set strategy 2");
  gMinuit->fUp = 2.3/2; // 90% in 1D
  gMinuit->Command("mncont 1 6 200");
  TGraph * ninty_1d = (TGraph*)((TGraph*)gMinuit->GetPlot())->Clone();
  gMinuit->fUp = 4.61/2; // 90%
  gMinuit->Command("mncont 1 6 200");
  TGraph * ninty_2d = (TGraph*)((TGraph*)gMinuit->GetPlot())->Clone();
  gMinuit->fUp = 11.83/2; // 99.73%
  gMinuit->Command("mncont 1 6 200");
  TGraph * ninty973_2d = (TGraph*)((TGraph*)gMinuit->GetPlot())->Clone();
  if(ninty973_2d) ninty973_2d->SetFillColor(kViolet);
  if(ninty973_2d) ninty973_2d->Draw("alf");
  if(ninty_2d) ninty_2d->SetFillColor(kBlue);
  if(ninty_2d) ninty_2d->Draw("lf");
  if(ninty_1d) ninty_1d->SetLineColor(kRed);
  if(ninty_1d) ninty_1d->Draw("l");

  if(ninty973_2d) ninty973_2d->GetYaxis()->SetRangeUser(0, 300/denominator);
  if(ninty973_2d) ninty973_2d->GetXaxis()->SetRangeUser(0, 120/denominator);


  //////////////////
 
  c1->cd();

  h2fit->Fit("ee2", "le", "", 0, 20); // for gMinuit
  e.SetParameter(0, ee2.GetParameter(0));
  e.SetParameter(1, 0.1783);
 
  const double N2raw = e.Integral(0, 10)/h2fit->GetBinWidth(1);
  const double N2rawup = e.Integral(0, 10)/h2fit->GetBinWidth(1) / ee2.GetParameter(0) * gMinuit->fErp[0];
  const double N2rawlo = e.Integral(0, 10)/h2fit->GetBinWidth(1) / ee2.GetParameter(0) * gMinuit->fErn[0];

  printf("%sJoint fit to number of raw selected Li-9 = %.2f +%f %f%s\n",
         RED, N2raw, N2rawup , N2rawlo, CLR);

  for(int i = 0; i <= 4; i+=2)
    eedisp->SetParameter(i, ee.GetParameter(i)*hdispg->GetBinWidth(1)/
                                               hfitb->GetBinWidth(1));

  TF1 ee2s("ee2s",
    Form("([0]*%f*%f*exp(-    x    *log(2)/0.1783)+[1])*(x < 10) +"
         "([2]*%f*%f*exp(-(x-9.999)*log(2)/0.1783)+[3])*(x >=10)",
         Heff, li9ebn, Geff, li9ebn), 0, 20);
  ee2s.SetParameters(1, 1, 1, 1);
  ee2s.SetLineColor(kRed);
  ee2s.SetParLimits(0, 0, 10);
  ee2s.SetParLimits(2, 0, 10);
  th.Draw("dt/1000 >> h2fit(1000, 0.001, 19.998)", Form("%sdist < 300  && miche < 12 && !earlymich", neutron?"n==1&&":" "));
  tg.Draw("dt/1000+9.999 >> +h2fit", Form("%sdist < 300  && miche < 12 && !earlymich", neutron?"n==1&&":" "));
  h2fit->Fit("ee2s", "l", "", 0, 20);

  const double bestlike = gMinuit->fAmin;

  ee2s.FixParameter(0, 0);
  ee2s.FixParameter(2, 0);
  h2fit->Fit("ee2s", "l");

  const double nulllike = gMinuit->fAmin;

  printf("%ssignificance = %f\n%s", RED, sqrt(2)*sqrt(nulllike - bestlike), CLR);

  TF1 ee("ee", "[0]*exp(-x*log(2)/[1]) + [2]*exp(-x*log(2)/[3]) + [4]", 0, 10);
  ee.SetParameter(0, 1);
  ee.SetParameter(1, 0.1783);
  ee.FixParameter(3, 0.1191);
  ee.FixParameter(2, 0);

  hfitb->Fit("ee", "le");

  printf("\n%s--> li9 half-life = %f +%f %f%s\n", RED,
          ee.GetParameter(1), 
          gMinuit->fErp[1], gMinuit->fErn[1], CLR);

  // Just for looks:
  hdispb->Draw("e");
  eedisp->SetLineColor(kRed);
  eedisp->SetNpx(400);
  eedisp->Draw("Same");

  printf("Candidates:\n");
  th.SetScanField(0);
  tg.SetScanField(0);
  th.Scan("run:trig:dt", Form("%sdist < 300  && miche < 12 && !earlymich && dt < 10000", neutron?"n==1&&":" "));
  tg.Scan("run:trig:dt", Form("%sdist < 300  && miche < 12 && !earlymich && dt < 10000", neutron?"n==1&&":" "));
}
