void li9finalfit(bool neutron = false)
{
  const char * const RED     = "\033[31;1m"; // bold red
  const char * const CLR      = "\033[m"    ; // clear

  ///////////////////////////////////////////////////////////////////
 
  const float dist = 300;

  const double livetime = 489.509, nstop = 367.;

  const double li9ebn = 0.5080, he8ebn = 0.1600;
  const double distcuteff = (dist == 400?0.948:dist == 300?0.852:
                             dist == 200?0.565:dist==159?0.376:100000);

  // 30s begin-of-run requirement taken into account here
  const double denominator = 0.99127*livetime*nstop*distcuteff;

  // First factor takes into account the prompt energy cut, second
  // as documented 
  const double Geff = (neutron?0.97*0.64:1)* // neutron eff
               0.996*0.9217, /* DC3rdPub product of muon, light noise, OV, multiplicity, neutron (E, t, R), FV and IV efficiencies */
               Heff = (neutron?0.97*0.93:1)* // neutron eff
               0.993*0.870 /* doc-5787 and others, but maybe misinterpreted */;

  char cut[1000];
  snprintf(cut, 999, "%sdist < %f  && miche < 12 && !earlymich", dist, neutron?"n==1&&":" ");

  ////////////////////////////////////////////////////////////////////

  TCanvas * c2 = new TCanvas("c2", "c2", 600, 350);
  TTree tg("t", "t");
  TTree th("t", "t");
  tg.ReadFile("li9-20140925-earlymich.Gd.ntuple");
  th.ReadFile("li9-20141119.H.ntuple");
  th.Draw("dt/1000 >> h2fit(1000, 0.001, 19.998)", cut);
  tg.Draw("dt/1000+9.999 >> +h2fit", cut);

  //////////////////

  TF1 ee2str("ee2str",
    Form("%f*%f*"
    "(%f*(%f*([2]*(1-[1]))/0.1783*log(2)*exp(-    x    *log(2)/0.1783) +" // H Li-9
         "%f*([3]*(1-[1]))/0.1191*log(2)*exp(-    x    *log(2)/0.1191) +" // H He-8
          "[0]*(1-[1]))*(x < 10) + "                                      // H bg
     "%f*(%f*  [2]*[1]    /0.1783*log(2)*exp(-(x-9.999)*log(2)/0.1783) +" // Gd Li-9
         "%f*  [3]*[1]    /0.1191*log(2)*exp(-(x-9.999)*log(2)/0.1191) +" // Gd He-8
          "[0]*[1])*(x >=10))",                                           // Gd bg
         h2fit->GetBinWidth(1), denominator, Heff, li9ebn, he8ebn, Geff, li9ebn, he8ebn), 0, 20);
  ee2str.SetParameters(0.0001, 0.38, 0.0001, 0.00001);
  ee2str.SetParNames(  "bg",  "gdfrac", "li9",  "he8");

  ee2str.SetParLimits(2, 0, 1e-3);
//  ee2str.FixParameter(1,  165799./7778371 * 8.32142e-05); // XXX
  ee2str.SetParLimits(1, 0, 1);
  ee2str.FixParameter(1, 0.87 /* H-n paper */ *139./367 * Geff/Heff);
  ee2str.SetParLimits(3, 0, 1e-3);

  ee2str.FixParameter(3, 0);
  h2fit->Fit("ee2str", "le", "", 0, 20);
  printf("%sLi-9 prob without He-8: %f +%f %f%s\n", RED,
    ee2str.GetParameter(2), gMinuit->fErp[2], gMinuit->fErn[2], CLR);

  ee2str.ReleaseParameter(3);
  ee2str.SetParLimits(3, 0, 1e-3);
  h2fit->Fit("ee2str", "le", "", 0, 20);
  printf("%sLi-9 prob with He-8: %f +%f %f%s\n", RED,
    ee2str.GetParameter(2), gMinuit->fErp[2], gMinuit->fErn[2], CLR);

  const double li9prob = ee2str.GetParameter(2),
               he8prob = ee2str.GetParameter(3),
               gdfrac = ee2str.GetParameter(1);
  double       bg = ee2str.GetParameter(0)* denominator;

  double minx = ee2str.GetParameter(2), miny = ee2str.GetParameter(3);

  gMinuit->Command("Set print 0");
  gMinuit->Command("Set strategy 2");

  gMinuit->fUp = 1.0/2; // 68% in 1D
  gMinuit->Command("mncont 3 4 200");
  TGraph * sigma_1d = gMinuit->GetPlot()?(TGraph*)((TGraph*)gMinuit->GetPlot())->Clone():NULL;

  gMinuit->fUp = 2.3/2; // 90% in 1D
  gMinuit->Command("mncont 3 4 200");
  TGraph * ninty_1d = gMinuit->GetPlot()?(TGraph*)((TGraph*)gMinuit->GetPlot())->Clone():NULL;

  gMinuit->fUp = 4.61/2; // 90%
  gMinuit->Command("mncont 3 4 200");
  TGraph * ninty_2d = gMinuit->GetPlot()?(TGraph*)((TGraph*)gMinuit->GetPlot())->Clone():NULL;

  gMinuit->fUp = 11.83/2; // 99.73%
  gMinuit->Command("mncont 3 4 200");
  TGraph * ninty973_2d = gMinuit->GetPlot()?(TGraph*)((TGraph*)gMinuit->GetPlot())->Clone():NULL;

  if(ninty973_2d) ninty973_2d->SetFillColor(kViolet);
  if(ninty973_2d) ninty973_2d->Draw("al");
  if(ninty973_2d) ninty973_2d->GetXaxis()->SetRangeUser(0, 0.00075);
  if(ninty973_2d) ninty973_2d->GetYaxis()->SetRangeUser(0, 0.001);
  if(ninty973_2d) ((TGaxis*)(ninty973_2d->GetXaxis()))->SetMaxDigits(2);
  if(ninty973_2d) ((TGaxis*)(ninty973_2d->GetYaxis()))->SetMaxDigits(2);
  if(ninty_2d) ninty_2d->SetFillColor(kBlue);
  if(ninty_2d) ninty_2d->Draw("l");
  if(ninty_1d) ninty_1d->SetLineColor(kRed);
  if(ninty_1d) ninty_1d->Draw("l");
  if(sigma_1d) sigma_1d->SetLineColor(kBlack);
  if(sigma_1d) sigma_1d->Draw("l");

  TMarker * best = new TMarker(minx, miny, kStar);
  best->Draw();
  
  //////////////////////////////////////////////////////////////////////
 
  TCanvas * c1 = new TCanvas;

  th.Draw("dt/1000 >> hdisp(60, 0, 30)", cut);
  tg.Draw("dt/1000 >> +hdisp", cut, "e");

  TF1 * eedisp = new TF1("eedisp",
    "[0]*exp(-x*log(2)/0.1783)+[1]*exp(-x*log(2)/0.1191)+[2]", 0, 30);
  eedisp->SetLineColor(kRed);
  eedisp->SetNpx(400);

  eedisp->SetParameter(0, li9prob*denominator*li9ebn*(Heff*(1-gdfrac)+Geff*gdfrac)/0.1783*log(2));
  eedisp->SetParameter(1, he8prob*denominator*he8ebn*(Heff*(1-gdfrac)+Geff*gdfrac)/0.1191*log(2));
  eedisp->SetParameter(2, bg * hdisp->GetBinWidth(1));

  // Just for looks:
  eedisp->Draw("Same");

  printf("Candidates:\n");
  th.SetScanField(0);
  tg.SetScanField(0);
  th.Scan("run:trig:dt:sqrt(dx*dx+dy*dy):dx:dy:dz", cut);
  tg.Scan("run:trig:dt:sqrt(dx*dx+dy*dy):dx:dy:dz", cut);
}
