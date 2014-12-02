void li9finalfit(bool neutron = false)
{
  const char * const RED     = "\033[31;1m"; // bold red
  const char * const CLR      = "\033[m"    ; // clear

  ///////////////////////////////////////////////////////////////////
 
  const float dist = 300;

  const double livetime = 489.509, nstop = 367.;
  const double nstoptarg = 139.0;

  const double li9ebn = 0.5080, he8ebn = 0.1600, n17ebn = 0.951;
  const double distcuteff = (dist == 400?0.948:dist == 300?0.852:
                             dist == 200?0.565:dist==159?0.376:100000);

  // 30s begin-of-run requirement taken into account here
  const double denominator = 0.99127*livetime*nstop*distcuteff;

  const double hcapfrac = 0.871;

  /* DC3rdPub product of muon, light noise, OV, multiplicity, neutron (E, t, R), FV and IV efficiencies */
  const double Geff_sans_prompt_or_mun = (1-4.49/100.)*
         (1-0.01/100.)*
         (1-0.06/100.)*
         (1-1.06/100.)*
         0.9829*
         (1-0.66/100.)*
         (1-0.04/100.);

  const double Heff_sans_prompt_or_mun =   (1-1.25*4.49/100.)* // muon - ok, straightforwards scaling
         (1-0.01/100.)* // ? LN - same cut, but not obviously same eff
                        // however, *very* small for Gd, so...
         (1-0.06/100.)* // OV - ok, same
         (1-((8.+9.)/(2.+6.))*1.06/100.)* // Mult - I think this is valid
         0.9512*       // neutron (ANN,E,t,R) -- plot on slide 6 of 
                      // doc-5863 -- I think this is right
         (1-0.806/100.)* // FV -- doc5480
         (1-0.025/100.)* // IV prompt -- doc5813
         (1-0.014/100.); // IV delayed -- doc5813


  // First factor takes into account the efficiency of selecting a
  // neutron after a muon, second is the prompt energy cut, third as
  // documented above
  const double Geff = (neutron?0.97*0.64:1)*0.996*Geff_sans_prompt_or_mun,
               Heff = (neutron?0.97*0.93:1)*0.993*Heff_sans_prompt_or_mun;

  char cut[1000];
  snprintf(cut, 999, "%sdist < %f  && miche < 12 && !earlymich", dist, neutron?"n==1&&":" ");

  ////////////////////////////////////////////////////////////////////

  TCanvas * c2 = new TCanvas("c2", "c2", 600, 350);
  TTree tg("t", "t");
  TTree th("t", "t");
  tg.ReadFile("li9-20141202.Gd.ntuple");
  th.ReadFile("li9-20141202.H.ntuple");
  th.Draw("dt/1000 >> h2fit(1000, 0.001, 59.998)", cut);
  tg.Draw("dt/1000+29.999 >> +h2fit", cut, "e");

  //////////////////

  TF1 ee2str("ee2str",
    Form("%f*%f*"
    "(%f*(%f*[2]*(1-[1])/0.1783*log(2)*exp(-     x    *log(2)/0.1783) +" // H Li-9
         "%f*[3]*(1-[1])/0.1191*log(2)*exp(-     x    *log(2)/0.1191) +" // H He-8
         "%f*[4]*(1-[5])/4.173 *log(2)*exp(-     x    *log(2)/4.173 ) +" // H N-17
          "[0]*(1-[1]))*(x < 30) + "                                    // H bg
     "%f*(%f*  [2]*[1]  /0.1783*log(2)*exp(-(x-29.999)*log(2)/0.1783) +" // Gd Li-9
         "%f*  [3]*[1]  /0.1191*log(2)*exp(-(x-29.999)*log(2)/0.1191) +" // Gd He-8
         "%f*  [4]*[5]  /4.173 *log(2)*exp(-(x-29.999)*log(2)/4.173 ) +" // Gd N-17
          "[0]*[1])*(x >=30))",                                         // Gd bg
         h2fit->GetBinWidth(1), denominator, Heff, li9ebn, he8ebn, n17ebn, Geff, li9ebn, he8ebn, n17ebn), 0, 60);
  ee2str.SetParameters(0.0001, 0.38, 0.0001, 0.00001, 0.00001, 0.1);
  ee2str.SetParNames(  "bg",  "gdfrac", "li9",  "he8", "n17", "n17gdfrac");
  ee2str.SetNpx(400);

  ee2str.SetParLimits(2, 0, 1e-3);
//  ee2str.FixParameter(1,  165799./7778371 * 8.32142e-05); // XXX
  ee2str.SetParLimits(1, 0, 1);
  ee2str.SetParLimits(5, 0, 1);
  ee2str.FixParameter(1, hcapfrac*nstoptarg/(nstop-hcapfrac*nstoptarg) * Geff/Heff);
  ee2str.SetParLimits(3, 0, 1e-3);
  ee2str.SetParLimits(4, 0, 1e-2);

  ee2str.FixParameter(3, 0);
  h2fit->Fit("ee2str", "le", "", 0, 60);
  printf("%sLi-9 prob without He-8: %f +%f %f%s\n", RED,
    ee2str.GetParameter(2), gMinuit->fErp[2], gMinuit->fErn[2], CLR);

  ee2str.ReleaseParameter(3);
  ee2str.SetParLimits(3, 0, 1e-3);
  h2fit->Fit("ee2str", "le", "", 0, 60);
  printf("%sLi-9 prob with He-8: %f +%f %f%s\n", RED,
    ee2str.GetParameter(2), gMinuit->fErp[2], gMinuit->fErn[2], CLR);

  TF1 * ee2str_save = (TF1 *)ee2str.Clone("ee2str_save");

  double minx = ee2str.GetParameter(2), miny = ee2str.GetParameter(3);

  TCanvas * c3 = new TCanvas("c3", "c3", 600, 350);
  gMinuit->Command("Set print 0");
  gMinuit->Command("Set strategy 2");

  gMinuit->fUp = 1.0/2; // 68% in 1D
  gMinuit->Command("mncont 3 4 100");
  TGraph * sigma_1d = gMinuit->GetPlot()?(TGraph*)((TGraph*)gMinuit->GetPlot())->Clone():NULL;

  gMinuit->fUp = 2.3/2; // 90% in 1D
//  gMinuit->Command("mncont 3 4 100");
  TGraph * ninty_1d = gMinuit->GetPlot()?(TGraph*)((TGraph*)gMinuit->GetPlot())->Clone():NULL;

  gMinuit->fUp = 4.61/2; // 90%
//  gMinuit->Command("mncont 3 4 100");
  TGraph * ninty_2d = gMinuit->GetPlot()?(TGraph*)((TGraph*)gMinuit->GetPlot())->Clone():NULL;

  gMinuit->fUp = 11.83/2; // 99.73%
  gMinuit->Command("mncont 3 4 100");
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
  c1->Divide(1, 3);
  c1->cd(1);

  tg.Draw("dt/1000 >> hdispg(120, 0, 30)", cut);
  th.Draw("dt/1000 >> hdisph(120, 0, 30)", cut);

  th.Draw("dt/1000 >> hdisp(120, 0, 30)", cut);
  tg.Draw("dt/1000 >> +hdisp", cut, "e");

  TF1 * eedisp = new TF1("eedisp",
    Form("%f*%f*"
    "([6]*%f*(%f*[2]*(1-[1])/0.1783*log(2)*exp(-x*log(2)/0.1783) +"// H Li-9
         "%f*[3]*(1-[1])/0.1191*log(2)*exp(-x*log(2)/0.1191) +"    // H He-8
         "%f*[4]*(1-[5])/4.173 *log(2)*exp(-x*log(2)/4.173 ) +"    // H N-17
          "[0]*(1-[1])) + "                                        // H bg
     "[7]*%f*(%f*  [2]*[1]  /0.1783*log(2)*exp(-x*log(2)/0.1783) +"// Gd Li-9
         "%f*  [3]*[1]  /0.1191*log(2)*exp(-x*log(2)/0.1191) +"    // Gd He-8
         "%f*  [4]*[5]  /4.173 *log(2)*exp(-x*log(2)/4.173 ) +"    // Gd N-17
          "[0]*[1]))",                                             // Gd bg
         h2fit->GetBinWidth(1), denominator, Heff, li9ebn, he8ebn,
         n17ebn, Geff, li9ebn, he8ebn, n17ebn), 0, 30);
  for(int i = 0; i < eedisp->GetNpar(); i++)
    eedisp->SetParameter(i, ee2str_save->GetParameter(i));

  eedisp->SetParameter(6, hdisp->GetBinWidth(1)/h2fit->GetBinWidth(1));
  eedisp->SetParameter(7, hdisp->GetBinWidth(1)/h2fit->GetBinWidth(1));
  eedisp->SetLineColor(kRed);
  eedisp->SetNpx(400);
  eedisp->Draw("Same");

  c1->cd(2);
  TF1 * eedisph = (TF1*)eedisp->Clone("eedisph");
  eedisph->SetParameter(6, hdisp->GetBinWidth(1)/h2fit->GetBinWidth(1));
  eedisph->SetParameter(7, 0);
  hdisph->Draw("e");
  eedisph->Draw("same");
  
  c1->cd(3);
  TF1 * eedispg = (TF1*)eedisp->Clone("eedispg");
  eedispg->SetParameter(6, 0);
  eedispg->SetParameter(7, hdisp->GetBinWidth(1)/h2fit->GetBinWidth(1));
  hdispg->Draw("e");
  eedispg->Draw("same");

return;
  printf("Candidates:\n");
  th.SetScanField(0);
  tg.SetScanField(0);
  th.Scan("run:trig:dt:sqrt(dx*dx+dy*dy):dx:dy:dz", cut);
  tg.Scan("run:trig:dt:sqrt(dx*dx+dy*dy):dx:dy:dz", cut);
}
