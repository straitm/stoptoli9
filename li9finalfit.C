void setupee2str(TF1 & ee2str, const double expectedgdfrac)
{
  ee2str.SetParLimits(0, 0, 1e-3);
  ee2str.SetParLimits(1, 0, 1);
  ee2str.FixParameter(1, expectedgdfrac);
  ee2str.SetParLimits(2, 0, 1e-3);
  ee2str.SetParLimits(3, 0, 1e-3);
  ee2str.SetParLimits(4, 0, 20); // allow 100% production, plus
                                // a factor of two for the O capture
                                // normalization, plus a factor of ten
                                // for other complications.
  ee2str.SetParLimits(5, 0, 0.5);
  ee2str.SetParLimits(6, 0, 1);
  ee2str.SetParLimits(7, 0, 1);
  ee2str.SetParLimits(8, 0, 1);
}

void li9finalfit(bool neutron = false)
{
  const char * const RED     = "\033[31;1m"; // bold red
  const char * const CLR      = "\033[m"    ; // clear

  ///////////////////////////////////////////////////////////////////
 
  const float dist = 300;

  const double livetime = 489.509, nstop = 367.;
  const double nstoptarg = 139.0;

  // bn decay probability multiplied by the number of captures relative
  // to C-12
  const double li9ebn = 0.5080,
               he8ebn = 0.1600,
               c16ebn = 0.99  * 7.0/nstop,
               n17ebn = 0.951 * 7.0/nstop,
               b13ebn = 0.0029 * 3.7/nstop,
               li11ebn = 0.789 * 3.7/nstop;
               
  const double distcuteff = (dist == 400?0.948:dist == 300?0.852:
                             dist == 200?0.565:dist==159?0.376:100000);

  // 30s begin-of-run requirement taken into account here
  const double denominator = 0.99127*livetime*nstop*distcuteff;

  const double gdcapfrac = 0.871;

  /* DC3rdPub product of muon, light noise, OV, multiplicity,
     neutron (E, t, R), FV and IV efficiencies */
  const double Geff_sans_prompt_or_mun = (1-4.49/100.)*
         (1-0.01/100.)*
         (1-0.06/100.)*
         (1-1.06/100.)*
         0.9829*
         (1-0.66/100.)*
         (1-0.04/100.);

  const double Heff_sans_prompt_or_mun =
         (1-1.25*4.49/100.)* // muon - ok, straightforwards scaling
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
  const double Geff=(neutron?0.97*0.64:1)*0.996*Geff_sans_prompt_or_mun,
               Heff=(neutron?0.97*0.93:1)*0.993*Heff_sans_prompt_or_mun;

  const double expectedgdfrac = gdcapfrac*nstoptarg/nstop;

  char cut[1000];
  snprintf(cut, 999, "%sdist < %f  && miche < 12 && !earlymich && dt < 30000",
           dist, neutron?"n==1&&":" ");

  ////////////////////////////////////////////////////////////////////

  TCanvas * c2 = new TCanvas("c2", "c2", 600, 350);
  TTree tg("t", "t");
  TTree th("t", "t");
  tg.ReadFile("li9-20141202.Gd.ntuple");
  th.ReadFile("li9-20141202.H.ntuple");
  th.Draw("dt/1000 >> h2fit(1000, 0.001, 59.999)", cut);
  tg.Draw("dt/1000+29.999 >> +h2fit", cut, "e");

  //////////////////

  TF1 ee2str("ee2str", Form("%f*"
    "(%f*("
       "%f*[2]*(1-[1])/0.257233*exp(-x/0.257233)+"  // H Li-9
       "%f*[3]*(1-[1])/0.171825*exp(-x/0.171825)+"  // H He-8
       "%f*[4]*(1-[5])/6.020366*exp(-x/6.020366)+"  // H N-17
       "%f*[6]*(1-[5])/1.077693*exp(-x/1.077693)+"  // H C-16
       "%f*[7]*(1-[1])/0.025002*exp(-x/0.025002)+"  // H B-13
       "%f*[8]*(1-[1])/0.012624*exp(-x/0.012624)+"  // H Li-11
        "[0]*(1-[1]))*(x < 30) + "                  // H bg
     "%f*("
       "%f*[2]*[1]/0.257233*exp(-(x-29.999)*(x >=30)/0.257233)+" // Gd Li-9
       "%f*[3]*[1]/0.171825*exp(-(x-29.999)*(x >=30)/0.171825)+" // Gd He-8
       "%f*[4]*[5]/6.020366*exp(-(x-29.999)*(x >=30)/6.020366)+" // Gd N-17
       "%f*[6]*[5]/1.077693*exp(-(x-29.999)*(x >=30)/1.077693)+" // Gd C-16
       "%f*[7]*[1]/0.025002*exp(-(x-29.999)*(x >=30)/0.025002)+" // Gd B-13
       "%f*[8]*[1]/0.012624*exp(-(x-29.999)*(x >=30)/0.012624)+" // Gd Li-11
        "[0]*[1])*(x >=30))",                           // Gd bg
     h2fit->GetBinWidth(1)*denominator,
     Heff, li9ebn, he8ebn, n17ebn, c16ebn, b13ebn, li11ebn,
     Geff, li9ebn, he8ebn, n17ebn, c16ebn, b13ebn, li11ebn
     ), 0, 60);
  ee2str.SetParameters(1e-4, 0.38, 1e-4, 1e-5, 0.1,  0.1, 0.01, 0.01, 0.01);
  ee2str.SetParNames("bg", "gdfrac","li9", "he8","n17","ocapgdfrac",
                     "c16", "b13", "li11");
  ee2str.SetNpx(400);

  setupee2str(ee2str, expectedgdfrac);

  const unsigned int nfixforpureli9 = 6;
  const double fixforpureli9[nfixforpureli9] = { 3, 4, 5, 6, 7, 8 };

  for(int i = 0; i < nfixforpureli9; i++)
    ee2str.FixParameter(fixforpureli9[i], 0);

  h2fit->Fit("ee2str", "le", "", 0, 60);

  printf("%sLi-9 prob without other isotopes: %f +%f %f%s\n", RED,
    ee2str.GetParameter(2), gMinuit->fErp[1], gMinuit->fErn[1], CLR);

  for(int i = 0; i < nfixforpureli9; i++)
    ee2str.ReleaseParameter(fixforpureli9[i]);
  setupee2str(ee2str, expectedgdfrac);

  h2fit->Fit("ee2str", "le", "", 0, 60);
  printf("%sLi-9 prob with other isotopes: %f +%f %f%s\n", RED,
    ee2str.GetParameter(2), gMinuit->fErp[1], gMinuit->fErn[1], CLR);

  TF1 * ee2str_save = (TF1 *)ee2str.Clone("ee2str_save");

  double minx = ee2str.GetParameter(2), miny = ee2str.GetParameter(3);

  //TCanvas * c3 = new TCanvas("c3", "c3", 600, 350);
  gMinuit->Command("Set print 0");
  gMinuit->Command("Set strategy 2");

  gMinuit->fUp = 1.0/2; // 68% in 1D
  gMinuit->Command("mncont 3 4 20");
  TGraph * sigma_1d =
    gMinuit->GetPlot()?(TGraph*)((TGraph*)gMinuit->GetPlot())->Clone():NULL;

  gMinuit->fUp = 2.3/2; // 90% in 1D
  gMinuit->Command("mncont 3 4 20");
  TGraph * ninty_1d =
    gMinuit->GetPlot()?(TGraph*)((TGraph*)gMinuit->GetPlot())->Clone():NULL;

  gMinuit->fUp = 4.61/2; // 90%
  gMinuit->Command("mncont 3 4 20");
  TGraph * ninty_2d =
    gMinuit->GetPlot()?(TGraph*)((TGraph*)gMinuit->GetPlot())->Clone():NULL;

  if(ninty_2d) ninty_2d->SetFillColor(kViolet);
  if(ninty_2d) ninty_2d->Draw("alf");
  if(ninty_2d) ninty_2d->GetXaxis()->SetRangeUser(0, 0.00079);
  if(ninty_2d) ninty_2d->GetYaxis()->SetRangeUser(0, 0.00149);
  if(ninty_2d) ((TGaxis*)(ninty_2d->GetXaxis()))->SetMaxDigits(2);
  if(ninty_2d) ((TGaxis*)(ninty_2d->GetYaxis()))->SetMaxDigits(2);
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

  tg.Draw("dt/1000 >> hdispg(60, 0, 30)", cut);
  th.Draw("dt/1000 >> hdisph(60, 0, 30)", cut);

  th.Draw("dt/1000 >> hdisp(60, 0, 30)", cut);
  tg.Draw("dt/1000 >> +hdisp", cut, "e");

  TF1 * eedisp = new TF1("eedisp",
    Form("%f*"
    "([9]*%f*("
       "%f*[2]*(1-[1])/0.257233*exp(-x/0.257233)+"  // H Li-9
       "%f*[3]*(1-[1])/0.171825*exp(-x/0.171825)+"  // H He-8
       "%f*[4]*(1-[5])/6.020366*exp(-x/6.020366)+"  // H N-17
       "%f*[6]*(1-[5])/1.077693*exp(-x/1.077693)+"  // H C-16
       "%f*[7]*(1-[1])/0.025002*exp(-x/0.025002)+"  // H B-13
       "%f*[8]*(1-[1])/0.012624*exp(-x/0.012624)+"  // H Li-11
       "[0]*(1-[1])) + "                            // H bg
     "[10]*%f*("
       "%f*[2]*[1]/0.257233*exp(-x/0.257233)+" // Gd Li-9
       "%f*[3]*[1]/0.171825*exp(-x/0.171825)+" // Gd He-8
       "%f*[4]*[5]/6.020366*exp(-x/6.020366)+" // Gd N-17
       "%f*[6]*[5]/1.077693*exp(-x/1.077693)+" // Gd C-16
       "%f*[7]*[1]/0.025002*exp(-x/0.025002)+" // Gd B-13
       "%f*[8]*[1]/0.012624*exp(-x/0.012624)+" // Gd Li-11
       "[0]*[1]))",                                              // Gd bg
         h2fit->GetBinWidth(1)*denominator,
         Heff, li9ebn, he8ebn, n17ebn, c16ebn, b13ebn, li11ebn,
         Geff, li9ebn, he8ebn, n17ebn, c16ebn, b13ebn, li11ebn), 0, 30);
  for(int i = 0; i < eedisp->GetNpar(); i++)
    eedisp->SetParameter(i, ee2str_save->GetParameter(i));

  eedisp->SetParameter(ee2str.GetNpar(),
                       hdisp->GetBinWidth(1)/h2fit->GetBinWidth(1));
  eedisp->SetParameter(ee2str.GetNpar()+1,
                       hdisp->GetBinWidth(1)/h2fit->GetBinWidth(1));
  eedisp->SetLineColor(kRed);
  eedisp->SetNpx(400);
  eedisp->Draw("Same");

  c1->cd(2);
  TF1 * eedisph = (TF1*)eedisp->Clone("eedisph");
  eedisph->SetParameter(ee2str.GetNpar(),
                        hdisp->GetBinWidth(1)/h2fit->GetBinWidth(1));
  eedisph->SetParameter(ee2str.GetNpar()+1, 0);
  hdisph->Draw("e");
  eedisph->Draw("same");
  
  c1->cd(3);
  TF1 * eedispg = (TF1*)eedisp->Clone("eedispg");
  eedispg->SetParameter(ee2str.GetNpar(), 0);
  eedispg->SetParameter(ee2str.GetNpar()+1,
                        hdisp->GetBinWidth(1)/h2fit->GetBinWidth(1));
  hdispg->Draw("e");
  eedispg->Draw("same");

return;
  printf("Candidates:\n");
  th.SetScanField(0);
  tg.SetScanField(0);
  th.Scan("run:trig:dt:sqrt(dx*dx+dy*dy):dx:dy:dz", cut);
  tg.Scan("run:trig:dt:sqrt(dx*dx+dy*dy):dx:dy:dz", cut);
}
