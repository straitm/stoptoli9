#include "consts.h"

void b12gammafinalfit(const int region = 1)
{
  const double lowt   = region == 0? 4500 :region == 1?  3008: 2750;
  const double highfq = region == 0? 1.8e6:region == 1? 1.8e6: 10e6;
 
  TFile * f = new TFile(rootfile0up, "read");
  TTree * t = (TTree *)f->Get("t");

  // C-12 -> B-12 included in raw count
  // C-13 -> B-13 also included
  // C-13 -> B-12+n mostly excluded by neutron cut
  // Conservatively take 0.5% syst normalization error due to these
  const double n_c12b12_withcuts = 11599.482876 * 162.837101/164.473163;

  // Relative to B-12 selection
  const double eff = 1
    // shared with B-12 * 0.962 // subsequent muons with 1ms veto
    // shared with B-12 * 0.977 // previous muons
    * 0.9999418/0.9709 // 200ms from end of run, better than B-12 selection
    // shared with B-12 * 0.82 // B-12 beta decay energy cut, eval'd for short tracks
    // shared with B-12 * 0.9177 // dist cut for B-12 beta decay
    * (exp(-lowt/1000*log(2)/2.028)-exp(-5.008*log(2)/2.028)) // B-12 capture time
    * (exp(-2.0*log(2)/20.2)-exp(-200.0*log(2)/20.2)) // B-12 beta decay time
  ;

  printf("%sEfficiency: %.1f%%%s\n", RED, 100*eff, CLR);
/*
  t->Draw("micht >> thist(50, 0, 6400)",
         Form("miche > 2.1 && miche < 3.1 && "
         "!earlymich && "
         "latennear==0 && "
         "ndecay == 0 && "
         "e> 4 && e < 15 && "
         "dt >2 && dt < 200 && "
         "timeleft>200 && "
         "dist < 400 && "
         "fq < %f  ", highfq)
      , "e"); */

  TCanvas * c2 = new TCanvas;
  c2->SetLogy();

  t->Draw("-0.026 + miche*1.011 - 0.0006*miche*miche >> ehist(600, 0.7, 60.7)",
         Form("!earlymich && "
         "latennear==0 && "
         "ndecay == 0 && "
         "e> 4 && e < 15 && "
         "dt >2 && dt < 200 && "
         "timeleft>200 && "
         "dist < 400 && "
         "fq < %f && "
         "micht >= %f && micht < 5500"
         , highfq, lowt)
      , "e");

  t->Draw("-0.026 + miche*1.011 - 0.0006*miche*miche >> bg(600, 0.7, 60.7)",
         Form("!earlymich && "
         "latennear==0 && "
         //"ndecay == 0 && " how to handle this?  Really need the first event in 
         // lots of windows, which isn't convenient
         "e> 4 && e < 15 && "
         "dt > 10000 && dt < 10000+19800*4 && "
         "timeleft>100e3 && "
         "dist < 400 && "
         "fq < %f && "
         "micht >= %f && micht < 5500 "
         , highfq, lowt)
      , "e");

  TH1 * ehist = gROOT->FindObject("ehist");
  ehist->Draw("e");
  TH1 * bg = gROOT->FindObject("bg");
  bg->SetLineColor(kRed);
  bg->SetMarkerColor(kRed);
  bg->Scale(0.0025/0.9709);
  bg->Draw("histsame");
  ehist->GetYaxis()->SetRangeUser(1e-3, 1e2);

  TF1 *bggg = new TF1("bggg",
    "[0]+gaus(1)+gaus(4)+[7]*(3*(x/52.8)^2-2*(x/52.8)^3)",0,40);
  bggg->SetNpx(400);
  bggg->SetLineColor(kRed);

  bggg->SetParameter(0,0.006);

  bggg->SetParameter(1,0.15);
  bggg->SetParameter(2,2.0);
  bggg->SetParameter(3,1.0);

  bggg->SetParameter(4,0.15);
  bggg->SetParameter(5,2.5);
  bggg->SetParameter(6,1.0);

  bggg->SetParameter(7,0.02);


  bggg->SetParLimits(0,0,1);

  bggg->SetParLimits(1,0,10);
  bggg->SetParLimits(2,0.7,2.3);
  bggg->SetParLimits(3,0.5,5);

  bggg->SetParLimits(4,0,10);
  bggg->SetParLimits(5,2.3,4);
  bggg->SetParLimits(6,0.5,5);

  bggg->SetParLimits(7,0,1);

  bg->Fit("bggg", "li", "", 0.7, 40);


  TF1 *gg = new TF1("gg",
      "[0] +"
      " [3]/[2]* exp(-(((x- [4]*[1])/[2])**2)/2) +"
      " [5]/[2]* exp(-(((x- [6]*[1])/[2])**2)/2) +"
      " [7]/[2]* exp(-(((x- [8]*[1])/[2])**2)/2) +"
      " [9]/[2]* exp(-(((x-[10]*[1])/[2])**2)/2) +"
      "[11]/[2]* exp(-(((x-[12]*[1])/[2])**2)/2) +"
      "[13]/[2]* exp(-(((x-[14]*[1])/[2])**2)/2) +"
 // bg
    "[15]+gaus(16)+gaus(19)+[22]*(3*(x/52.8)^2-2*(x/52.8)^3)"
      , 0,15);

  char * ggpars[23] = { "accidentals", "energyscale", "resolution",
  "n1", "e1",
  "n2", "e2",
  "n3", "e3",
  "n4", "e4",
  "n5", "e5",
  "n6", "e6",
  "bgconst",
  "bggaus1norm",
  "bggaus1mean",
  "bggaus1sig",
  "bggaus2norm",
  "bggaus2mean",
  "bggaus2sig",
  "bgmich" };

  for(int i = 0; i < gg->GetNpar(); i++) gg->SetParName(i, ggpars[i]);

  gg->SetNpx(400);

  for(int i = 3; i <= 13; i+=2) gg->SetParLimits(i, 0, 100);
  for(int i = 15; i <= 22; i++) 
    gg->FixParameter(i, bggg->GetParameter(i-15));

  gg->FixParameter(0, 0);
  gg->FixParameter(1, 1);
  gg->FixParameter(2, 0.22);

  gg->FixParameter(4, 0.95314);
  gg->FixParameter(6, 1.67365);
  gg->FixParameter(8, 2.6208);
  gg->FixParameter(10,2.723);
  gg->FixParameter(12,3.759);
  gg->FixParameter(14,9.040);


  gg->FixParameter(9, 0);

  ehist->Fit("gg", "lq", "", 0, 15);
  gg->ReleaseParameter(2);
  gg->SetParLimits(2, 0.16, 0.27);
  ehist->Fit("gg", "limq", "", 0, 15);

  bggg->Draw("same");


  gMinuit->Command("Set strategy 2");
  gMinuit->Command("MINOS 10000 8");
  gMinuit->Command("MINOS 10000 12");
  gMinuit->Command("MINOS 10000 14");

  gMinuit->Command("showmin");

  if(region == 1) {
    TF1 * mypeak = new TF1("mypeak",
      "[2]/[1]* exp(-(((x- [3]*[0])/[1])**2)/2)", 0, 15);
    mypeak->SetNpx(400);

    for(int i = 0; i < 2; i++) mypeak->SetParameter(i, gg->GetParameter(i+1));
    for(int i = 2; i < 4; i++) mypeak->SetParameter(i, gg->GetParameter(i+5));

    mypeak->Draw("same");

    const double ncorr = mypeak->Integral(0, 15)/eff/ehist->GetBinWidth(1);
    const double ncorreup = gMinuit->fErp[4]/gg->GetParameter(7) * ncorr;
    const double ncorrelo = gMinuit->fErn[4]/gg->GetParameter(7) * ncorr;

    printf("%s2.6208MeV events, eff corrected: %.1f%.1f+%.1f%s\n",
           RED, ncorr, ncorrelo, ncorreup, CLR);

    printf("%s2.6208MeV per B-12, eff corrected: (%.1f%.1f+%.1f)%%%s\n",
           RED, 100*ncorr/n_c12b12_withcuts,
                100*ncorrelo/n_c12b12_withcuts,
                100*ncorreup/n_c12b12_withcuts, CLR);

    printf("%s2.6208MeV rate, eff corrected: (%.2f%.2f+%.2f)e-3%s\n",
           RED, ncorr/n_c12b12_withcuts*37.9,
                ncorrelo/n_c12b12_withcuts*37.9,
                ncorreup/n_c12b12_withcuts*37.9, CLR);
    printf("Looks from timing that it is <1/4 neutrons at 90%%CL\n");
    printf("More convincingly, there's no big peak of Gd-n\n");
  }

  if(region == 1){
    TF1 * mypeak2 = new TF1("mypeak2",
      "[2]/[1]* exp(-(((x- [3]*[0])/[1])**2)/2)", 0, 15);
    mypeak2->SetNpx(400);

    for(int i = 0; i < 2; i++) mypeak2->SetParameter(i, gg->GetParameter(i+1));
    for(int i = 2; i < 4; i++) mypeak2->SetParameter(i, gg->GetParameter(i+9));

    mypeak2->Draw("same");

    const double ncorr = mypeak2->Integral(0, 15)/eff/ehist->GetBinWidth(1);
    const double ncorreup = gMinuit->fErp[5]/gg->GetParameter(11) * ncorr;
    const double ncorrelo = gMinuit->fErn[5]/gg->GetParameter(11) * ncorr;

    printf("%s3.7590MeV events, eff corrected: %.1f%.1f+%.1f%s\n",
           RED, ncorr, ncorrelo, ncorreup, CLR);

    printf("%s3.7590MeV per B-12, eff corrected: (%.2f%.2f+%.2f)%%%s\n",
           RED, 100*ncorr/n_c12b12_withcuts,
                100*ncorrelo/n_c12b12_withcuts,
                100*ncorreup/n_c12b12_withcuts, CLR);

    printf("%s3.7590MeV rate, eff corrected: %.3f%.3f+%.3f e-3%s\n",
           RED, ncorr/n_c12b12_withcuts*37.9,
                ncorrelo/n_c12b12_withcuts*37.9,
                ncorreup/n_c12b12_withcuts*37.9, CLR);
  }

  if(region == 2){
    TF1 * mypeak3 = new TF1("mypeak3",
      "[2]/[1]* exp(-(((x- [3]*[0])/[1])**2)/2)", 0, 15);
    mypeak3->SetNpx(400);

    for(int i = 0; i < 2; i++) mypeak3->SetParameter(i, gg->GetParameter(i+1));
    for(int i = 2; i < 4; i++) mypeak3->SetParameter(i, gg->GetParameter(i+11));

    mypeak3->Draw("same");

    const double ncorr = mypeak3->Integral(0, 15)/eff/ehist->GetBinWidth(1);
    const double ncorreup = gMinuit->fErp[6]/gg->GetParameter(13) * ncorr;
    const double ncorrelo = gMinuit->fErn[6]/gg->GetParameter(13) * ncorr;

    printf("%s9.0400MeV events, eff corrected: %.1f%.1f+%.1f%s\n",
           RED, ncorr, ncorrelo, ncorreup, CLR);

    printf("%s9.0400MeV per B-12, eff corrected: (%.2f%.2f+%.2f)%%%s\n",
           RED, 100*ncorr/n_c12b12_withcuts,
                100*ncorrelo/n_c12b12_withcuts,
                100*ncorreup/n_c12b12_withcuts, CLR);

    printf("%s9.0400MeV rate, eff corrected: %.3f%.3f+%.3f e-3%s\n",
           RED, ncorr/n_c12b12_withcuts*37.9,
                ncorrelo/n_c12b12_withcuts*37.9,
                ncorreup/n_c12b12_withcuts*37.9, CLR);
  }

  ehist->GetXaxis()->SetRangeUser(0, 12);

  printf("Efficiency not redone!!!!!!!!!!!!!!!!\n");
}
