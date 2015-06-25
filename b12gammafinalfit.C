#include "consts.h"

double logis(const double x, const double p0, const double p1, const double p2)
{
  return p0/(1+ exp(-(x - p1)/p2));
}

double corrmiche(const double e, const double me, const double mt)
{
  const double nh = 2.224573;
  const double early = e /(logis(me, 0.402/nh, 203., 25.0) + 1);
  const double late  = e /(logis(me, 0.559/nh, 286., 73.0) + 1);
  return early + (late - early) * (mt - 3625.)/1250.;
}

void b12gammafinalfit(const int region = 1, const double lowt1 = 3008.)
{
  const double lowt   = region == 0? 4500 :region == 1?  lowt1: 2750;
  const double highfq = region == 0? 215:region == 1? 215: 1000;
 
  TFile * f = new TFile(rootfile3up, "read");
  TTree * t = (TTree *)f->Get("t");

  // C-12 -> B-12 included in raw count
  // C-13 -> B-13 also included
  // C-13 -> B-12+n mostly excluded by neutron cut
  // Conservatively take 0.5% syst normalization error due to these
  const double n_c12b12_withcuts = (region < 2? 8895.779311:22871.08493)
    * 162.837101/164.473163;

  // Relative to B-12 selection
  const double eff = 1
    // shared with B-12 * 0.962 // subsequent muons with 1ms veto
    * 0.9999418/0.9709 // 200ms from end of run, better than B-12 selection
    // shared with B-12 * 0.82 // B-12 beta decay energy cut, eval'd for short tracks
    // shared with B-12 * 0.9177 // dist cut for B-12 beta decay
    * (exp(-lowt/1000*log(2)/2.028)-exp(-5.5*log(2)/2.028)) // B-12 capture time
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
         "fq < %f  ", 8300*highfq)
      , "e"); */

  TCanvas * c2 = new TCanvas;
//  c2->SetLogy();

  const double b12lowt = 2, b12hight = 60;

  const double acclowt = 10e3, acchight = 100e3;

  t->Draw("corrmiche(miche, fq/8300, micht) >> ehist(600, 0.7, 60.7)",
         Form("!earlymich && "
         "latennear==0 && "
         "ndecay == 0 && "
         "e> 4 && e < 15 && "
         "dt >%f && dt < %f && "
         "timeleft>%f && "
         "dist < 400 && "
         "fq < %f && "
         "micht >= %f && micht < 5500"
         , b12lowt, b12hight, b12hight, 8300*highfq, lowt)
      , "e");

  t->Draw("corrmiche(miche, fq/8300, micht) >> bg(600, 0.7, 60.7)",
         Form("!earlymich && "
         "latennear==0 && "
         //"ndecay == 0 && " how to handle this?  Really need the first event in 
         // lots of windows, which isn't convenient -- I think this is close enough
         "e> 4 && e < 15 && "
         "dt > %f && dt < %f && "
         "timeleft>%f && "
         "dist < 400 && "
         "fq < %f && "
         "micht >= %f && micht < 5500 "
         , acclowt, acchight, acchight, 8300*highfq, lowt)
      , "e");


  // Possible background from li-8 gammas, particularly at 980.8keV
  const double li8lowt = 300, li8hight = 5*839.9;

  t->Draw("corrmiche(miche, fq/8300, micht) >> corrbg(600, 0.7, 60.7)",
         Form("!earlymich && "
         "latennear==0 && "
         "e> 4 && e < 15 && "
         "dt > %f && dt < %f && "
         "timeleft>100e3 && "
         "dist < 400 && "
         "fq < %f && "
         "micht >= %f && micht < 5500 "
         , li8lowt, li8hight, 8300*highfq, lowt)
      , "e");

  TH1 * ehist = gROOT->FindObject("ehist");
  ehist->Draw("e");
  TH1 * bg = gROOT->FindObject("bg");
  TH1 * corrbg = gROOT->FindObject("corrbg");
  bg->SetLineColor(kRed);
  bg->SetMarkerColor(kRed);
  corrbg->SetLineColor(kViolet);
  corrbg->SetMarkerColor(kViolet);
  bg->Scale((b12hight - b12lowt)/(acchight-acclowt)/0.9709);
  bg->Draw("histsame");
  //ehist->GetYaxis()->SetRangeUser(1e-3, 1e2);

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

  // corrbg has accidentals in it too
  corrbg->Add(bggg, -(li8hight-li8lowt)/(b12hight-b12lowt));

  // And scale the accidental-substracted version to the 
  // expected amount in the signal window
  corrbg->Scale((exp(-b12lowt/839.9) - exp(-b12hight/839.9))/
                (exp(-li8lowt/839.9) - exp(-li8hight/389.9))/0.9709);

  corrbg->Draw("histsame");

  TF1 * corrbgfit = new TF1("corrbgfit", "gaus(0)", 0.7, 2);
  corrbgfit->SetLineColor(kViolet);
  const double li8gamma = 0.9808;
  corrbgfit->FixParameter(1, li8gamma);
  corrbgfit->FixParameter(2, sqrt(0.077*0.077*li8gamma + pow(0.018*li8gamma, 2) + 0.167*0.167));
  corrbg->Fit("corrbgfit", "li");
  corrbgfit->SetNpx(300);

  TF1 *gg = new TF1("gg",
   "[0] +"
   " [3]/sqrt([23]^2*[4]+([24]*[4])**2+[25]**2)* exp(-(((x- [4]*[1])/sqrt([23]^2*[4]+([24]*[4])**2+[25]**2))**2)/2) +"
   " [5]/sqrt([23]^2*[6]+([24]*[6])**2+[25]**2)* exp(-(((x- [6]*[1])/sqrt([23]^2*[6]+([24]*[6])**2+[25]**2))**2)/2) +"
   " [7]/sqrt([23]^2*[8]+([24]*[8])**2+[25]**2)* exp(-(((x- [8]*[1])/sqrt([23]^2*[8]+([24]*[8])**2+[25]**2))**2)/2) +"
   " [9]/sqrt([23]^2*[10]+([24]*[10])**2+[25]**2)* exp(-(((x-[10]*[1])/sqrt([23]^2*[10]+([24]*[10])**2+[25]**2))**2)/2) +"
   "[11]/sqrt([23]^2*[12]+([24]*[12])**2+[25]**2)* exp(-(((x-[12]*[1])/sqrt([23]^2*[12]+([24]*[12])**2+[25]**2))**2)/2) +"
   "[13]/sqrt([23]^2*[14]+([24]*[14])**2+[25]**2)* exp(-(((x-[14]*[1])/sqrt([23]^2*[14]+([24]*[14])**2+[25]**2))**2)/2) +"
 // bg
   "[15]+gaus(16)+gaus(19)+[22]*(3*(x/52.8)^2-2*(x/52.8)^3) + gaus(26) + gaus(29)"
      , 0,15);

  char * ggpars[32] = { "accidentals", "energyscale", "unused",
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
  "bgmich", "er_a", "er_b", "er_c", "li8n", "li8m", "li8s",
  "hn_n", "hn_m", "hn_s"
  };

  for(int i = 0; i < gg->GetNpar(); i++) gg->SetParName(i, ggpars[i]);

  gg->SetNpx(400);

  for(int i = 3; i <= 13; i+=2) gg->SetParLimits(i, 0, 100);
  for(int i = 15; i <= 22; i++) 
    gg->FixParameter(i, bggg->GetParameter(i-15));

  gg->FixParameter(23, 0.077);
  gg->FixParameter(24, 0.018);
  gg->FixParameter(25, 0.017);

  gg->FixParameter(0, 0);
  gg->FixParameter(1, 1);
  gg->FixParameter(2, 0.22);

  gg->FixParameter(4, 0.95314);
  gg->FixParameter(6, 1.67365);
  gg->FixParameter(8, 2.6208);

  gg->FixParameter(9, 0);
  gg->FixParameter(10,0);

  gg->FixParameter(12,3.759);
  gg->FixParameter(14,9.040);

  gg->FixParameter(26, corrbgfit->GetParameter(0));
  gg->FixParameter(27, corrbgfit->GetParameter(1));
  gg->FixParameter(28, corrbgfit->GetParameter(2));

  // contamination by Hn events.
  {
    const double a = gg->GetParameter("er_a");
    const double b = gg->GetParameter("er_b");
    const double c = gg->GetParameter("er_c");
    const double nomEerr = sqrt(a*a*2.223 + b*b*2.223*2.223 + c*c);
    const double normNN = 23.3 * (exp(-lowt/1000*log(2)/179.)-exp(-5.5*log(2)/179))
                         + 1.5 * (exp(-lowt/1000*log(2)/28.)-exp(-5.5*log(2)/28.));
    gg->FixParameter(29, normNN/nomEerr/sqrt(2*3.14159265258));
    gg->FixParameter(30, 2.223);
    gg->FixParameter(31, nomEerr);
  }

  ehist->Fit("gg", "lq", "", 0.7, 15);
  gg->ReleaseParameter(25);
  gg->SetParLimits(25, 0, 1.92530e-01 + 1.81760e-02);
  ehist->Fit("gg", "limq", "", 0.7, 15);

  bggg->Draw("same");
  corrbgfit->Draw("same");


  gMinuit->Command("Set strategy 2");
  if(region == 0){
    gMinuit->Command("MINOS 10000 4");
    gMinuit->Command("MINOS 10000 6");
  }
  else{
    gMinuit->Command("MINOS 10000 8");
    gMinuit->Command("MINOS 10000 12");
    gMinuit->Command("MINOS 10000 14");
  }

  gMinuit->Command("showmin");

  if(region == 0) {
    {
      TF1 * mypeak0a = new TF1("mypeak0a",
        "[2]/[1]* exp(-(((x- [3]*[0])/[1])**2)/2)", 0, 15);
      mypeak0a->SetNpx(400);

      mypeak0a->SetParameter(0, gg->GetParameter(1));
      const double E = gg->GetParameter("e1");
      const double a = gg->GetParameter("er_a");
      const double b = gg->GetParameter("er_b");
      const double c = gg->GetParameter("er_c");
      mypeak0a->SetParameter(1, sqrt(a*a*E + b*b*E*E + c*c));
 
      for(int i = 2; i < 4; i++) mypeak0a->SetParameter(i, gg->GetParameter(i+1));

      mypeak0a->Draw("same");

      const double ncorr = mypeak0a->Integral(0, 15)/eff/ehist->GetBinWidth(1);
      const double ncorreup = gMinuit->fErp[0]/gg->GetParameter("n1") * ncorr;
      const double ncorrelo = gMinuit->fErn[0]/gg->GetParameter("n1") * ncorr;

      printf("\n%s953keV per B-12, eff corrected: (%.2f%.2f+%.2f)%%%s\n",
             RED, 100*ncorr/n_c12b12_withcuts,
                  100*ncorrelo/n_c12b12_withcuts,
                  100*ncorreup/n_c12b12_withcuts, CLR);

      printf("%s953keV rate, eff corrected: %.3f%.3f+%.3f e-3%s\n",
             RED, ncorr/n_c12b12_withcuts*7.05,
                  ncorrelo/n_c12b12_withcuts*7.05,
                  ncorreup/n_c12b12_withcuts*7.05, CLR);
    
    }
    {
      TF1 * mypeak0b = new TF1("mypeak0b",
        "[2]/[1]* exp(-(((x- [3]*[0])/[1])**2)/2)", 0, 15);
      mypeak0b->SetNpx(400);

      mypeak0b->SetParameter(0, gg->GetParameter(1));
      const double E = gg->GetParameter("e2");
      const double a = gg->GetParameter("er_a");
      const double b = gg->GetParameter("er_b");
      const double c = gg->GetParameter("er_c");
      mypeak0b->SetParameter(1, sqrt(a*a*E + b*b*E*E + c*c));

      for(int i = 2; i < 4; i++) mypeak0b->SetParameter(i, gg->GetParameter(i+gg->GetParNumber("n2")-2));

      mypeak0b->Draw("same");

      const double ncorr = mypeak0b->Integral(0, 15)/eff/ehist->GetBinWidth(1);
      const double ncorreup = gMinuit->fErp[1]/gg->GetParameter("n2") * ncorr;
      const double ncorrelo = gMinuit->fErn[1]/gg->GetParameter("n2") * ncorr;

      printf("\n%s1674keV per B-12, eff corrected: (%.2f%.2f+%.2f)%%%s\n",
             RED, 100*ncorr/n_c12b12_withcuts,
                  100*ncorrelo/n_c12b12_withcuts,
                  100*ncorreup/n_c12b12_withcuts, CLR);

      printf("%s1674keV rate, eff corrected: %.3f%.3f+%.3f e-3%s\n",
             RED, ncorr/n_c12b12_withcuts*7.05,
                  ncorrelo/n_c12b12_withcuts*7.05,
                  ncorreup/n_c12b12_withcuts*7.05, CLR);
    }
  }
  else if(region == 1){
    {
      TF1 * mypeak1a = new TF1("mypeak1a",
        "[2]/[1]* exp(-(((x- [3]*[0])/[1])**2)/2)", 0, 15);
      mypeak1a->SetNpx(400);

      mypeak1a->SetParameter(0, gg->GetParameter(1));
      const double E = gg->GetParameter("e3");
      const double a = gg->GetParameter("er_a");
      const double b = gg->GetParameter("er_b");
      const double c = gg->GetParameter("er_c");
      mypeak1a->SetParameter(1, sqrt(a*a*E + b*b*E*E + c*c));

      for(int i = 2; i <= 3; i++)
        mypeak1a->SetParameter(i, gg->GetParameter(i+gg->GetParNumber("n3")-2));

      mypeak1a->Draw("same");

      const double ncorr = mypeak1a->Integral(0, 15)/eff/ehist->GetBinWidth(1);
      const double ncorreup = gMinuit->fErp[2]/gg->GetParameter("n3") * ncorr;
      const double ncorrelo = gMinuit->fErn[2]/gg->GetParameter("n3") * ncorr;

      printf("\n%s2621keV per B-12, eff corrected: (%.2f%.2f+%.2f)%%%s\n",
             RED, 100*ncorr/n_c12b12_withcuts,
                  100*ncorrelo/n_c12b12_withcuts,
                  100*ncorreup/n_c12b12_withcuts, CLR);

      printf("%s2621keV rate, eff corrected: %.3f%.3f+%.3f e-3%s\n",
             RED, ncorr/n_c12b12_withcuts*7.05,
                  ncorrelo/n_c12b12_withcuts*7.05,
                  ncorreup/n_c12b12_withcuts*7.05, CLR);
    
    }
    {
      TF1 * mypeak1b = new TF1("mypeak1b",
        "[2]/[1]* exp(-(((x- [3]*[0])/[1])**2)/2)", 0, 15);
      mypeak1b->SetNpx(400);

      mypeak1b->SetParameter(0, gg->GetParameter(1));
      const double E = gg->GetParameter("e5");
      const double a = gg->GetParameter("er_a");
      const double b = gg->GetParameter("er_b");
      const double c = gg->GetParameter("er_c");
      mypeak1b->SetParameter(1, sqrt(a*a*E + b*b*E*E + c*c));

      for(int i = 2; i < 4; i++) mypeak1b->SetParameter(i, gg->GetParameter(i+gg->GetParNumber("n5")-2));

      mypeak1b->Draw("same");

      const double ncorr = mypeak1b->Integral(0, 15)/eff/ehist->GetBinWidth(1);
      const double ncorreup = gMinuit->fErp[3]/gg->GetParameter("n5") * ncorr;
      const double ncorrelo = gMinuit->fErn[3]/gg->GetParameter("n5") * ncorr;

      printf("\n%s3759keV per B-12, eff corrected: (%.2f%.2f+%.2f)%%%s\n",
             RED, 100*ncorr/n_c12b12_withcuts,
                  100*ncorrelo/n_c12b12_withcuts,
                  100*ncorreup/n_c12b12_withcuts, CLR);

      printf("%s3759keV rate, eff corrected: %.3f%.3f+%.3f e-3%s\n",
             RED, ncorr/n_c12b12_withcuts*7.05,
                  ncorrelo/n_c12b12_withcuts*7.05,
                  ncorreup/n_c12b12_withcuts*7.05, CLR);
    }
  }
  else if(region == 2){
    TF1 * mypeak2 = new TF1("mypeak2",
      "[2]/[1]* exp(-(((x- [3]*[0])/[1])**2)/2)", 0, 15);
    mypeak2->SetNpx(400);

    mypeak2->SetParameter(0, gg->GetParameter(1));
    const double E = gg->GetParameter("e6");
    const double a = gg->GetParameter("er_a");
    const double b = gg->GetParameter("er_b");
    const double c = gg->GetParameter("er_c");
    mypeak2->SetParameter(1, sqrt(a*a*E + b*b*E*E + c*c));

    for(int i = 2; i < 4; i++) mypeak2->SetParameter(i, gg->GetParameter(i+gg->GetParNumber("n6")-2));

    mypeak2->Draw("same");

    const double ncorr = mypeak2->Integral(0, 15)/eff/ehist->GetBinWidth(1);
    const double ncorreup = gMinuit->fErp[4]/gg->GetParameter("n6") * ncorr;
    const double ncorrelo = gMinuit->fErn[4]/gg->GetParameter("n6") * ncorr;

    printf("\n%s9040keV per B-12, eff corrected: (%.2f%.2f+%.2f)%%%s\n",
           RED, 100*ncorr/n_c12b12_withcuts,
                100*ncorrelo/n_c12b12_withcuts,
                100*ncorreup/n_c12b12_withcuts, CLR);

    printf("%s9040keV rate, eff corrected: %.3f%.3f+%.3f e-3%s\n",
           RED, ncorr/n_c12b12_withcuts*7.05,
                ncorrelo/n_c12b12_withcuts*7.05,
                ncorreup/n_c12b12_withcuts*7.05, CLR);
  }

  ehist->GetXaxis()->SetRangeUser(0.7, 12);
}
