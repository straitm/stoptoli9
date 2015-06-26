#include "consts.h"
#include "math.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TMinuit.h"
#include "TTree.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TROOT.h"
#include <stdio.h>

const double muClife = 2028.; // muon lifetime in C-12, ns
const double hn_e = 2.224573;

const double mulife = 2196.9811;
const double muClife_err = 2.;
const double capprob12 = 1-muClife/mulife;
const double errcapprob12 = (1-(muClife+muClife_err)/mulife)/2
                           -(1-(muClife-muClife_err)/mulife)/2;

#ifdef HP
  const double No16cap_betan = n_o16cap_betan_hp*livetime;
  const double Nc12captarget = n_c12captarget_hp*livetime;
  const double Nc12cap       = n_c12cap_hp*livetime;
  const double Nc13cap       = n_c13cap_hp*livetime;
#else
  const double No16cap_betan = n_o16cap_betan*livetime;
  const double Nc12captarget = n_c12captarget*livetime;
  const double Nc12cap       = n_c12cap*livetime;
  const double Nc13cap       = n_c13cap*livetime;
#endif

/*
 * Prints the message once with the requested floating point precision
 * and in RED, then again with all digits in the default color, starting
 * with the first floating point number.
 *
 * If a precision is already given in the msg, the number is printed
 * with that precision both times.
 */
void printtwice(const char * const msg, const int digits, ...)
{
  char * bmsg = (char *)malloc(strlen(msg)+100); // Ha!
  char * pmsg = (char *)malloc(strlen(msg)+100); // Ha!

  // Just for fun...
  char * pmp = pmsg;
  char * bmp = bmsg;
  bool gotone = false;
  for(unsigned int i = 0; i <= strlen(msg); i++){
    switch(msg[i]){
      case '\0':
        *pmp++ = '\0';
        *bmp++ = '\0';
        break;
      case '%':
        switch(msg[i+1]){
          case 'e': case 'E': case 'f': case 'F':
          case 'g': case 'G': case 'a': case 'A':
            gotone = true;
            *pmp++ = '%';
            *bmp++ = '%';
            *pmp++ = '.';
            *pmp++ = digits+'0';
            break;
          default:
            *pmp++ = msg[i];
            if(gotone) *bmp++ = msg[i];
        }
        break;
      default:
        *pmp++ = msg[i];
        if(gotone) *bmp++ = msg[i];
        break;
    }
  }

  va_list ap;
  va_start(ap, digits);
  printf(RED);
  vprintf(pmsg, ap);
  printf(CLR);

  va_start(ap, digits);
  vprintf(bmsg, ap);
}

double logis(const double x,
             const double p0, const double p1, const double p2)
{
  return p0/(1+ exp(-(x - p1)/p2));
}

double corrmiche(const double e, const double me, const double mt)
{
  const double early = e /(logis(me, 0.402/hn_e, 203., 25.0) + 1);
  const double late  = e /(logis(me, 0.559/hn_e, 286., 73.0) + 1);
  return early + (late - early) * (mt - 3625.)/1250.;
}

void print_results(const double eff, const double energy,
                   const double n, const double nelo, const double neup)
{
  printtwice("\n%.0fkeV fitted number of events, raw: "
             "%f %f +%f\n", 1, energy, n, nelo, neup);

  const double n_ec = n/eff, nelo_ec = nelo/eff, neup_ec = neup/eff;

  printtwice("\n%.0fkeV fitted number of events, eff corrected: "
             "%f %f +%f\n", 1, energy, n_ec, nelo_ec, neup_ec);

  printtwice("\n%.0fkeV per C-12 nuclear capture: "
             "(%f %f +%f)%%\n", 2, energy,
             100*n_ec   /Nc12cap,
             100*nelo_ec/Nc12cap,
             100*neup_ec/Nc12cap);

  const double ratemult = capprob12/(Nc12cap*muClife)*1e6;

  printtwice("%.0fkeV rate: %f %f +%f e-3\n", 3, energy,
             n_ec*ratemult,
             nelo_ec*ratemult,
             neup_ec*ratemult);
}

void b12gammafinalfit(const int region = 1, const double lowt1 = 3008.)
{
  const double lowt   = region == 0? 4500 :region == 1?  lowt1: 2750;
  const double hight = 5500.; // ns
  const double highfq = region == 0? 215:region == 1? 215: 1000; // MeV


#ifdef HP
  // From muon counting.  Exact.
  const double fq_eff = region == 0?481414./1628874
                       :region == 1?481414./1628874
                       :1;
#else
  // From B-12-like event counting, with a bit over 1% stat error. Since
  // the most precise output of this fit has 10% stat errors, I won't
  // obsess over that.
  const double fq_eff = region == 0?0.3857
                       :region == 1?0.3857
                       :1;
#endif

  const double li8life = 839.9;

  TFile * f = new TFile(rootfile3up, "read");
  TTree * t = (TTree *)f->Get("t");

  const double b12hl = 20.2; // b12 half life, ms

  const double b12lowt = 2, b12hight = 60;
  const double acclowt = 10e3, acchight = 100e3;

  const double eff_eor_100s  = 0.9709;
  const double eff_eor = 1-(1-eff_eor_100s)*b12hight/100e3;

  const double fq_per_mev = 8300;

  const double b12ecutlow = 4;

  const double distcut = 400;

  // Relative to B-12 selection
  const double eff = 1
    * 0.981  // Subsequent muon veto efficiency 
    * eff_eor // timeleft cut
    * 0.8504 // B-12 energy cut
    * wholedet_dist400eff
    * (exp(-lowt/muClife)-exp(-hight/muClife)) // C-12 capture time
    * (exp(-b12lowt *log(2)/b12hl)
      -exp(-b12hight*log(2)/b12hl)) // B-12 beta decay time
    * fq_eff
  ;

  printtwice("Efficiency: %f%%\n", 1, 100*eff);

  TCanvas * c2 = new TCanvas;

  TH1D * ehist  = new TH1D("ehist" , "", 600, 0.7, 60.7);
  TH1D * bg     = new TH1D("bg"    , "", 600, 0.7, 60.7);
  TH1D * corrbg = new TH1D("corrbg", "", 600, 0.7, 60.7);

  t->Draw(Form("corrmiche(miche, fq/%f, micht) >> ehist",
               fq_per_mev),
         Form("!earlymich && "
         "latennear==0 && "
         "ndecay == 0 && "
         "e > %f && e < 15 && "
         "dt >%f && dt < %f && "
         "timeleft>%f && "
         "dist < %f && "
         "fq < %f && "
         "micht >= %f && micht < %f"
         , b12ecutlow, b12lowt, b12hight, b12hight, distcut,
         fq_per_mev*highfq, lowt, hight)
      , "e");

  t->Draw("corrmiche(miche, fq/8300, micht) >> bg",
         Form("!earlymich && "
         "latennear==0 && "
         //"ndecay == 0 && " how to handle this? Really need the first
         //event in lots of windows, which isn't convenient -- I think
         //this is close enough
         "e> %f && e < 15 && "
         "dt > %f && dt < %f && "
         "timeleft>%f && "
         "dist < %f && "
         "fq < %f && "
         "micht >= %f && micht < %f "
         , b12ecutlow, acclowt, acchight, acchight, distcut,
         fq_per_mev*highfq, lowt, hight)
      , "e");


  // Possible background from li-8 gammas, particularly at 980.8keV
  const double li8lowt = 300, li8hight = 5*li8life;

  t->Draw(Form("corrmiche(miche, fq/%f, micht) >> corrbg",
               fq_per_mev),
         Form("!earlymich && "
         "latennear==0 && "
         "e> %f && e < 15 && "
         "dt > %f && dt < %f && "
         "timeleft>100e3 && "
         "dist < %f && "
         "fq < %f && "
         "micht >= %f && micht < %f "
         , b12ecutlow, li8lowt, li8hight, distcut,
         fq_per_mev*highfq, lowt, hight)
      , "e");

  ehist->Draw("e");
  bg->SetLineColor(kRed);
  bg->SetMarkerColor(kRed);
  corrbg->SetLineColor(kViolet);
  corrbg->SetMarkerColor(kViolet);
  bg->Scale((b12hight - b12lowt)/(acchight-acclowt)/eff_eor_100s);
  bg->Draw("histsame");

  TF1 *bggg = new TF1("bggg",
    //                Simple Michel spectrum, stretched a little to
    //                account for the ~50% of positrons. It doesn't
    //                really matter if this is precise.
    "[0]+gaus(1)+gaus(4)+[7]*(3*(x/53.6)^2-2*(x/53.6)^3)",0,40);
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
  corrbg->Scale(
    (exp(-b12lowt/li8life) - exp(-b12hight/li8life))/
    (exp(-li8lowt/li8life) - exp(-li8hight/li8life))/eff_eor_100s);

  corrbg->Draw("histsame");

  TF1 * corrbgfit = new TF1("corrbgfit", "gaus(0)", 0.7, 2);
  corrbgfit->SetLineColor(kViolet);
  const double li8gamma = 0.9808;
  corrbgfit->FixParameter(1, li8gamma);

  // energy resolution
  corrbgfit->FixParameter(2, sqrt(pow(0.077, 2)*li8gamma +
                                  pow(0.018*li8gamma, 2) +
                                  pow(0.167, 2)));
  corrbg->Fit("corrbgfit", "li");
  corrbgfit->SetNpx(300);

  TF1 *gg = new TF1("gg",
  "[0] +"
  " [3]/sqrt([23]^2*[4]+([24]*[4])**2+[25]**2)*"
   "exp(-(((x- [4]*[1])/sqrt([23]^2*[4]+([24]*[4])**2+[25]**2))**2)/2)+"
  " [5]/sqrt([23]^2*[6]+([24]*[6])**2+[25]**2)*"
   "exp(-(((x- [6]*[1])/sqrt([23]^2*[6]+([24]*[6])**2+[25]**2))**2)/2)+"
  " [7]/sqrt([23]^2*[8]+([24]*[8])**2+[25]**2)*"
   "exp(-(((x- [8]*[1])/sqrt([23]^2*[8]+([24]*[8])**2+[25]**2))**2)/2)+"
  " [9]/sqrt([23]^2*[10]+([24]*[10])**2+[25]**2)*"
   "exp(-(((x-[10]*[1])/sqrt([23]^2*[10]+([24]*[10])**2+[25]**2))**2)/2)+"
  "[11]/sqrt([23]^2*[12]+([24]*[12])**2+[25]**2)*"
   "exp(-(((x-[12]*[1])/sqrt([23]^2*[12]+([24]*[12])**2+[25]**2))**2)/2)+"
  "[13]/sqrt([23]^2*[14]+([24]*[14])**2+[25]**2)*"
   "exp(-(((x-[14]*[1])/sqrt([23]^2*[14]+([24]*[14])**2+[25]**2))**2)/2)+"
 // bg
   "[15]+gaus(16)+gaus(19)+[22]*(3*(x/52.8)^2-2*(x/52.8)^3) +"
   "gaus(26) + gaus(29)"
    , 0,15);

  const char * ggpars[32] = { "accidentals", "energyscale", "unused",
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
    const double hn_t = 179., gn_t = 28.;
    const double nomEerr = sqrt(a*a*hn_e + b*b*hn_e*hn_e + c*c);
    const double normNN =
        23.3 * (exp(-lowt/1000/hn_t)-exp(-5.5/hn_t))
       + 1.5 * (exp(-lowt/1000/gn_t)-exp(-5.5/gn_t));
    gg->FixParameter(29, normNN/nomEerr/sqrt(2*M_PI));
    gg->FixParameter(30, hn_e);
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

  const double a = gg->GetParameter("er_a");
  const double b = gg->GetParameter("er_b");
  const double c = gg->GetParameter("er_c");

  if(region == 0) {
    {
      TF1 * mypeak0a = new TF1("mypeak0a",
        "[2]/[1]* exp(-(((x- [3]*[0])/[1])**2)/2)", 0, 15);
      mypeak0a->SetNpx(400);

      mypeak0a->SetParameter(0, gg->GetParameter(1));
      const double E = gg->GetParameter("e1");
      mypeak0a->SetParameter(1, sqrt(a*a*E + b*b*E*E + c*c));

      for(int i = 2; i < 4; i++)
        mypeak0a->SetParameter(i, gg->GetParameter(i+1));

      mypeak0a->Draw("same");

      // The number of events in the peak, corrected for efficiency.
      const double nev = mypeak0a->Integral(0, 15)/ehist->GetBinWidth(1);
      const double neveup = gMinuit->fErp[0]/gg->GetParameter("n1")*nev;
      const double nevelo = gMinuit->fErn[0]/gg->GetParameter("n1")*nev;

      print_results(eff, 953, nev, nevelo, neveup);
    }
    {
      TF1 * mypeak0b = new TF1("mypeak0b",
        "[2]/[1]* exp(-(((x- [3]*[0])/[1])**2)/2)", 0, 15);
      mypeak0b->SetNpx(400);

      mypeak0b->SetParameter(0, gg->GetParameter(1));
      const double E = gg->GetParameter("e2");
      mypeak0b->SetParameter(1, sqrt(a*a*E + b*b*E*E + c*c));

      for(int i = 2; i < 4; i++)
        mypeak0b->SetParameter(i,
          gg->GetParameter(i+gg->GetParNumber("n2")-2));

      mypeak0b->Draw("same");

      const double nev = mypeak0b->Integral(0, 15)/ehist->GetBinWidth(1);
      const double neveup = gMinuit->fErp[1]/gg->GetParameter("n2")*nev;
      const double nevelo = gMinuit->fErn[1]/gg->GetParameter("n2")*nev;

      print_results(eff, 1674, nev, nevelo, neveup);
    }
  }
  else if(region == 1){
    {
      TF1 * mypeak1a = new TF1("mypeak1a",
        "[2]/[1]* exp(-(((x- [3]*[0])/[1])**2)/2)", 0, 15);
      mypeak1a->SetNpx(400);

      mypeak1a->SetParameter(0, gg->GetParameter(1));
      const double E = gg->GetParameter("e3");
      mypeak1a->SetParameter(1, sqrt(a*a*E + b*b*E*E + c*c));

      for(int i = 2; i <= 3; i++)
        mypeak1a->SetParameter(i,
          gg->GetParameter(i+gg->GetParNumber("n3")-2));

      mypeak1a->Draw("same");

      const double nev = mypeak1a->Integral(0, 15)/ehist->GetBinWidth(1);
      const double neveup = gMinuit->fErp[2]/gg->GetParameter("n3")*nev;
      const double nevelo = gMinuit->fErn[2]/gg->GetParameter("n3")*nev;

      print_results(eff, 2621, nev, nevelo, neveup);
    }
    {
      TF1 * mypeak1b = new TF1("mypeak1b",
        "[2]/[1]* exp(-(((x- [3]*[0])/[1])**2)/2)", 0, 15);
      mypeak1b->SetNpx(400);

      mypeak1b->SetParameter(0, gg->GetParameter(1));
      const double E = gg->GetParameter("e5");
      mypeak1b->SetParameter(1, sqrt(a*a*E + b*b*E*E + c*c));

      for(int i = 2; i < 4; i++)
        mypeak1b->SetParameter(i,
          gg->GetParameter(i+gg->GetParNumber("n5")-2));

      mypeak1b->Draw("same");

      const double nev = mypeak1b->Integral(0, 15)/ehist->GetBinWidth(1);
      const double neveup = gMinuit->fErp[3]/gg->GetParameter("n5")*nev;
      const double nevelo = gMinuit->fErn[3]/gg->GetParameter("n5")*nev;

      print_results(eff, 3759, nev, nevelo, neveup);
    }
  }
  else if(region == 2){
    TF1 * mypeak2 = new TF1("mypeak2",
      "[2]/[1]* exp(-(((x- [3]*[0])/[1])**2)/2)", 0, 15);
    mypeak2->SetNpx(400);

    mypeak2->SetParameter(0, gg->GetParameter(1));
    const double E = gg->GetParameter("e6");
    mypeak2->SetParameter(1, sqrt(a*a*E + b*b*E*E + c*c));

    for(int i = 2; i < 4; i++)
      mypeak2->SetParameter(i,
        gg->GetParameter(i+gg->GetParNumber("n6")-2));

    mypeak2->Draw("same");

    const double nev = mypeak2->Integral(0, 15)/ehist->GetBinWidth(1);
    const double neveup = gMinuit->fErp[4]/gg->GetParameter("n6")*nev;
    const double nevelo = gMinuit->fErn[4]/gg->GetParameter("n6")*nev;

    print_results(eff, 9040, nev, nevelo, neveup);
  }

  ehist->GetXaxis()->SetRangeUser(0.7, 12);
}
