#include "unistd.h"
#include "consts.h"
#include "carbondenominators_finalfit.out.h"
#include "li8_finalfit.out.h"
#include "fullb12_finalfit.out.h"


#include <string>
#include "math.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TMinuit.h"
#include "TMath.h"
#include "Math/QuantFuncMathCore.h" // for ROOT::Math
#include "Math/ProbFuncMathCore.h" // for ROOT::Math
#include "TTree.h"
#include "TString.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TROOT.h"
#include <stdio.h>

using std::string;

const int ggnnpars = 5;

TCanvas * c2 = new TCanvas("c2", "", 0, 0, 500, 400);

TMinuit * mn = NULL;

TH1D * bg     = new TH1D("bg"    , "", 600, 0.7, 60.7);
TH1D * bgn    = new TH1D("bgn"   , "", 600, 0.7, 60.7);
TH1D * corrbg = new TH1D("corrbg", "", 600, 0.7, 60.7);
TH1D * corrbgn= new TH1D("corrbgn","", 600, 0.7, 60.7);
TH1D * ehist  = new TH1D("ehist" , "", 600, 0.7, 60.7);
TH1D * ehistn = new TH1D("ehistn", "", 600, 0.7, 60.7);

TF1 *gg = NULL, *ggn = NULL;

double n_to_0_bg_rat = 0; // set later

const double neff = 0.858;
const double neff_e = 0.01;

// Gamma energies, plus 6keV for the B-12 recoil
const double G1E = 0.95314 + 0.006;
const double G2E = 1.67365 + 0.006;
const double G3E = 2.6208 + 0.006;
const double G4E = 3.759 + 0.006;

// Measured probablity of getting one accidental neutron.  These
// are *detected* neutrons, so don't apply efficiency to them.
const double paccn = 1.1e-4;

const double hn_e = 2.224573;

const double mulife = 2196.9811;
const double capprob12 = 1-lifetime_c12/mulife;
const double errcapprob12 = (1-(lifetime_c12+lifetime_c12_err)/mulife)/2
                           -(1-(lifetime_c12-lifetime_c12_err)/mulife)/2;

const double capprob13 = 1-lifetime_c13/mulife;
const double errcapprob13 = (1-(lifetime_c12+lifetime_c13_err)/mulife)/2
                           -(1-(lifetime_c12-lifetime_c13_err)/mulife)/2;

const double Nc12cap = n_c12cap*livetime;
const double Nc13cap = n_c13cap*livetime;

const double li8hl = 839.9;

// Possible background from li-8 gammas, particularly at 980.8keV
const double li8lowt = 300, li8hight = 5*li8hl;

const double b12hl = 20.20; // b12 half life, ms

const double b12lowt = 2, b12hight = 60;
const double acclowt = 10e3, acchight = 100e3;

const double eff_eor_100s  = 0.9709;
const double eff_eor_b12 = 1-(1-eff_eor_100s)*b12hight/100e3;
const double eff_eor_li8 = 1-(1-eff_eor_100s)*li8hight/100e3;
const double eff_eor_acc = 1-(1-eff_eor_100s)*acchight/100e3;

const double fq_per_mev = 8300;

const double b12ecutlow = 4;

const double distcut = 400;

// Will be set by the corrbgfit (a.k.a. Li-8 fit)
// and used as pull terms
double corrbgnorm_nom = -1, corrbgnorm_err = -1;

double max(const double a, const double b)
{
  return a > b ? a : b;
}

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
        gotone = true;
        switch(msg[i+1]){
          case 'e': case 'E': case 'f': case 'F':
          case 'g': case 'G': case 'a': case 'A':
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

#ifndef __CINT__
  va_list ap;
  va_start(ap, digits);
  printf(RED);
  vprintf(pmsg, ap);
  printf(CLR);

  va_start(ap, digits);
  vprintf(bmsg, ap);
#endif

  free(bmsg);
  free(pmsg);
}

double logis(const double x,
             const double p0, const double p1, const double p2)
{
  return p0/(1+ exp(-(x - p1)/p2));
}

int whichcorr = 0;

// Return a michel/gamma/neutron energy adjusted for the baseline
// shift (presumably) after a muon.
double corrmiche(const double e, const double me)
{
  double pars[7][3] ={{ 0.583, 287., 94.2}, // nominal
                      { 0.609, 307., 109.}, // shifted 1sigma each way
                      { 0.562, 273., 83.6},
                      { 0.567, 269., 83.0},
                      { 0.606, 310., 109.},
                      { 0.568, 273., 79.5},
                      { 0.604, 306., 112.}};
         
  return e  - logis(me, pars[whichcorr][0],
                        pars[whichcorr][1],
                        pars[whichcorr][2]);
}

// Returns a pretty darn good approximation to the plot in 
// doc5608-v4, figure 5, giving the time of Gd captures
TF1 * gdtime()
{
  static TF1 *f = new TF1("f",
   "[0]*exp(-x/[1])/(1+exp(-(x-[2])/[3]))/(1+exp(-(x-[4])/[5]))",0,600);
  f->SetParameter(0,1194.59);
  f->SetParameter(1,26.117);
  f->SetParameter(2,3.66788);
  f->SetParameter(3,2.99152);
  f->SetParameter(4,1.62591);
  f->SetParameter(5,0.678265);
  return f;
}

// Find the bayesian limit given a central value, lower error and upper
// error, on a process that is bounded from 0-1, with a flat prior. This
// assumes the likelihood has a gaussian shape on both sides, which of
// course is only an approximation, and really you should do it properly
// by evaluating the whole likelihood space.
double bayeslimit(const double central, const double elo,
                  const double ehi)
{
  const double lim = 0.9;
  const double a = ROOT::Math::gaussian_cdf(   central /fabs(elo), 1)-0.5;
  const double b = ROOT::Math::gaussian_cdf((1-central)/fabs(ehi), 1)-0.5;
  const double c = lim*(a+b) - a + 0.5;

  if(a/(a+b) < lim)
    return central + ROOT::Math::gaussian_quantile(c, 1) * fabs(ehi);
  else {
    fprintf(stderr, "unhandled case of limit below central value\n");
    return 0;
  }
}

// For Li-8 gammas
// eff: the efficiency
// energy: the gamma line energy
// n: the raw number of fitted events
// nelo: the lower error on that
// neup: the upper error on that
// addsystf: additional symmetric fractional systematic error
void print_results8(const double eff, const double energy,
                    const double n, const double nelo,
                    const double neup, const double addsystf)
{
  printtwice("\n%.0fkeV Li-8 fitted number of events, raw: "
             "%f %f +%f\n", 1, energy, n, nelo, neup);

  printf("Additional systematic: %f%%\n", 100*addsystf);

  const double n_ec = n/eff, nelo_ec = nelo/eff,
                             neup_ec = neup/eff;

  // sic: lo goes with up and vice versa
  const double li8syst_fup=
    sqrt(pow(probEightLiFromTwelveC_statlo/probEightLiFromTwelveC, 2) +pow(addsystf,2)),
               li8syst_flo=
   -sqrt(pow(probEightLiFromTwelveC_statup/probEightLiFromTwelveC, 2) +pow(addsystf,2));

  const double
    li8stat_lo = 100*nelo_ec/Nc12cap/probEightLiFromTwelveC,
    li8stat_up = 100*neup_ec/Nc12cap/probEightLiFromTwelveC,
    li8syst_lo = 100*n_ec/Nc12cap/probEightLiFromTwelveC * li8syst_flo,
    li8syst_up = 100*n_ec/Nc12cap/probEightLiFromTwelveC * li8syst_fup;

  const double percentval = 100*n_ec/Nc12cap/probEightLiFromTwelveC;

  printtwice("\nTECHNOTE results.tex probLiEightGammaFrac, "
             "probLiEightGammaFracCent: %.0fkeV per Li-8 production: "
             "(%f %f +%f (stat) %f +%f (syst) %f %f (tot))%%\n", 0, energy,
             percentval,
             li8stat_lo, li8stat_up,
             li8syst_lo, li8syst_up,
             -sqrt(pow(li8stat_lo,2) + pow(li8syst_lo,2)),
             sqrt(pow(li8stat_up,2) + pow(li8syst_up,2))
            );

  const double fraclimit =
         bayeslimit(percentval/100, li8stat_lo/100, li8stat_up/100)*100;

  // To be applied to the rate limit below
  const double bayescorrectionforrate = fraclimit/
    (percentval + li8stat_up*ROOT::Math::gaussian_quantile(0.9, 1));

  printtwice("\nTECHNOTE results.tex probLiEightGammaFracSimple: "
             "%.0fkeV 90%% upper limit per Li-8 production: "
             "%f%%\n", 0, energy,
             /* See comments for rate limit, below */
             fraclimit*(1 + li8syst_up/percentval*0.17)
            );

  const double ratemult = capprob12/(Nc12cap*lifetime_c12)*1e6;

  const double ratesyst_f = sqrt(
     pow(n_c12cap_forb12gamma_additional_ferr, 2)
   + pow(mum_count_e/mum_count, 2)
   + pow(lifetime_c12_err/lifetime_c12, 2)
   + pow(addsystf, 2));

  const double
    ratestat_lo = nelo_ec*ratemult,
    ratestat_up = neup_ec*ratemult,
    ratesyst = ratesyst_f*n_ec*ratemult;

  const double rateval = n_ec*ratemult;
  printtwice("\nTECHNOTE results.tex probLiEightGammaRate, "
             "probLiEightGammaRateCent: %.0fkeV rate: "
             "%f %f +%f (stat) %f +%f (syst) "
             "%f +%f (tot) e-3\n", 2, energy,
             rateval,
             ratestat_lo, ratestat_up,
             ratesyst, ratesyst,
             -sqrt(pow(ratestat_lo,2) + pow(ratesyst,2)),
             sqrt(pow(ratestat_up,2) + pow(ratesyst,2)));

  printtwice("\nTECHNOTE results.tex probLiEightGammaRateSimple: "
             "%.0fkeV rate 90%% upper limit: %fe-3\n", 2, energy,
  /* See comments in C-13 function below */
       (rateval + ratestat_up * ROOT::Math::gaussian_quantile(0.9, 1))
          * bayescorrectionforrate
             *(1 + ratesyst/rateval*0.17));

  puts("");
}

// For B12+n
// eff: the efficiency
// energy: the gamma line energy
// n: the raw number of fitted events
// nelo: the lower error on that
// neup: the upper error on that
// addsystf: additional symmetric fractional systematic error
void print_results13(const double eff, const double energy,
                   const double n, const double nelo, const double neup,
                   const double addsystf)
{
  printtwice("\n%.0fkeV B12 + n fitted number of events, raw: "
             "%f %f +%f\n", 1, energy, n, nelo, neup);

  printf("Additional systematic: %f%%\n", 100*addsystf);

  const double n_ec = n/eff, nelo_ec = nelo/eff,
                             neup_ec = neup/eff;

  // sic: up goes with lo and vice versa because these go on the bottom
  // of a ratio. i.e. if the total rate is higher, the fraction of a
  // gamma line is lower.
  const double b12syst_fup= sqrt(pow(probTwelveBFromThirteenC_statlo/
                                     probTwelveBFromThirteenC, 2)
                               + pow(addsystf,2)),
               b12syst_flo=-sqrt(pow(probTwelveBFromThirteenC_statup/
                                     probTwelveBFromThirteenC, 2)
                               + pow(addsystf,2));

  const double
    b12stat_lo = 100*nelo_ec/Nc13cap/probTwelveBFromThirteenC,
    b12stat_up = 100*neup_ec/Nc13cap/probTwelveBFromThirteenC,
    b12syst_lo = 100*n_ec/Nc13cap/probTwelveBFromThirteenC * b12syst_flo,
    b12syst_up = 100*n_ec/Nc13cap/probTwelveBFromThirteenC * b12syst_fup;

  const char * const linenumber = fabs(energy - 953) < 10? "One":
                                  fabs(energy - 1674) < 10? "Two":
                                  fabs(energy - 2621) < 10? "Three":
                                  fabs(energy - 3759) < 10? "Four":"??";

  const double percentval = 100*n_ec/Nc13cap/probTwelveBFromThirteenC;
  printtwice("\nTECHNOTE results.tex nprobGamma%sFrac, nprobGamma%sFracCent: "
             "%.0fkeV per B-12 + n production: "
             "(%f %f +%f (stat) %f +%f (syst) %f %f (tot))%%\n",
             0, linenumber, linenumber, energy,
             percentval,
             b12stat_lo, b12stat_up,
             b12syst_lo, b12syst_up,
             -sqrt(pow(b12stat_lo,2) + pow(b12syst_lo,2)),
             sqrt(pow(b12stat_up,2) + pow(b12syst_up,2))
            );

  const double fraclimit =
         bayeslimit(percentval/100, b12stat_lo/100, b12stat_up/100)*100;

  // To be applied to the rate limit below
  const double bayescorrectionforrate = fraclimit/
    (percentval + b12stat_up*ROOT::Math::gaussian_quantile(0.9, 1));

  printtwice("\nTECHNOTE results.tex nprobGamma%sFracSimple: "
             "%.0fkeV 90%% upper limit per B-12 + n production: "
             "%f%%\n", 0, linenumber, energy,
             /* See comments for systematic on the rate limit, below */
             fraclimit*(1 + b12syst_up/percentval*0.17)
            );

  const double ratemult = capprob13/(Nc13cap*lifetime_c13)*1e6;

  const double ratesyst_f = sqrt(
     pow(n_c12cap_forb12gamma_additional_ferr, 2)
   + pow(mum_count_e/mum_count, 2)
   + pow(lifetime_c12_err/lifetime_c12, 2)
   + pow(addsystf, 2));

  const double
    ratestat_lo = nelo_ec*ratemult,
    ratestat_up = neup_ec*ratemult,
    ratesyst = ratesyst_f*n_ec*ratemult;

  const double rateval = n_ec*ratemult;
  printtwice("\nTECHNOTE results.tex nprobGamma%sRate, nprobGamma%sRateCent: "
             "%.0fkeV rate: %f %f +%f (stat) %f +%f (syst) "
             "%f +%f (tot) e-3\n", 1, linenumber, linenumber, energy,
             rateval,
             ratestat_lo, ratestat_up,
             ratesyst, ratesyst,
             -sqrt(pow(ratestat_lo,2) + pow(ratesyst,2)),
             sqrt(pow(ratestat_up,2) + pow(ratesyst,2)));

  printtwice("\nTECHNOTE results.tex nprobGamma%sRateSimple: "
             "%.0fkeV rate 90%% upper limit: %fe-3\n", 0,
             linenumber, energy,
  /* Find upper limit as though there were no bounds, then correct that
     given the bounds of 0-100% production, as found above. */
       (rateval + ratestat_up * ROOT::Math::gaussian_quantile(0.9, 1))
          * bayescorrectionforrate
  /* Kludgy accounting for rate systematic. I know that with 33%
  denominator error, a limit expands by 5.5239%. I don't think this
  scales linearly, but it is also a *tiny* effect, so just fudge it. */
             *(1 + ratesyst/rateval*0.17));

  puts("");
}

// eff: the efficiency
// energy: the gamma line energy
// n: the raw number of fitted events
// nelo: the lower error on that
// neup: the upper error on that
// addsystf: additional symmetric fractional systematic error
void print_results(const double eff, const double energy,
                   const double n, const double nelo, const double neup,
                   const double addsystf)
{
  printtwice("\n%.0fkeV fitted number of events, raw: "
             "%f %f +%f\n", 1, energy, n, nelo, neup);

  printf("Additional systematic: %f%%\n", 100*addsystf);

  const double n_ec = n/eff, nelo_ec = nelo/eff,
                             neup_ec = neup/eff;

  /*printtwice("\n%.0fkeV fitted number of events, eff corrected: "
             "%f %f +%f\n", 1, energy, n_ec, nelo_ec, neup_ec); */

  /*printtwice("\n%.0fkeV per C-12 nuclear capture: "
             "(%f %f +%f)%%\n", 2, energy,
             100*n_ec   /Nc12cap,
             100*nelo_ec/Nc12cap,
             100*neup_ec/Nc12cap);*/

  // sic: up goes with lo and vice versa because these go on the bottom
  // of a ratio. i.e. if the total rate is higher, the fraction of a
  // gamma line is lower.
  const double b12syst_fup = 
    sqrt(pow(probTwelveBFromTwelveC_statlo/probTwelveBFromTwelveC, 2)
       + pow(addsystf,2)),
               b12syst_flo =
   -sqrt(pow(probTwelveBFromTwelveC_statup/probTwelveBFromTwelveC, 2)
       + pow(addsystf,2));

  const double
    b12stat_lo = 100*nelo_ec/Nc12cap/probTwelveBFromTwelveC,
    b12stat_up = 100*neup_ec/Nc12cap/probTwelveBFromTwelveC,
    b12syst_lo = 100*n_ec/Nc12cap/probTwelveBFromTwelveC * b12syst_flo,
    b12syst_up = 100*n_ec/Nc12cap/probTwelveBFromTwelveC * b12syst_fup;

  const char * const linenumber = fabs(energy - 953) < 10? "One":
                                  fabs(energy - 1674) < 10? "Two":
                                  fabs(energy - 2621) < 10? "Three":
                                  fabs(energy - 3759) < 10? "Four":"??";

  printtwice("\nTECHNOTE results.tex probGamma%sFrac, probGamma%sFracCent: "
             "%.0fkeV per bound B-12 production: "
             "(%f %f +%f (stat) %f +%f (syst) %f %f (tot))%%\n", 1, 
             linenumber, linenumber, energy,
             100*n_ec   /Nc12cap/probTwelveBFromTwelveC,
             b12stat_lo, b12stat_up,
             b12syst_lo, b12syst_up,
             -sqrt(pow(b12stat_lo,2) + pow(b12syst_lo,2)),
             sqrt(pow(b12stat_up,2) + pow(b12syst_up,2))
            );

  const double ratemult = capprob12/(Nc12cap*lifetime_c12)*1e6;

  const double ratesyst_f = sqrt(
     pow(n_c12cap_forb12gamma_additional_ferr, 2)
   + pow(mum_count_e/mum_count, 2)
   + pow(lifetime_c12_err/lifetime_c12, 2)
   + pow(addsystf, 2));

  const double
    ratestat_lo = nelo_ec*ratemult,
    ratestat_up = neup_ec*ratemult,
    ratesyst = ratesyst_f*n_ec*ratemult;

  printtwice("\nTECHNOTE results.tex probGamma%sRate, probGamma%sRateCent: "
             "%.0fkeV rate: %f %f +%f (stat) %f +%f (syst) "
             "%f +%f (tot) e-3\n", 2, linenumber, linenumber, energy,
             n_ec*ratemult,
             ratestat_lo, ratestat_up,
             ratesyst, ratesyst,
             -sqrt(pow(ratestat_lo,2) + pow(ratesyst,2)),
             sqrt(pow(ratestat_up,2) + pow(ratesyst,2)));

  printf("const double probGamma%sRate = %f;\n", linenumber, n_ec*ratemult);
  printf("const double probGamma%sRate_statup = %f;\n", linenumber, ratestat_up);
  printf("const double probGamma%sRate_statlo = %f;\n", linenumber, ratestat_lo);

  puts("");
}

double getpar(int i) // zero indexed, really!
{
  double val, err;
  mn->GetParameter(i, val, err);
  return val;
}

double geterr(int i) // zero indexed!
{
  double val, err;
  mn->GetParameter(i, val, err);
  return err;
}

/* Print the TFormula-style expression for a Gaussian, arranged so that
 * one of the parameters gives the integral, another the mean, and
 * the rest take care of the energy scale and resolution as is useful
 * for this fit. If arepar is false, instead of the strings being
 * parameters, they are constants. */
string gaus(const string integral, const string mean,
            const bool arepar = true)
{
  const string INT  = (arepar?"[":"") + integral + (arepar?"]":"");
  const string MEAN = (arepar?"[":"") + mean     + (arepar?"]":"");

  const string width =
    "sqrt([23]**2* " + MEAN + "+([24]*" + MEAN + ")**2+[25]**2)";

  return INT + "/(" + width + "*sqrt(2*TMath::Pi()))"
    "*exp(-(((x-" + MEAN + "*[1])/" + width + ")**2)/2)";
}

/* Print the TFormula-style expression for a Gaussian, arranged so that
 * the first parameter gives the integral, the second the mean, and
 * the third sigma (as opposed to gaus(0), which is norm/mean/sigma).
 * "starting" gives the number of the first parameter. Give, e.g.
 * gaus("[0]", "[1]", "[2]") */
string plaingaus(const string integral, const string mean,
                 const string sigma)
{
  return integral + "/(" + sigma + "*sqrt(2*TMath::Pi()))"
    "*exp(-(((x-" + mean + ")/" + sigma + ")**2)/2)";
}

TF1 * drawpeak(const char * const peak)
{
  const double a = gg->GetParameter("er_a");
  const double b = gg->GetParameter("er_b");
  const double c = gg->GetParameter("er_c");

  static int uniq = 0;

  TF1 * mypeak = new TF1(Form("mypeak%d", uniq++),
    "[2]/([1]*sqrt(2*TMath::Pi()))*"
    "exp(-(((x- [3]*[0])/[1])**2)/2)", 0, 15);
  mypeak->SetNpx(400);

  mypeak->SetParameter(0, gg->GetParameter(1));
  const double E = gg->GetParameter(Form("e%s", peak));
  mypeak->SetParameter(1, sqrt(a*a*E + b*b*E*E + c*c));

  mypeak->SetParameter(2,
     gg->GetParameter(gg->GetParNumber(Form("n%s", peak)))
     // zero if there is no b12_nX parameter
     -ggn->GetParameter(ggn->GetParNumber(Form("b12n_n%s", peak)))
     * (1-neff)/neff);

  mypeak->SetParameter(3,
    gg->GetParameter(1+gg->GetParNumber(Form("n%s", peak))));

  mypeak->Draw("same");
  return mypeak;
}

TF1 * drawnpeaks(TF1 * gg)
{
  TF1 * f = (TF1 *)gg->Clone("gdpeak");
  for(int i = 0; i < f->GetNpar(); i++){
    if(
       !strcmp(f->GetParName(i), "er_a") || 
       !strcmp(f->GetParName(i), "er_b") || 
       !strcmp(f->GetParName(i), "er_c") || 
       !strcmp(f->GetParName(i), "energyscale") || 
       !strcmp(f->GetParName(i), "n_totn") || 
       !strcmp(f->GetParName(i), "e_hn") || 
       !strcmp(f->GetParName(i), "frac_hn")) continue;
    f->SetParameter(i, 0);
  }
  f->Draw("same");
  return f;
}

void fcn(int & npar, double * gin, double & chi2, double *par, int flag)
{
  chi2 = 0;

  // B-12 everything
  for(int i = 0; i < gg->GetNpar(); i++)
    gg ->SetParameter(i, par[i]);

  const int ncorr = 4;
  const int corrections[ncorr] = {3, 5, 7, 11};
  const int corrpar[ncorr] = {31, 32, 33, 34};

  // Add the feed-down from B-12n to B-12 to the expected number of
  // B-12-gamma-like events.
  for(int i = 0; i < ncorr; i++)
    gg->SetParameter(corrections[i],
      // B-12 signal without accidental neutron
      (1-paccn)*par[corrections[i]]

      // B-12n signal with inefficiency and without accidental neutron
      + (1-par[2])*(1-paccn)*par[corrpar[i]]); // par[2] ~= neff

  // B-12n gamma probabilities
  for(int i = 0; i < ncorr; i++)
    ggn->SetParameter(i+2,

      // B-12n signal with efficiency and without accidental neutron 
      (1-paccn)*par[2]*par[corrpar[i]]

      // B-12n signal with accidental neutron
     + paccn*par[corrections[i]]
    );

  ggn->SetParameter(7, par[35]); // through-Hn-n-B12
  ggn->SetParameter(30, par[36]); // through-Gdn-n-B12

  // B-12n energy scale.
  ggn->SetParameter(1, par[1]);

  // bg normalization
  ggn->SetParameter(6, n_to_0_bg_rat);


  // accidental and Michel background in the B-12n sample
  // and resolution variables
  for(int i = 15; i <= 25; i++) ggn->SetParameter(i, par[i]);

#if 1
  static int call = 0;
  if(call++%100 == 0){
    gg->Draw("same"); ggn->Draw("same");
    c2->SetLogy(); c2->Update();  c2->Modified();
  }
#endif
  for(int i = 0; i < ehist->GetNbinsX(); i++){
    const double bincenter = ehist->GetBinCenter(i);
    const double binlo = ehist->GetBinLowEdge(i);
    const double binhi = ehist->GetBinLowEdge(i+1);
    if(bincenter < 0.7 || bincenter > 15) continue;
 
    {
      const double data = ehist->GetBinContent(i);

      const double theory = gg->Integral(binlo, binhi)/(binhi-binlo);

      chi2 += theory - data;
      if(data > 0 && theory > 0) chi2 += data*log(data/theory);
    }
 
    {
      const double data = ehistn->GetBinContent(i);

      const double theory = ggn->Integral(binlo, binhi)/(binhi-binlo);

      chi2 += theory - data;
      if(data > 0 && theory > 0) chi2 += data*log(data/theory);
    }
  }
  chi2 *= 2;

  // pull term for neutron efficiency
  chi2 += pow((par[2] - neff)/neff_e, 2);

  // pull term for "correlated background" (a.k.a. Li-8) component
  chi2 += pow((par[26] - corrbgnorm_nom)/corrbgnorm_err, 2);
}

// call the function at the current value for the side effects, e.g.
// setting the parameters of gg and ggn
void callfcn()
{
  double chi2, dummy;
  int npar = mn->GetNumPars();
  double par[mn->GetNumPars()];

  for(int i = 0; i < npar; i++) mn->GetParameter(i, par[i], dummy);

  fcn(npar, NULL, chi2, par, 0);
}

double getlimup(TMinuit * mn, int i)
{
  double answer, dum;
  int idum;
  TString sdum;
  mn->mnpout(i, sdum, dum, dum, dum, answer, idum);
  return answer;
}


void fixat(TMinuit * mn, int i, float v)
{
  mn->Command(Form("REL %d", i));
  if(getlimup(mn, i-1)) mn->Command(Form("SET LIM %d", i));
  mn->Command(Form("SET PAR %d %g", i, v));
  mn->Command(Form("FIX %d", i));
}

// Given two chi2 (or loglikelihoods*2), return the significance
// of their separation
double lratsig(const double l1, const double l2)
{
  if(l1 < l2)
    fprintf(stderr, "lratsig: Reporting reversed significance\n");

  const double dll = fabs(l1-l2)/2;

  // Too far out in the tail to do the quantile.
  if(dll*2 > 64) return sqrt(2*dll) - 0.290756;

  const double rat = exp(dll);
  return TMath::NormQuantile(1-1/(2*rat));
}

void findsigforarate(const int fnum)
{
  return; // sometimes disable because it is slow
  mn->Command("MINIMIZE");
  mn->Command("show par");
  const double with = mn->fAmin;
  fixat(mn, fnum, 0);
  mn->Command("MINIMIZE");
  mn->Command("show par");
  const double without = mn->fAmin;
  printtwice("%s significance = %f\n", 2, mn->fCpnam[fnum-1].Data(),
             lratsig(without, with));
  mn->Command(Form("REL %d", fnum));
  mn->Command(Form("SET LIM %d 0. 100.", fnum));
  mn->Command("MINIMIZE");
}

void setlinemarkercolor(TH1 * h, int c)
{
  h->SetLineColor(c);
  h->SetMarkerColor(c);
}

void b12gamma_finalfit(const int region = 1, const int whichcorr_ = 0, double targfrac = 0)
{
  if(region < 0){
    puts("You gave me a negative region, so I am just compiling");
    return;
  }

  whichcorr = whichcorr_;
  const double lowt   = region == 0? 4000 :region == 1?  3008: 2016;
  const double hight = 5024.; // ns
  const double highfq = 215; // MeV

  // From B-12-like event counting, with a bit over 1% stat error. Since
  // the most precise output of this fit has 10% stat errors, I won't
  // obsess over that.
  const double fq_eff = n_c12cap_forb12gamma/n_c12cap;

  const double b12eff = 1
    * light_noise_eff
    * mich_eff // Yes, this makes sense, even here. It
               // is the probability of not having an
               // accidental muon or a radiative capture.
    * sub_muon_eff05  // Subsequent muon veto efficiency 
    * eff_eor_b12 // timeleft cut
    * (b12ecutlow == 3?0.9251:
       b12ecutlow == 4?b12energyeff:
       (exit(1),1)) // B-12 energy cut
    * (distcut == 400?wholedet_dist400eff:(exit(1),1))
    * (exp(-b12lowt *log(2)/b12hl)
      -exp(-b12hight*log(2)/b12hl)) // B-12 beta decay time
    * fq_eff
  ;

  const double li8eff = 1
    * sub_muon_eff05  // Subsequent muon veto efficiency 
    * eff_eor_li8 // timeleft cut
    * (b12ecutlow == 3?0.9057:
       b12ecutlow == 4?0.8227:
       (exit(1),1)) // energy cut
    * (distcut == 400?wholedet_dist400eff:(exit(1),1))
    * (exp(-li8lowt *log(2)/li8hl)
      -exp(-li8hight*log(2)/li8hl)) // Li-8 beta decay time
    * fq_eff
  ;
  const double gammatimecut_eff = 
    (exp(-lowt/lifetime_c12)-exp(-hight/lifetime_c12));
  const double b12geff = b12eff * gammatimecut_eff;
  const double li8geff = li8eff * gammatimecut_eff;

  printtwice("Gamma time efficiency: %f%%\n", 1, 100*gammatimecut_eff);
  printtwice("Efficiency for B-12 gammas: %f%%\n", 1, 100*b12geff);
  printtwice("Efficiency for Li-8 gammas: %f%%\n", 1, 100*li8geff);

  const string gdndist = 
    // Gd-n distribution, normalized to 1, with a 
    // parameterization derived from the histogram in doc-5593-v3
    "("

    // normalization of the gaussian peak
    "0.760129/(sqrt([23]**2*7.983+([24]*7.983)**2+[25]**2 + 0.1091)*"
    "sqrt(2*TMath::Pi()))"

    // main gaussian.  The width is the quadrature sum of the intrinsic
    // width (sqrt(0.1091)) and the resolution
     "*exp(-(((x-7.983*[1])/"
     "sqrt([23]**2* 7.983+([24]*7.983)**2+[25]**2 + 0.1091))**2)/2)"

    // Then we have a parametrization for the low-energy tail
    // First a double logisitic that goes up from zero at 5MeV, then
    // down in the middle of the gaussian peak
    "+0.0372863/(1+exp(-(x-5.0757)/0.332407))/(1+exp((x-8)/0.1))"

    // Then a logistic that handles the stuff below 5, and goes
    // away in the middle of the gaussian peak.
    "+0.0149869/(1+exp((x-8)/0.1))"

    ")";

  gg = new TF1("gg", (
    "[0] +"
    // Six gaussians in which one of the fit parameters is the integral,
    // with a parameter for the energy scale, in for which the width is
    // set by the energy and three resolution parameters, all of which
    // are the same for all gaussians.
    + gaus( "3",  "4") + "+" + gaus( "5",  "6") + "+"
    + gaus( "7",  "8") + "+" + gaus( "9", "10") + "+"
    + gaus("11", "12") + "+" + gaus("13", "14") + "+"
    // bg
     "[15]+gaus(16)+gaus(19)+[22]*(3*(x/52.8)^2-2*(x/52.8)^3) +"
    + gaus("26", "27") + "+" + gaus("([30]*[28])", "[29]", false)

    // [30] is the normalization of the all n captures, and 28 is the
    // fraction that are Hn
    + " + ([30]*(1-[28]))*" + gdndist 
     ).c_str(), 0,15);

  ggn = new TF1("ggn", (gaus("[2]", Form("%.9f", G1E), false) + "+" + // sic
                        gaus("[3]", Form("%.9f", G2E), false) + "+" +
                        gaus("[4]", Form("%.9f", G3E), false) + "+" +
                        gaus("[5]", Form("%.9f", G4E), false) + "+" + 
                        gaus("[7]", "2.224573",false) +  // sic
      // accidental bg and michels.  Right for the accidentals,
      // wrong for the Michels, but I think it is ok.
      "+  [6]*("
       "[15]+gaus(16)+gaus(19)+[22]*(3*(x/52.8)^2-2*(x/52.8)^3)"
          ")"
      " + [30]*" + gdndist
  
      ).c_str(), 0, 15);

  ggn->SetParNames("unused0", "unused1", "b12n_n1",
                   "b12n_n2", "b12n_n3", "b12n_n5",
                   "bgnorm", "thrub12hnn"); // incomplete
  TFile * f = 
//#define FAST
#ifdef FAST
  new TFile(Form("b12gamma.sel%d.root", region), "read");
  targfrac = 0.181098;
  if(!f || f->IsZombie()) f =
#endif
  new TFile(rootfile3up, "read");
  TTree * t = (TTree *)f->Get("t");

  bg->Sumw2();
  bgn->Sumw2();
  corrbg->Sumw2();
  ehist->Sumw2();

  printf("Pre-selecting...\n");
  TTree * seltree = t->CopyTree(Form(
         "!earlymich && "
         "e > %f && e < 15 && "
         "dist < %f && "
         "fq < %f && "
         "micht >= %f && micht < %f"
         , b12ecutlow, distcut, fq_per_mev*highfq, lowt, hight));
  seltree->Write();

  printf("%lld events in t, %lld in seltree\n",
         t->GetEntries(), seltree->GetEntries());
  if(seltree->GetEntries() == 0){ printf("No events in seltree\n"); exit(1);}

  int ndecay, latennear;
  float dt, timeleft, miche, micht, fq; // yes, all floats
  #define SBA(x) seltree->SetBranchAddress(#x, &x);
  SBA(ndecay); 
  SBA(latennear); 
  SBA(dt); 
  SBA(timeleft); 
  SBA(miche); 
  SBA(micht); 
  SBA(fq); 

  printf("Drawing...\n");
  for(int i = 0; i < seltree->GetEntries(); i++){
    seltree->GetEntry(i);

    if(latennear == 0){
      if(timeleft > b12hight && dt > b12lowt &&
         dt < b12hight && ndecay == 0)
        ehist ->Fill(corrmiche(miche, fq/fq_per_mev));

      // "ndecay == 0" how to handle this? Really need the first event
      // in lots of windows, which isn't convenient -- I think this is
      // close enough
      if(timeleft > acchight && dt > acclowt && dt < acchight)
        bg->Fill(corrmiche(miche, fq/fq_per_mev));

      if(timeleft > li8hight && dt > li8lowt &&
         dt < li8hight && ndecay == 0)
        corrbg->Fill(corrmiche(miche, fq/fq_per_mev));
    }
    else if(latennear == 1){
      if(timeleft > b12hight && dt > b12lowt &&
         dt < b12hight && ndecay == 0)
        ehistn->Fill(corrmiche(miche, fq/fq_per_mev));

      if(timeleft > acchight && dt > acclowt && dt < acchight)
        bgn->Fill(corrmiche(miche, fq/fq_per_mev));

      if(timeleft > li8hight && dt > li8lowt &&
         dt < li8hight && ndecay == 0)
        corrbgn->Fill(corrmiche(miche, fq/fq_per_mev));
    }
  }

  if(ehist->GetEntries() == 0){ printf("No signal events\n"); exit(1); }

  setlinemarkercolor(bgn, kRed-7);
  setlinemarkercolor(bg, kRed);
  setlinemarkercolor(corrbg, kViolet);
  setlinemarkercolor(corrbgn, kViolet-8);
  setlinemarkercolor(ehistn, kOrange);

  const double bgscaleb12 = (b12hight - b12lowt)/(acchight-acclowt)
           *eff_eor_b12/eff_eor_acc;
  const double bgscaleli8 = (li8hight - li8lowt)/(acchight-acclowt)
           *eff_eor_li8/eff_eor_acc;

  n_to_0_bg_rat = bgn->Integral()/bg->Integral();

  ggn->SetParameter(6, n_to_0_bg_rat);

  bg->Scale(bgscaleb12);
  bgn->Scale(bgscaleb12);

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

  printf("Fitting accidental background...\n");
  bg->Fit("bggg", "li", "", 0.7, 40);

  // corrbg has accidentals in it too
  corrbg->Add(bggg, -bgscaleli8/bgscaleb12);

  // And scale the accidental-substracted version to the
  // expected amount in the signal window
  const double corrbgscale =
    (exp(-b12lowt/li8hl) - exp(-b12hight/li8hl))/
    (exp(-li8lowt/li8hl) - exp(-li8hight/li8hl))/eff_eor_li8;

  TF1 * corrbgfit = new TF1("corrbgfit",
    plaingaus("[0]", "[1]", "[2]").c_str(),
    0.7, 2);
  corrbgfit->SetLineColor(kViolet);
  const double li8gamma = 0.9808
     + 4.40 / (13.0*0.2/*NT*/ + 17.4*0.8/*GC*/) /* quenched alpha */;
  corrbgfit->FixParameter(1, li8gamma);

  // energy resolution
  corrbgfit->FixParameter(2, sqrt(pow(0.077, 2)*li8gamma +
                                  pow(0.018*li8gamma, 2) +
                                  pow(0.167, 2)));
  for(int i = 1; i <= corrbg->GetNbinsX(); i++)
    if(corrbg->GetBinContent(i) < 0)
      corrbg->SetBinError(i, 1);
  corrbg->Scale(corrbgscale);

  printf("Fitting correlated background...\n");
  corrbg->Fit("corrbgfit", "ie");
  gMinuit->Command("show min");
  corrbgfit->SetNpx(300);

  const char * ggpars[31] = { "accidentals", "energyscale", "neff",
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
  "bgmich",
  "er_a", "er_b", "er_c",
  "li8n", "li8m",
  "frac_hn", "e_hn",
  "n_totn",
  };

  const int NMINUITPAR = gg->GetNpar() + ggnnpars;
  mn = new TMinuit(NMINUITPAR);
  mn->SetPrintLevel(-1);
  mn->fGraphicsMode = false;
  mn->SetFCN(fcn);

  for(int i = 0; i < gg->GetNpar(); i++) gg->SetParName(i, ggpars[i]);

  {
    int err = 0;
    for(int mni = 0; mni < gg->GetNpar(); mni++)
      mn->mnparm(mni, gg->GetParName(mni), 0, 0.1, 0, 0, err);

    mn->mnparm(gg->GetNpar()+0, "b12n_n1", 0, 0.1, 0, 0, err);
    mn->mnparm(gg->GetNpar()+1, "b12n_n2", 0, 0.1, 0, 0, err);
    mn->mnparm(gg->GetNpar()+2, "b12n_n3", 0, 0.1, 0, 0, err);
    mn->mnparm(gg->GetNpar()+3, "b12n_n5", 0, 0.1, 0, 0, err); // sic
    mn->mnparm(gg->GetNpar()+4, "thrub12hnn",0,0.1, 0, 0, err);
    mn->mnparm(gg->GetNpar()+5, "thrub12gdnn",0,0.1,0, 0, err);
  }

  gg->SetNpx(400);
  ggn->SetNpx(400);
  ggn->SetLineColor(kOrange);

  for(int i = 15; i <= 22; i++)
    fixat(mn, 1+i, bggg->GetParameter(i-15));

  fixat(mn, 1+23, 0.077);
  fixat(mn, 1+24, 0.018);
  fixat(mn, 1+25, 0.1);

  fixat(mn, 1+0, 0);
  fixat(mn, 1+1, 1);
  mn->Command(Form("SET PAR 3 %f", neff));

  fixat(mn, 1+4, G1E);
  fixat(mn, 1+6, G2E);
  fixat(mn, 1+8, G3E);

  fixat(mn, 1+9, 0);
  fixat(mn, 1+10,0);

  fixat(mn, 1+12,G4E);
  fixat(mn, 1+13, 0);
  fixat(mn, 1+14,9.040);


  corrbgnorm_nom = max(0, corrbgfit->GetParameter(0));
  corrbgnorm_err = corrbgfit->GetParError(0);
  mn->Command(Form("SET PAR 27 %f", corrbgnorm_nom));
  mn->Command(Form("set limits 27 0. %f",
                   corrbgnorm_nom + 5*corrbgnorm_err));
  fixat(mn, 1+27, corrbgfit->GetParameter(1));
  fixat(mn, 1+28, corrbgfit->GetParameter(2));

  // contamination by Hn and Gdn events.
  {
    const double hn_t = 179.;

    // from my own measurements of B-12n and Li-8n from C-13
    const double Pn = 0.516 + 0.049 * b12hl/li8hl;
    const double ntrueb12n = Nc13cap*Pn;
    const double nobsb12n = ntrueb12n*b12eff;

    printf("Getting muon stop positions...\n");
    if(targfrac == 0) targfrac = 
      double(t->GetEntries(Form(
          "mx**2+my**2 < 1154**2 && "
          "abs(mz) < 1229+4+0.03*(1154-sqrt(mx**2+my**2)) && "
         "ndecay == 0 && fq < %f", fq_per_mev*highfq)))/
      t->GetEntries(Form(
         "ndecay == 0 && fq < %f", fq_per_mev*highfq));
    printf("neutron Target fraction: %f\n", targfrac);

    // Rough number for early Gd-h from my own study.
    // Will be less earlier and approach 0.875 (?) later.
    const double gdfrac = 0.78;

    // Hn in the GC, easy.  It's just the total times the probability
    // of an exponential during the time window
    const double n_hn_gc = (1-targfrac) * nobsb12n *
      (exp(-lowt/1000./hn_t)-exp(-hight/1000./hn_t));

    // Gd-n [in the Target].  Complicated, so integral a function
    // to get the probability.
    const double inwindowgdprob =
      gdtime()->Integral(lowt/1000, hight/1000)/gdtime()->Integral(0,600);
    const double n_gn = gdfrac * targfrac * nobsb12n *
      inwindowgdprob;
    printf("In-window Gd probability: %.3f\n", inwindowgdprob);


    // H-n in the Target.  It is obviously proportional to the Gd-n
    // by the fraction defined above.
    const double n_hn_targ = (1 - gdfrac)/gdfrac * n_gn;

    const double n_hn = n_hn_gc + n_hn_targ;

    printf("Fraction of muon stops in the target: %.3f\n", targfrac);
    printf("Number of Hn = %.1f (T %.3f, GC %.1f), Gdn = %.3f\n",
           n_hn, n_hn_targ, n_hn_gc, n_gn);

    fixat(mn, 1+gg->GetParNumber("frac_hn"), n_hn/(n_hn+n_gn));
    fixat(mn, 1+gg->GetParNumber("e_hn"), hn_e);

    mn->Command(Form("set limits %d 0. %f",
                     1+gg->GetParNumber("n_totn"),
                     10*(n_gn+n_hn)*ehist->GetBinWidth(1)));
    fixat(mn, 1+gg->GetParNumber("n_totn"),
           (n_gn+n_hn)*ehist->GetBinWidth(1));
  }

  // Limit the number of events to be positive
  for(int i = 3; i <= 13; i+=2)
    mn->Command(Form("set limits %d 0. 100.", i+1));
  for(int i = 31; i <= 36; i++)
    mn->Command(Form("set limits %d 0. 100.", i+1));

  mn->Command("show par");

  ehist->Draw("e");
  ehist->GetXaxis()->SetRangeUser(0.7, 10);
  ehist->GetYaxis()->SetRangeUser(0.0001, 30);


  ehistn->Draw("samee");
  bg->Draw("histsame");
  bgn->Draw("histsame");
  corrbg->Draw("samee");
  corrbgn->Draw("samee");
  bggg->Draw("same");
  corrbgfit->Draw("same");

  printf("Fitting signal, step 1/2...\n");
  mn->Command("SIMPLEX");
  mn->Command("MIGRAD");
  mn->Command("Release 26"); // er_c
  mn->Command("Release 2"); // energy scale
  mn->Command(Form("Set limits 26 0 %f", 1.92530e-01 + 1.81760e-02));
  if(ehist->GetEntries() < 50) fixat(mn, 26, 0.186);

  printf("Fitting signal, step 2/2...\n");
  mn->Command("Set strategy 2");
  mn->Command("MIGRAD");
  mn->Command("HESSE"); // we use the covariance matrix
                        // for the ground state calculation

  if(region == 0 || region == 1){
    double errmat[NMINUITPAR][NMINUITPAR];
    mn->mnemat(&errmat[0][0], NMINUITPAR);

    // translate gamma lines to MINUIT *free* parameters
    const int l2fp[5] = { -1, 2, 3, 4, 5 };
    // translate gamma lines to MINUIT *parameters*
    const int l2p[5] = { -1, 3, 5, 7, 11 };
    for(int i = 1; i <= 4; i++){
      for(int j = i+1; j <= 4; j++){
        // Use the conditions under which we fit lines 1 and 2 to give
        // us the correlation coefficient for them. Use the conditions
        // under which we fit lines 3 and 4 to give us all the rest,
        // including those involving lines 1 and 2.
        if((region == 1) ^ (i == 1 && j == 2))
          printf("const double b12cc_%d%d = %g;\n", i, j,
            errmat[l2fp[i]][l2fp[j]]/geterr(l2p[i])/geterr(l2p[j]));
      }
    }
  }

  printf("Finding MINOS errors...\n");
  if(region == 0){

    // Can't let the energy scale float when the 2621 line is zeroed,
    // since it mostly sets the energy scale. For consistency, keep it
    // fixed for the other tests too (but not for MINOS).
    mn->Command("fix 2");
    findsigforarate(4);
    findsigforarate(6);
    findsigforarate(32);
    findsigforarate(33);
    mn->Command("release 2");

    mn->Command("MINOS 10000 4");
    mn->Command("MINOS 10000 6");
    mn->Command("MINOS 10000 32");
    mn->Command("MINOS 10000 33");
  }
  else{
    mn->Command("fix 2");
    findsigforarate(8);
    findsigforarate(12);
    findsigforarate(34);
    findsigforarate(35);
    mn->Command("release 2");

    mn->Command("MINOS 10000 8");
    mn->Command("MINOS 10000 12");
    mn->Command("MINOS 10000 34");
    mn->Command("MINOS 10000 35");
  }

  callfcn(); // get gg and ggn set right

  mn->Command("show min");

  TF1 * nh_peak = drawpeak("_hn");
  nh_peak->SetLineStyle(7);

  TF1 * gh_peak = drawnpeaks(gg);
  gh_peak->SetLineStyle(7);

  if(region == 0) {
    {
      drawpeak("1");

      // The number of events in the peak, corrected for efficiency.
      const double nev = gg->GetParameter("n1")/ehist->GetBinWidth(1);
      const double neveup = mn->fErp[2]/gg->GetParameter("n1")*nev;
      const double nevelo = mn->fErn[2]/gg->GetParameter("n1")*nev;

      print_results(b12geff, 953, nev, nevelo, neveup, 
        sqrt(pow(b12lineEsyst[0], 2) + pow(b12lineNsyst[0], 2)));
    }
    {
      drawpeak("2");

      const double nev = gg->GetParameter("n2")/ehist->GetBinWidth(1);
      const double neveup = mn->fErp[3]/gg->GetParameter("n2")*nev;
      const double nevelo = mn->fErn[3]/gg->GetParameter("n2")*nev;

      print_results(b12geff, 1674, nev, nevelo, neveup,
        sqrt(pow(b12lineEsyst[1], 2) + pow(b12lineNsyst[1], 2)));
    }
    {
      const double nev = getpar(31)/ehist->GetBinWidth(1);
      const double neveup = mn->fErp[8]/getpar(31)*nev;
      const double nevelo = mn->fErn[8]/getpar(31)*nev;

      print_results13(neff*b12geff, 953, nev, nevelo, neveup, 0);
    }
    {
      const double nev = getpar(32)/ehist->GetBinWidth(1);
      const double neveup = mn->fErp[9]/getpar(32)*nev;
      const double nevelo = mn->fErn[9]/getpar(32)*nev;

      print_results13(neff*b12geff, 1674, nev, nevelo, neveup, 0);
    }
    {
      const double nev = getpar(26)/ehist->GetBinWidth(1)/corrbgscale;
      const double neveup = geterr(26)/getpar(26)*nev;
      const double nevelo = geterr(26)/getpar(26)*nev;

      print_results8(li8geff, 981, nev, nevelo, neveup, 0);
    }
  }
  else if(region == 1){
    {
      drawpeak("3");

      const double nev = gg->GetParameter("n3")/ehist->GetBinWidth(1);
      const double neveup = mn->fErp[4]/gg->GetParameter("n3")*nev;
      const double nevelo = mn->fErn[4]/gg->GetParameter("n3")*nev;

      print_results(b12geff, 2621, nev, nevelo, neveup,
        sqrt(pow(b12lineEsyst[2], 2) + pow(b12lineNsyst[2], 2)));
    }
    {
      drawpeak("5");

      const double nev = gg->GetParameter("n5")/ehist->GetBinWidth(1);
      const double neveup = mn->fErp[5]/gg->GetParameter("n5")*nev;
      const double nevelo = mn->fErn[5]/gg->GetParameter("n5")*nev;

      print_results(b12geff, 3759, nev, nevelo, neveup,
        sqrt(pow(b12lineEsyst[3], 2) + pow(b12lineNsyst[3], 2)));
    }
    {
      const double nev = getpar(33)/ehist->GetBinWidth(1);
      const double neveup = mn->fErp[10]/getpar(33)*nev;
      const double nevelo = mn->fErn[10]/getpar(33)*nev;

      print_results13(neff*b12geff, 2621, nev, nevelo, neveup, 0);
    }
    {
      const double nev = getpar(34)/ehist->GetBinWidth(1);
      const double neveup = mn->fErp[11]/getpar(34)*nev;
      const double nevelo = mn->fErn[11]/getpar(34)*nev;

      print_results13(neff*b12geff, 3759, nev, nevelo, neveup, 0);
    }
  }
  else if(region == 2 && gg->GetParameter("n6") != 0){
    drawpeak("6");

    const double nev = gg->GetParameter("n6")/ehist->GetBinWidth(1);
    const double neveup = mn->fErp[6]/gg->GetParameter("n6")*nev;
    const double nevelo = mn->fErn[6]/gg->GetParameter("n6")*nev;

    puts("THE FOLLOWING WILL BE WRONG UNLESS YOU FIX THE CODE");
    
    print_results(b12geff, 9040, nev, nevelo, neveup, 0);
  }
  gg->Draw("same");

  mn->Command("show cov");
}
