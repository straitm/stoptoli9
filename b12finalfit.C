#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TH1.h"
#include "TROOT.h"
#include "TRandom3.h"
#include <stdio.h>
#include <vector>
using std::vector;

const char * const RED     = "\033[31;1m"; // bold red
const char * const CLR      = "\033[m"    ; // clear

static void printfr(const char * const msg, ...)
{
  va_list ap;
  va_start(ap, msg);
  printf(RED);
  vprintf(msg, ap);
  printf(CLR);
}

static const double energyeff = 0.83;  // B-12 energy cut
static const double othereff = 0.977 // previous muons
                             * 0.981 // subsequent muons
                             * 0.9709; // time until end of run
static const double errorecut = 0.02;
static const double eff = othereff * energyeff;
static const double ferr_energy = errorecut/energyeff;
static const double livedays = 489.509;

static const double capprob12 = (37.9/(37.9 + 1/2197e-6));
static const double ferrcapprob12 = 5./379.;

static const double capprob13 = (35.0/(35.0 + 1/2197e-6));
static const double ferrcapprob13 = 2./35.0;

static const double nefftarg = 0.63*0.97;
// innermost first
static const double neffsmaller[3] = {0.5685*0.97, 0.6437*0.97, 0.6634*0.97 };

TRandom3 ran;

static const double f13 = 0.0107;
static const double f13_hi = 0.01147, f13_lo = 0.00963;

// ratio of probability of capture on C-13 / C-12
static const double caprat = capprob13 / capprob12;
static const double ferrorp13op12 =
  sqrt(pow(ferrcapprob13,2)+pow(ferrcapprob12,2));

// Convert between observed B-12 decays and muon captures given
// some parameters
double conversion(const bool nominal = false)
{
  const double e0 = 0.99983;
  const double e1 = 0.8162*0.97;

  const double p1212 = 0.18602; // NOM
  const double p1212_hi = 0.1981, p1212_lo = 0.1739;


  // my guesses for the 13->13 and 13-12 reactions
  const double p1313 = 0.2, p1312 = 0.5;
  const double p1313_lo = 0.1, p1313_hi = 0.3;
  const double p1312_lo = 0.3, p1312_hi = 0.7;

  double caprat_now = caprat + ran.Gaus(0, ferrorp13op12*caprat);
  if(caprat_now < 0) caprat_now = 0;
  double f13_now   = ran.Rndm()*(ran.Rndm() > 0.5?
                                (f13_hi - f13   ) + f13:
                                (f13    - f13_lo) + f13_lo);
  double p1212_now = ran.Rndm()*(p1212_hi-p1212_lo) + p1212_lo;
  double p1313_now = ran.Rndm()*(p1313_hi-p1212_lo) + p1313_lo;
  double p1312_now = ran.Rndm()*(p1312_hi-p1312_lo) + p1312_lo;

  // If nominal, throw all that out and set to nominal values.
  // (const-correct, what's that?)
  if(nominal){
    f13_now = f13;
    caprat_now = caprat;
    p1212_now = p1212;
    p1313_now = p1313;
    p1312_now = p1312;
  }


  // The first three lines take care of the very slight underestimation 
  // of the amount of B-12 caused by the shorter livetime of B-13 
  const double dem = 1/(1 + (0.2020 - 0.1733)/0.2020 * 
                          f13_now * p1313_now / ( (1-f13_now)*p1212_now +
                                                    f13_now  *p1312_now))
                  *
  // The rest pulls the B-12 number out given the fractions and efficiencies
                  ( (1-caprat_now*f13_now)*p1212_now*e0
                   +   caprat_now*f13_now  *p1313_now*e0
                   +   caprat_now*f13_now  *p1312_now*(1-e1));
  return 1/dem;
}

static double mean(const vector<double> & vals)
{
  if(vals.empty()) return 0;
  double tot = 0;
  for(unsigned int i = 0; i < vals.size(); i++) tot += vals[i];
  return tot/vals.size();
}

static double rms(const vector<double> & vals)
{
  if(vals.empty()) return 0;

  const double m = mean(vals);
  double sum = 0;
  for(unsigned int i = 0; i < vals.size(); i++)
    sum += pow(vals[i] - m,2);
  return sqrt(sum / vals.size());
  
}

// Find the fractional systematic error associated with the uncertainties
// on the fraction of C-13, the probability of C-13 -> B-13 and the 
// probability of C-13 -> B-12
static double fsysterr()
{
  vector<double> vals;
  for(int i = 0; i < 10000; i++) vals.push_back(conversion(false));
  rms(vals);

  return rms(vals)/conversion(true);
}

const int nbin = 1000;
 
double all(TTree * t, TF1 * ee, const char * const addcut)
{
  t->Draw(Form("dt/1000 - 2.028e-6 >> h(%d, 0.001, 100)", nbin),
          Form("(%s) && "
          "latennear == 0 && timeleft > 100e3 && miche < 12 && "
          "e > 4 && e < 15 && !earlymich", addcut), "e");
  TH1D * h = (TH1D*)gROOT->FindObject("h");
  h->Fit("ee", "lq");
  h->Fit("ee", "li");
  TF1 e("e", "[0]*exp(-x*log(2)/0.0202)", 0, 100);
  e.SetParameter(0, ee->GetParameter(0));

  const double ferrorfit = ee->GetParError(0)/ee->GetParameter(0);
  const double rawintegral = e.Integral(0, 10)/h->GetBinWidth(1);

  printf("b12 raw %f +- %f\n", rawintegral, ferrorfit * rawintegral);

  const double integral_pd_oec = rawintegral/livedays/eff;
  printf("b12 raw per day with overall eff corrected\n\t%f +- %f\n",
         integral_pd_oec, ferrorfit * integral_pd_oec);

  const double fsyst = sqrt(pow(fsysterr(),2) + pow(ferr_energy,2));
  printf("fractional stat err from fit: %f\n", ferrorfit);
  printf("fractional systematic error:  %f\n", fsyst);

  const double totferr = sqrt(pow(ferrorfit,2)+pow(fsyst,2));

  const double caprate = integral_pd_oec * conversion(true);

  printfr("Capture rate for both C-12 and C-13:\n"
          "%f +- %f (stat) +- %f (syst) --- +- %f total\n"
          "Total fractional error: %.2f%%\n",
         caprate,
         caprate * ferrorfit,
         caprate * fsyst,
         caprate * totferr,
         totferr*100);

  const double c13capfrac = 
   f13*capprob13/( f13*capprob13 + (1-f13)*capprob12 );
  
  const double c13capfrac_err = 
   sqrt(pow((f13_hi-f13_lo)/sqrt(12)/f13,2) + pow(ferrorp13op12,2));

  const double fsyst13 = sqrt(pow(fsyst,2)+pow(c13capfrac_err,2));
  const double totferr13 = sqrt(pow(totferr,2)+pow(c13capfrac_err,2));
  
  printfr("\nCapture rate for C-13:\n"
          "%f +- %f (stat) +- %f (syst) --- +- %f total\n"
          "Total fractional error: %.2f%%\n",
          c13capfrac*caprate,
          c13capfrac*caprate * ferrorfit,
          c13capfrac*caprate * fsyst13,
          c13capfrac*caprate * totferr13,
          totferr13*100);

  const double c12capfrac = 
   (1-f13)*capprob12/( f13*capprob13 + (1-f13)*capprob12 );
  
  const double c12capfrac_err = 
   ((1-f13_hi)-(1-f13_lo))/sqrt(12)/(1-f13);

  const double fsyst12 = sqrt(pow(fsyst,2)+pow(c12capfrac_err,2));
  const double totferr12 = sqrt(pow(totferr,2)+pow(c12capfrac_err,2));
  
  printfr("\nCapture rate for C-12:\n"
          "%f +- %f (stat) +- %f (syst) --- +- %f total\n"
          "Total fractional error: %.2f%%\n",
          c12capfrac*caprate,
          c12capfrac*caprate * ferrorfit,
          c12capfrac*caprate * fsyst12,
          c12capfrac*caprate * totferr12,
          totferr12*100);


  const double capprob = capprob12*(1-f13) + capprob13*f13;
  const double fsystmu = sqrt(pow(fsyst,2)+pow(ferrcapprob12,2));
  const double totferrmu = sqrt(pow(ferrorfit,2)+pow(fsystmu,2));

  printfr("Mu- stop rate:\n"
          "%f +- %f (stat) +- %f (syst) --- +- %f total\n"
          "Total fractional error: %.2f%%\n",
         caprate/capprob,
         caprate/capprob * ferrorfit,
         caprate/capprob * fsystmu,
         caprate/capprob * totferrmu,
         totferrmu*100);

  return rawintegral;
}

void targ(const int nn, TTree * t, TF1 * ee, const char * const addcut,
          const double norm)
{
  printf("\nIn the target only:\n");
  t->Draw(Form("dt/1000 - 2.028e-6 >> h(%d, 0.001, 100)", nbin),
          Form("(%s) && latennear == %d && !(dx**2 + dy**2 < 1150**2 && "
          "abs(dz) < 1229+0.03*(1150-sqrt(dx**2+dy**2))) &&"
          "timeleft > 100e3 && miche < 12 && "
          "e > 4 && e < 15 && !earlymich", addcut, nn), "e");
  TH1D * h = (TH1D*)gROOT->FindObject("h");
  h->Fit("ee", "li");
  TF1 e("e", "[0]*exp(-x*log(2)/0.0202)", 0, 100);
  e.SetParameter(0, ee->GetParameter(0));

  const double ferrorfit = ee->GetParError(0)/ee->GetParameter(0);
  const double rawintegral = e.Integral(0, 10)/h->GetBinWidth(1);

  printfr("b12 relative to whole %f +- %f\n", rawintegral/norm,
         ferrorfit/norm * rawintegral);
}

void innertarg(const int nn, TTree * t, TF1 * ee, const double limr,
               const double limz, const double norm)
{
  printf("\nIn the inner %.1fx%.1f target only:\n", limr, limz);
  t->Draw(Form("dt/1000 - 2.028e-6 >> h%.0f(%d, 0.001, 100)", limr, nbin),
          Form("latennear == %d && dx**2 + dy**2 < %f**2 && "
          "abs(dz) < %f &&  timeleft > 100e3 && miche < 12 && "
          "e > 4 && e < 15 && !earlymich", nn, limr, limz), "e");
  TH1D * h = (TH1D*)gROOT->FindObject(Form("h%.0f", limr));
  h->Fit("ee", "li");
  TF1 e("e", "[0]*exp(-x*log(2)/0.0202)", 0, 100);
  e.SetParameter(0, ee->GetParameter(0));

  const double ferrorfit = ee->GetParError(0)/ee->GetParameter(0);
  const double rawintegral = e.Integral(0, 10)/h->GetBinWidth(1);

  printfr("b12 relative to whole %f +- %f\n", rawintegral/norm,
          ferrorfit * rawintegral/norm);
}

void b12finalfit(const char * const addcut = "1")
{
  printf("%sB-12 selection efficiency is %.1f%%%s\n", RED, eff*100, CLR);

  TFile *_file0 = TFile::Open(
    "/cp/s4/strait/fullfido-300s-3-25MeV-20150219.root");
  TTree * t = (TTree *)_file0->Get("t");
  TF1 * ee = new TF1("ee", "[0]*exp(-x*log(2)/0.0202) + "
     "[1]*exp(-x*log(2)/0.8399) + "
     "[2]*exp(-x*log(2)/[3]) + "
     "[4]*exp(-x*log(2)/[5]) + "
     "[6]*exp(-x*log(2)/[7]) + "
     "[8]", 0, 100);
  ee->SetParLimits(1, 0, 1e5/nbin);

  ee->SetParLimits(2, 0, 3e5/nbin);
  ee->SetParameter(3, 7.13);
  ee->SetParLimits(3, 1.5, 10);

  ee->SetParLimits(4, 0, 3e5/nbin);
  ee->SetParameter(5, 15);
  ee->SetParLimits(5, 10, 20);

  ee->SetParLimits(6, 0, 3e5/nbin);
  ee->SetParameter(7, 0.2);
  ee->SetParLimits(7, 0.05, 0.4);
  ee->FixParameter(6, 0);
  ee->FixParameter(7, 0);

  const double norm = all(t, ee, addcut);
return; // <------
  ee->FixParameter(4, 0);
  ee->FixParameter(5, 15);
  targ(0, t, ee, addcut, norm);

  innertarg(0, t, ee, 1045, 1127, norm);
  ee->FixParameter(2, 0);
  ee->FixParameter(3, 7.13);
  innertarg(0, t, ee, 913, 985, norm);
  innertarg(0, t, ee, 724, 781, norm);
}
