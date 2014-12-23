#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TH1.h"
#include "TROOT.h"
#include <stdio.h>

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

double eff = 0.82  // B-12 energy cut
           * 0.977 // previous muons
           * 0.962 // subsequent muons
           * 0.9709; // time until end of run
double ferrorecut = 0.01;
double livedays = 489.509;

double capprob12 = (37.9/(37.9 + 1/2197e-6));
double ferrcapprob12 = 5./379.;

double capprob13 = (35.0/(35.0 + 1/2197e-6));
double ferrcapprob13 = 2./35.0;

double nefftarg = 0.63*0.97;
// innermost first
double neffsmaller[3] = {0.5685*0.97, 0.6437*0.97, 0.6634*0.97 };

// Convert between observed B-12 decays and muon captures given
// some parameters
double conversion(const bool isef13hi, const bool isp1212hi,
                  const bool isp1313hi, const bool isp1312hi,
                  const bool nomoverride = false)
{
  const double e0 = 0.99983;
  const double e1 = 0.8162*0.97;

  const double p1212 = 0.18602; // NOM
  const double p1212_hi = 0.1981, p1212_lo = 0.1739;

  // ratio of probability of capture on C-13 / C-12
  const double caprat = capprob13 / capprob12;
  const double ferrorp13op12=sqrt(pow(ferrcapprob13,2)+pow(ferrcapprob12,2));

  const double f13 = 0.0107;
  const double f13_hi = 0.01147, f13_lo = 0.00963;

  // the gaussianized f13 error
  const double gf13err = (f13_hi-f13_lo)/sqrt(12);

  // the gaussian f13 error added in quadrature with the capture ratio
  // error to produce the "effective f13 fraction error"
  const double gef13err = sqrt(pow(gf13err,2) + pow(ferrorp13op12*f13,2));
  const double ef13 = f13 * caprat; // NOM

  // converted to a range.
  const double ef13_hi = ef13 + sqrt(12)/2 * gef13err,
               ef13_lo = ef13 - sqrt(12)/2 * gef13err;

  // my guesses for the 13->13 and 13-12 reactions
  const double p1313 = 0.2, p1312 = 0.5;
  const double p1313_lo = 0.1, p1313_hi = 0.3;
  const double p1312_lo = 0.3, p1312_hi = 0.7;

  double ef13_now = isef13hi?ef13_hi:ef13_lo;
  double p1212_now = isp1212hi?p1212_hi:p1212_lo;
  double p1313_now = isp1313hi?p1313_hi:p1313_lo;
  double p1312_now = isp1312hi?p1312_hi:p1312_lo;

  if(nomoverride){
    ef13_now = ef13;
    p1212_now = p1212;
    p1313_now = p1313;
    p1312_now = p1312;
  }

  const double dem = (1-ef13_now)*p1212_now*e0
                   +   ef13_now  *p1313_now*e0
                   +   ef13_now  *p1312_now*(1-e1);
  return 1/dem;
}

// Find the fractional systematic error associated with the uncertainties
// on the fraction of C-13, the probability of C-13 -> B-13 and the 
// probability of C-13 -> B-12
static double fsysterr()
{
  double min = 1e38, max = 0;
  for(int i = 0; i < 0x10; i++){
    const double val = conversion(i & 0x01, i & 0x02, i & 0x04, i & 0x08);
    if(val > max) max = val;
    if(val < min) min = val;
  }

  const double range = max - min;
  return range/sqrt(12)/conversion(0,0,0,0,1);
}
 
void all(TTree * t, TF1 * ee)
{
  t->Draw("dt/1000 - 2.028e-6 >> h(1000, 0.001, 100)",
          "latennear == 0 && timeleft > 100e3 && miche < 12 && "
          "e > 4 && e < 14.5 && !earlymich", "e");
  TH1D * h = (TH1D*)gROOT->FindObject("h");
  h->Fit("ee", "li");
  TF1 e("e", "[0]*exp(-x*log(2)/0.0202)", 0, 100);
  e.SetParameter(0, ee->GetParameter(0));

  const double ferrorfit = ee->GetParError(0)/ee->GetParameter(0);
  const double rawintegral = e.Integral(0, 10)/h->GetBinWidth(1);

  double ferror = sqrt(pow(ferrorfit, 2) + pow(ferrorecut, 2));
  printf("b12 raw %f +- %f\n", rawintegral, ferrorfit * rawintegral);

  const double integral_pd_oec = rawintegral/livedays/eff;
  printf("b12 raw per day with overall eff corrected\n\t%f +- %f\n",
         integral_pd_oec, ferrorfit * integral_pd_oec);

  const double fsyst = fsysterr();
  printf("fractional stat err from fit: %f\n", ferrorfit);
  printf("fractional systematic error:  %f\n", fsyst);

  const double totferr = sqrt(pow(ferrorfit,2)+pow(fsyst,2));

  const double caprate = integral_pd_oec * conversion(0,0,0,0,1);

  printfr("Capture rate: %f +- %f (stat) +- %f (syst)\n+- %f total\n",
         caprate, caprate * ferrorfit, caprate * fsyst,
         caprate * totferr);

  // Haven't yet put in efficiency errors

        
}

void targ(const int nn, TTree * t, TF1 * ee)
{
  printf("\nIn the target only:\n");
  t->Draw("dt/1000 - 2.028e-6 >> h(10000, 0.001, 100)",
          Form("latennear == %d && dx**2 + dy**2 < 1150**2 && "
          "abs(dz) < 1229+0.03*(1150-sqrt(dx**2+dy**2)) &&"
          "timeleft > 100e3 && miche < 12 && "
          "e > 4 && e < 14.5 && !earlymich", nn), "e");
  TH1D * h = (TH1D*)gROOT->FindObject("h");
  h->Fit("ee", "l");
  TF1 e("e", "[0]*exp(-x*log(2)/0.0202)", 0, 100);
  e.SetParameter(0, ee->GetParameter(0));

  const double ferrorfit = ee->GetParError(0)/ee->GetParameter(0);
  const double rawintegral = e.Integral(0, 10)/h->GetBinWidth(1);

  double ferror = sqrt(pow(ferrorfit, 2) + pow(ferrorecut, 2));
  printf("b12 raw %f +- %f\n", rawintegral, ferrorfit * rawintegral);
}

void innertarg(const int nn, TTree * t, TF1 * ee, const double limr,
               const double limz)
{
  printf("\nIn the inner %.1fx%.1f target only:\n", limr, limz);
  t->Draw(Form("dt/1000 - 2.028e-6 >> h%.0f(10000, 0.001, 100)", limr),
          Form("latennear == %d && dx**2 + dy**2 < %f**2 && "
          "abs(dz) < %f &&  timeleft > 100e3 && miche < 12 && "
          "e > 4 && e < 14.5 && !earlymich", nn, limr, limz), "e");
  TH1D * h = (TH1D*)gROOT->FindObject(Form("h%.0f", limr));
  h->Fit("ee", "l");
  TF1 e("e", "[0]*exp(-x*log(2)/0.0202)", 0, 100);
  e.SetParameter(0, ee->GetParameter(0));

  const double ferrorfit = ee->GetParError(0)/ee->GetParameter(0);
  const double rawintegral = e.Integral(0, 10)/h->GetBinWidth(1);

  double ferror = sqrt(pow(ferrorfit, 2) + pow(ferrorecut, 2));
  printf("b12 raw %f +- %f\n", rawintegral, ferrorfit * rawintegral);
}

void b12finalfit()
{
  TFile *_file0 = TFile::Open(
    "/cp/s4/strait/fullfido-300s-3-25MeV-20141117.root");
  TTree * t = (TTree *)_file0->Get("t");
  TF1 * ee = new TF1("ee", "[0]*exp(-x*log(2)/0.0202) + "
     "[1]*exp(-x*log(2)/0.8399) + "
     "[2]*exp(-x*log(2)/[3]) + "
     "[4]*exp(-x*log(2)/[5]) + "
     "[6]*exp(-x*log(2)/[7]) + "
     "[8]", 0, 100);
  ee->SetParLimits(1, 0, 100);

  ee->SetParLimits(2, 0, 300);
  ee->SetParameter(3, 7.13);
  ee->SetParLimits(3, 1.5, 10);

  ee->SetParLimits(4, 0, 300);
  ee->SetParameter(5, 15);
  ee->SetParLimits(5, 10, 20);

  ee->SetParLimits(6, 0, 300);
  ee->SetParameter(7, 0.2);
  ee->SetParLimits(7, 0.05, 0.4);
  ee->FixParameter(6, 0);
  ee->FixParameter(7, 0);

  all(t, ee);
  return;
  ee->FixParameter(4, 0);
  ee->FixParameter(5, 15);
  targ(0, t, ee);

/*
  innertarg(0, t, ee, 1045, 1127);
  ee->FixParameter(2, 0);
  ee->FixParameter(3, 7.13);
  innertarg(0, t, ee, 913, 985);
  innertarg(0, t, ee, 724, 781);
*/
}
