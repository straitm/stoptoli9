#include <fstream>
#include "consts.h"
#include "TFile.h"
#include "TMinuit.h"
#include "TTree.h"
#include "TF1.h"
#include "TH1.h"
#include "TROOT.h"
#include "TRandom3.h"
#include <stdio.h>
#include <vector>
#include <algorithm>
using std::vector;

// Usually set to zero to get the capture rates, but set to 1 
// to do the C-13 -> B-12+n analysis
static const int NNEUTRON = 0;

static void printfr(const char * const msg, ...)
{
  va_list ap;
  va_start(ap, msg);
  printf(RED);
  vprintf(msg, ap);
  printf(CLR);
}

int reactorpowerbin(const int run)
{
  static bool inited = false;
  static vector<int> on_off, off_off;
  if(!inited){
    inited = true;
    // From doc-5095-v2. I see that most are included in the on-off
    ifstream offofffile("offoff.h");
    // From doc-5341
    ifstream onofffile ("onoff.h");
    int r;
    while(offofffile >> r) off_off.push_back(r);
    while(onofffile  >> r) on_off.push_back(r);
  }

  if(std::binary_search(off_off.begin(), off_off.end(), run)) return 0;
  if(std::binary_search( on_off.begin(),  on_off.end(), run)) return 1;
  return 2;
}

static const double energyeff = 0.83;  // B-12 energy cut
static const double othereff = 0.977 // previous muons
                             * 0.981 // subsequent muons
                             * 0.9709; // time until end of run
static const double errorecut = 0.02;
static const double eff = othereff * energyeff;
static const double ferr_energy = errorecut/energyeff;

static const double capprob12 = (37.9/(37.9 + 1/2197e-6));
static const double ferrcapprob12 = 5./379.;

static const double capprob13 = (35.0/(35.0 + 1/2197e-6));
static const double ferrcapprob13 = 2./35.0;

TRandom3 ran;

// Weighted average of T and GC measurements
static const double f13 = 0.010918;

// Conservatively assign a range 10x the quoted (gaussian?) error,
// i.e. about 3x the quoted error if it is gaussian. This makes *NO*
// difference since the other errors are so much larger.
static const double f13_hi = 0.010928, f13_lo = 0.010908;

// ratio of probability of capture on C-13 / C-12
static const double caprat = capprob13 / capprob12;
static const double ferrorp13op12 =
  sqrt(pow(ferrcapprob13,2)+pow(ferrcapprob12,2));

// Convert between observed B-12 decays and muon captures given
// some parameters
double conversion(const bool nominal = false)
{
  const double e0 = 0.99983;

  // only valid for full detector analysis
  const double e1 = neff_dr_800_avg*neff_dt_avg;
  
  // For the (whole) target
  //const double e1 = neff_dr_800_targ*neff_dt_targ;

  const double p1212 = 0.18602; // NOM
  const double p1212_hi = 0.1981, p1212_lo = 0.1739;


  // My guess for C-13->B-13. There is no data other than my own, which
  // only constrains it to < 85%. I am putting in a guess somewhat
  // higher on average than for B-12->B-12 since the neutron separation
  // energy of B-13 is larger. Measday says that the (pi, gamma) reaction
  // feeds the ground state quite strongly, supporting this.
  const double p1313 = 0.2;
  const double p1313_lo = 0.1, p1313_hi = 0.3;

  // my own measurement for C-13->B-12, which circularly depends on the
  // results gotten here, so since I don't have joint fit set up, must
  // iterate a little. Fortunately, this is not a dominant error, so it
  // pretty much doesn't matter.
  const double p1312 = 0.4738;
  const double p1312_err = sqrt(0.0247*0.0247 + 0.0336*0.0336);

  double caprat_now = caprat + ran.Gaus(0, ferrorp13op12*caprat);
  if(caprat_now < 0) caprat_now = 0;
  double f13_now   = ran.Rndm()*(ran.Rndm() > 0.5?
                           (f13_hi - f13   ) + f13:
                           (f13    - f13_lo) + f13_lo);

  // NOTE: nearly all the spread in the result comes from this one
  // This one is not relevant for the B-12 gamma analysis
  double p1212_now = ran.Rndm()*(p1212_hi-p1212_lo) + p1212_lo;

  double p1313_now = ran.Rndm()*(p1313_hi-p1313_lo) + p1313_lo;
  double p1312_now = ran.Gaus(0, p1312_err) + p1312;

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
  return rms(vals)/conversion(true);
}

const int nbin = 4000;
 
double all(TTree * t, TF1 * ee, const char * const addcut)
{
  t->Draw(Form("dt/1000 - 2.028e-6 >> h(%d, 0.001, 100)", nbin),
          Form("(%s) && "
          "latennear == %d && timeleft > 100e3 && miche < 12 && "
          "e > 4 && e < 15 && !earlymich", addcut, NNEUTRON), "e");
  TH1D * h = (TH1D*)gROOT->FindObject("h");
  h->Fit("ee", "lq");
  h->Fit("ee", "li");
  gMinuit->Command("MINOS 10000 1");
  gMinuit->Command("SHOW MIN");
  double err = (-gMinuit->fErn[0]+gMinuit->fErp[0])/2;
  if(err == 0){
    h->Fit("ee", "li");
    gMinuit->Command("MINOS 20000");
    err = (-gMinuit->fErn[0]+gMinuit->fErp[0])/2;
    if(err == 0){
      err = ee->GetParError(0);
      fprintf(stderr, "warning, fell back to parabolic error\n");
    }
  }

  ee->SetParError(0, err);

  const double ferrorfit = ee->GetParError(0)/ee->GetParameter(0);
  const double rawintegral = ee->GetParameter(0)*0.0202/log(2)/h->GetBinWidth(1);

  const double integral_pd_oec = rawintegral/livetime/eff;

  const double fsyst = sqrt(pow(fsysterr(),2) + pow(ferr_energy,2));

  printf("stuff with b12 lifetime raw %f +- %f\n",
         rawintegral, ferrorfit * rawintegral);

  printf("stuff with b12 lifetime raw with overall eff "
         "corrected\n\t%f +- %f (just the fit error)\n",
         integral_pd_oec*livetime, ferrorfit * integral_pd_oec*livetime);

  printf("stuff with b12 lifetime raw per day with overall eff "
         "corrected\n\t%f +- %f (just the fit error)\n",
         integral_pd_oec, ferrorfit * integral_pd_oec);

  printf("fractional stat err from fit: %f\n", ferrorfit);
  printf("fractional systematic error from B-13 alone: %.2f%%\n", fsysterr()*100);
  printf("fractional systematic energy cut error:      %.2f%%\n", ferr_energy*100);
  printf("Quadrature sum of those, FWIW:               %.2f%%\n", fsyst*100);
  printf("fractional systematic neutron eff error:     %.2f%%\n", f_neff_dt_error*100);

  if(NNEUTRON == 1){
    const double f_b12nsysterr = sqrt(pow(ferr_energy, 2)+pow(f_neff_dt_error, 2));
    printf("Since you are looking at B-12+n production,\nyour fractional "
      "systematic error is %f\n", f_b12nsysterr);

    printf("And the rate of B-12 like events with a neutron, corrected "
      "for\nefficiencies, is %f +- %f +- %f\n\n",
       integral_pd_oec/neff_dt_avg/neff_dr_800_avg,
       integral_pd_oec/neff_dt_avg/neff_dr_800_avg * ferrorfit,
       integral_pd_oec/neff_dt_avg/neff_dr_800_avg * f_b12nsysterr);

    printf("And now correct for accidentals: %f +- %f +- %f\n\n",
       integral_pd_oec/neff_dt_avg/neff_dr_800_avg - 0.007,
       integral_pd_oec/neff_dt_avg/neff_dr_800_avg * ferrorfit,
       integral_pd_oec/neff_dt_avg/neff_dr_800_avg * f_b12nsysterr);

    const double f_b12n_finalsyst = sqrt(pow(n_c13cap_err/n_c13cap, 2)
      +pow(f_neff_dt_error, 2)
      -pow(ferr_energy, 2)); // cancels, so remove from n_c13cap error!

    printf("For the probablility per capture, the B-12 energy cut error\n"
      "*cancels*, but there is the error on the C-13 denominator, so\n"
      "the fractional systematic is %f\n\n", f_b12n_finalsyst);

    printf("And so the probability per capture is %f +- %f +- %f\n\n", 
       (integral_pd_oec/neff_dt_avg/neff_dr_800_avg - 0.007)/n_c13cap,
       (integral_pd_oec/neff_dt_avg/neff_dr_800_avg * ferrorfit)/n_c13cap,
       (integral_pd_oec/neff_dt_avg/neff_dr_800_avg * f_b12n_finalsyst)/n_c13cap);
  }

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

  TFile *_file0 = TFile::Open(rootfile3up);
  TTree * t = (TTree *)_file0->Get("t");
  TF1 * ee = new TF1("ee",
     "[0]*exp(-x*log(2)/0.0202) + " // B-12
     "[1]*exp(-x*log(2)/0.8399) + " // Li-8
     "[2]*exp(-x*log(2)/[3]) + "    // N-16
     "[4]*exp(-x*log(2)/[5]) + "    // Li-9/He-8
     "[6]*exp(-x*log(2)/[7]) + "    // Something long-lived
     "[8]", 0, 100);
  ee->SetParameter(0, 1e5);

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

  // 
  ee->FixParameter(4, 0);
  ee->FixParameter(5, 15);
  ee->FixParameter(6, 0);
  ee->FixParameter(7, 15);
  // 

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
