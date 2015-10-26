#include <fstream>
#include "consts.h"
#include "TFile.h"
#include "TMinuit.h"
#include "TTree.h"
#include "TF1.h"
#include "TH1.h"
#include "TROOT.h"
#include <stdio.h>
#include <vector>
#include <algorithm>
using std::vector;

const double b12life = 20.20/log(2.);
const double li8life = 839.9/log(2.);
const double n16life = 7130./log(2.);


// Weighted average of T and GC measurements for the HP region
const double f13 = 0.010921;

const double mulife = 2196.9811e-6;

const double c_atomic_capture_prob = 0.998;
const double c_atomic_capture_prob_err = 0.001;

const double lifetime_c12 = 2028.e-6;
const double lifetime_c13 = 2037.e-6;
const double lifetime_c12_err = 2.e-6;
const double lifetime_c13_err = 8.e-6;

const double lifetime_c = lifetime_c12*(1-f13)+lifetime_c13*f13;
const double lifetime_c_err = sqrt(pow(lifetime_c12_err*(1-f13),2)
                                  +pow(lifetime_c13_err*   f13 ,2));

const double capprob12 = 1-lifetime_c12/mulife;
const double errcapprob12 = (1-(lifetime_c12+lifetime_c12_err)/mulife)/2
                           -(1-(lifetime_c12-lifetime_c12_err)/mulife)/2;

const double capprob13 = 1-lifetime_c13/mulife;
const double errcapprob13 = (1-(lifetime_c13+lifetime_c13_err)/mulife)/2
                           -(1-(lifetime_c13-lifetime_c13_err)/mulife)/2;

const double capprob = c_atomic_capture_prob *
                      (capprob12*(1-f13) + capprob13*f13);

const double err_capprob = sqrt(pow(errcapprob12,2)*(1-f13)
                              + pow(errcapprob13,2)*f13 +
  pow(c_atomic_capture_prob_err/c_atomic_capture_prob * capprob, 2));

struct ev{
  double t; // time
  ev(double t_){
    t = t_;
  }
};

int reactorpowerbin(const int run)
{
  bool inited = false;
  vector<int> on_off, off_off;
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

TMinuit * mn = NULL;
TTree * selt = NULL;

vector<ev> events;

static void printfr(const char * const msg, ...)
{
  va_list ap;
  va_start(ap, msg);
  printf(RED);
  vprintf(msg, ap);
  printf(CLR);
}

/*
 * Prints the message once with the requested precision and in RED, then
 * again with all digits in the default color, starting with the first
 * number.
 *
 * The message must only have floating point substitutions.
 */
static void printtwice(const char * const msg, const int digits, ...)
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
        *pmp++ = '%';
        *bmp++ = '%';
        *pmp++ = '.';
        *pmp++ = digits+'0';
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

/**********************************************************************/
#include "mucountfinalfit.C"

// The number of stopping mu-.
const double mum_count = mucountfinalfit_cut(
  "ndecay == 0 && mx**2+my**2 < 1050**2 && mz > -1175 && "
  "abs(fez + 62*ivdedx/2 - 8847.2) < 1000 && rchi2 < 2", false);
const double mum_count_e = mucountfinalfit_cut(
  "ndecay == 0 && mx**2+my**2 < 1050**2 && mz > -1175 && "
  "abs(fez + 62*ivdedx/2 - 8847.2) < 1000 && rchi2 < 2", true);

// The number of mu- stopping and entering atomic orbitals of
// any isotope of carbon.  Subtly different from the number of
// stopping mu-!
const double mumc_count = c_atomic_capture_prob*mum_count;
const double mumc_count_e =  c_atomic_capture_prob*
sqrt(
  pow(mum_count_e, 2)
  + mumc_count*pow(c_atomic_capture_prob_err/c_atomic_capture_prob,2)
);
/**********************************************************************/

// Will subtract mean muon lifetime, 2028ns, and mean transit time for
// light from B-12, 12ns. Doesn't make a real difference.
const double offset = 2.040e-6;

const double lowtime = 1.0 - offset;
const double hightime = 100e3;
const double totaltime = hightime - lowtime;

static const double energyeff = 0.8504;  // B-12 energy cut
static const double energyeff_e = 0.0065;

// time until end of run
static const double eor_eff = 1-(1-0.9709)*hightime/100e3;

// Subsequent muon veto efficiency (an efficiency on the isotope decay,
// NOT on the muon), for the hard cut imposed on events in order to get
// into the ntuples of 0.5ms.
double sub_muon_eff = 0; // set by main function argument

static double eff = 0;
static const double ferr_energy = energyeff_e/energyeff;

void fcn(int & npar, double * gin, double & like, double *par, int flag)
{
  like = 2*(par[0]*exp(-lowtime/b12life) + par[1]*exp(-lowtime/li8life)
          + par[2]*exp(-lowtime/n16life) + totaltime*par[3]);

  printf(".");  fflush(stdout);

  for(unsigned int i = 0; i < events.size(); i++){
    const double mt = -events[i].t;
    const double f = par[0]/b12life*exp(mt/b12life) +
                     par[1]/li8life*exp(mt/li8life) +
                     par[2]/n16life*exp(mt/n16life) +
                     par[3];
    if(f > 0) like -= 2*log(f);
  }
}

static double getpar(TMinuit * mn, int i)
{
  double answer, dum;
  mn->GetParameter(i, answer, dum);
  return answer;
}

struct ve{
  double val, err;
};

ve b12like_finalfit(const char * const cut =
"mx**2+my**2 < 1050**2 && mz > -1175 && "
"abs(fez + 62*ivdedx/2 - 8847.2) < 1000 && rchi2 < 2 && "
"timeleft > %f && miche < 12 && !earlymich && "
"e > 4 && e < 15 && dt < %f",
const bool verbose = true,
const double sub_muon_eff_in = sub_muon_eff05,
const double mylivetime = -1.0)
{
  if(verbose){
    printtwice("TECHNOTE 3.5: Number of mu- stops is %f +- %f\n", 0,
               mum_count, mum_count_e);
    printtwice("TECHNOTE 3.5: Number of mu- atomic captures, any C "
               "isotope is %f +- %f\n", 0, mumc_count, mumc_count_e);
  }

  sub_muon_eff = sub_muon_eff_in;
  eff = light_noise_eff * mich_eff * eor_eff * sub_muon_eff * energyeff;
  if(verbose)
    printtwice("TECHNOTE 4.1.2: B-12 selection efficiency: %f +- %f %%\n",
      2, eff*100, ferr_energy*eff*100);

  // XXX kludgy!
  TFile *_file0 = TFile::Open(mylivetime<0?rootfile3up:rootfile3up_extended);

  TTree * t = (TTree *)_file0->Get("t");
 
  const int npar = 4;
  mn = new TMinuit(npar);
  mn->SetPrintLevel(-1);
  mn->Command("SET STRATEGY 2");
  mn->SetFCN(fcn);
  int err;
  mn->mnparm(0, "n_b12", 1e5,  1e2, 0, 1e7, err);
  mn->mnparm(1, "n_li8", 1e3,  1e1, 0, 1e4, err);
  mn->mnparm(2, "n_n16", 1e3,  1e1, 0, 1e3, err);
  mn->mnparm(3, "acc",   10 ,    1, 0, 1e3, err);
  printf("Making cuts...\n");
  TFile tmpfile("/tmp/b12tmp.root", "recreate");
  selt = t->CopyTree(Form(cut, hightime+offset, hightime+offset));
  events.clear();

  float dt;
  selt->SetBranchAddress("dt", &dt);

  printf("Filling in data array...\n");
  for(int i = 0; i < selt->GetEntries(); i++){
    selt->GetEntry(i);
    events.push_back(ev(dt-offset));
  }

  printf("MIGRAD");
  mn->Command("MIGRAD");
  puts(""); mn->Command("show par");
  printf("MINOS");
  mn->Command("MINOS 10000 1");
  puts(""); mn->Command("show min");

  const double ferrorfit = (mn->fErp[0]-mn->fErn[0])/2/getpar(mn, 0);

  // The "stuff with b12 lifetime raw" result is the one we need
  // to find the ratio between high-purity and loose selections
  // that is used to find the denominator for the loose selection.
  ve result_for_ratio;
  result_for_ratio.val = getpar(mn, 0);
  result_for_ratio.err = ferrorfit * getpar(mn, 0);

  // Always print since this is needed for both HP and anti-HP samples
  printtwice("\nTECHNOTE 5.1: Stuff with b12 lifetime raw %f +- %f\n", 0,
             result_for_ratio.val, result_for_ratio.err);


  if(mylivetime > 0){
    const double b12like_central = 1000*getpar(mn, 0)/eff/mylivetime;
    const double staterr = ferrorfit*b12like_central;
    const double muerr = mumc_count_e/mumc_count*b12like_central;
    const double b12err = ferr_energy*b12like_central;
    const double toterr = sqrt(pow(ferrorfit,2) +
                                 pow(mumc_count_e/mumc_count,2) +
                                 pow(ferr_energy, 2))*b12like_central;

    printtwice("\nStuff with b12 lifetime raw with overall eff "
         "corrected, per live kilosecond\n"
         "%f +-%f(fit) +-%f(mu count) +-%f(B-12 eff), %f(total)\n",
         3, b12like_central, staterr, muerr, b12err, toterr);
  }

  const double b12like_central = getpar(mn, 0)/eff/mumc_count * 100;
  const double staterr = ferrorfit*b12like_central;
  const double muerr = mumc_count_e/mumc_count*b12like_central;
  const double b12err = ferr_energy*b12like_central;
  const double toterr = sqrt(pow(ferrorfit,2) +
                               pow(mumc_count_e/mumc_count,2) +
                               pow(ferr_energy, 2))*b12like_central;

  if(verbose)
    printtwice("\nTECHNOTE 4.2: Stuff with b12 lifetime raw with overall eff "
         "corrected, percent per mu- stop: "
         "%f +-%f(fit) +-%f(mu count) +-%f(B-12 eff), %f(total)\n",
         3, b12like_central, staterr, muerr, b12err, toterr);


  const double b12like_central_percap = b12like_central/capprob;
  const double staterr_percap = staterr/capprob;
  const double muerr_percap = muerr/capprob;
  const double b12err_percap = b12err/capprob;
  const double capfracerr_percap =
    b12like_central_percap * err_capprob/capprob;
  const double toterr_percap = sqrt(pow(staterr_percap,2)+
                                    pow(muerr_percap,2)+
                                    pow(b12err_percap,2)+
                                    pow(capfracerr_percap,2));
  if(verbose)
    printtwice("\nTECHNOTE 4.2: Or percent per mu- capture: "
         "%f +-%f(fit) +-%f(mu count) +-%f(B-12 eff) +-%f(cap frac), "
         "%f(total)\n", 2, 
         b12like_central_percap, staterr_percap, muerr_percap,
         b12err_percap, capfracerr_percap, toterr_percap);

  const double b12like_central_rate = b12like_central/lifetime_c/100;
  const double staterr_rate = staterr/lifetime_c/100;
  const double muerr_rate = muerr/lifetime_c/100;
  const double b12err_rate = b12err/lifetime_c/100;
  const double lifetimeerr_rate = b12like_central_rate
    * lifetime_c_err/lifetime_c;
  const double toterr_rate = sqrt(pow(staterr_rate,2)+
                                  pow(muerr_rate,2)+
                                  pow(b12err_rate,2)+
                                  pow(lifetimeerr_rate,2));

  if(verbose)
    printtwice("\nTECHNOTE 4.2: Or 10^3/s: %f +- %f(fit) "
         "+- %f(mu count) +- %f(B-12 eff), "
         "+- %f(lifetime) %f(total)\n", 2,
         b12like_central_rate, staterr_rate, muerr_rate,
         b12err_rate, lifetimeerr_rate, toterr_rate);
  puts("");

  return result_for_ratio;
}

void loosecaptures_finalfit()
{
  const char * const HPcut =
    "mx**2+my**2 < 1050**2 && mz > -1175 && "
    "abs(fez + 62*ivdedx/2 - 8847.2) < 1000 && rchi2 < 2";

  const char * const othercuts =
    "timeleft > %f && miche < 12 && !earlymich && "
    "e > 4 && e < 15 && dt < %f";

  const ve hpresult =
    b12like_finalfit(Form("  %s  && %s", HPcut, othercuts), true);

  const ve antihpresult =
    b12like_finalfit(Form("!(%s) && %s", HPcut, othercuts), false);

  const double rat = antihpresult.val/hpresult.val;
  const double answer = (1+rat)*mumc_count;
  const double error = sqrt(
    pow(rat*mumc_count/hpresult.val *     hpresult.err,2) +
    pow(    mumc_count/hpresult.val * antihpresult.err,2)
  );

  puts("The following go into consts.h, AND into section 5.1+5.2");
  puts("of the tech note, AND into dcfluids.ods");
  printtwice("TECHNOTE 5.1: Atomic carbon captures in the loose sample: %g +- %g\n",
    4, answer,
    answer * sqrt(
      pow(error/answer,2)
      +pow(mumc_count_e/mumc_count,2)
    )
  );

  printf("TECHNOTE 5.2: Error due to B-12-like statistics: %.1f%%\n",
         error/answer * 100);

  puts("");

  printtwice("TECHNOTE 5.2: Atomic captures/day on C-12: %f +- %f\n",
    0, answer/livetime*(1-f13),
        error/livetime*(1-f13));

  const double n_c12cap = answer/livetime*(1-f13) * capprob12;
  const double n_c12cap_err = 
    answer*sqrt(
      pow(error/answer,2) // fractional stat error
      +pow(errcapprob12/capprob12,2) // fractional capture error
      +pow(mumc_count_e/mumc_count,2)
    )/livetime*(1-f13)*capprob12;

  printtwice("TECHNOTE 5.2: *Nuclear* captures/day on C-12: %f +- %f\n",
    1, n_c12cap, n_c12cap_err);

  printf("const double n_c12cap     = %f\n", n_c12cap);
  printf("const double n_c12cap_err = %f\n", n_c12cap_err);

  puts("");

  printtwice("TECHNOTE 5.2: Atomic captures/day on C-13: %f +- %f\n",
    2, answer/livetime*f13, error/livetime*f13);

  const double n_c13cap = answer/livetime*f13 * capprob13;
  const double n_c13cap_err = answer*sqrt(
      pow(error/answer,2) // fractional stat error
      +pow(errcapprob13/capprob13,2) // fractional capture error
      +pow(mumc_count_e/mumc_count,2)
    )/livetime*f13*capprob13;

  printtwice("TECHNOTE 5.2: *Nuclear* captures/day on C-13: %f +- %f\n",
    2, n_c13cap, n_c13cap_err);

  printf("const double n_c13cap     = %f\n", n_c13cap);
  printf("const double n_c13cap_err = %f\n", n_c13cap_err);

  puts("");
  puts("");

  puts("These numbers go into section 6.1 in the tech note AND ");
  puts("consts.h AND dcfluids.ods.");
  puts("They don't depend on the above fits, but I'm putting them");
  puts("in this order in the output to mirror the technote.");
  printtwice("TECHNOTE 6.1: Atomic captures/day on C-12: %f +- %f\n",
    0, mumc_count/livetime*(1-f13),
    mumc_count_e/livetime*(1-f13));

  const double n_c12cap_hp = mumc_count/livetime*(1-f13) * capprob12;
  const double n_c12cap_hp_err = mumc_count*sqrt(
      pow(mumc_count_e/mumc_count,2) // fractional stat error
      +pow(errcapprob12/capprob12,2) // fractional capture error
    )/livetime*(1-f13)*capprob12;

  printtwice("TECHNOTE 6.1: *Nuclear* captures/day on C-12: %f +- %f\n",
    1, n_c12cap_hp, n_c12cap_hp_err);

  printf("const double n_c12cap_hp     = %f\n", n_c12cap_hp);
  printf("const double n_c12cap_hp_err = %f\n", n_c12cap_hp_err);

  puts("");

  printtwice("TECHNOTE 6.1: Atomic captures/day on C-13: %f +- %f\n",
    2, mumc_count/livetime*f13, mumc_count_e/livetime*f13);

  const double n_c13cap_hp = mumc_count/livetime*f13 * capprob13;
  const double n_c13cap_hp_err = mumc_count*sqrt(
      pow(mumc_count_e/mumc_count,2) // fractional stat error
      +pow(errcapprob13/capprob13,2) // fractional capture error
    )/livetime*f13*capprob13;

  printtwice("TECHNOTE 6.1: *Nuclear* captures/day on C-13: %f +- %f\n",
    3, n_c13cap_hp, n_c13cap_hp_err);

  printf("const double n_c13cap_hp     = %f\n", n_c13cap_hp);
  printf("const double n_c13cap_hp_err = %f\n", n_c13cap_hp_err);
}
