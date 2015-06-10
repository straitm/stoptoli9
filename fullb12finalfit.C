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
const double b13life = 17.33/log(2.);
const double li8life = 839.9/log(2.);
const double n16life = 7130./log(2.);

const double b12life_err = 0.02/log(2.);
const double b13life_err = 0.17/log(2.);
const double li8life_err = 0.9/log(2.);
const double n16life_err = 20./log(2.);

struct ev{
  double t; // time
  int n; // number of neutrons
  ev(double t_, int n_){
    t = t_;
    n = n_;
  }
};

TMinuit * mn = NULL;
TTree * selt = NULL;

vector<ev> events;

void printfr(const char * const msg, ...)
{
  va_list ap;
  va_start(ap, msg);
  printf(RED);
  vprintf(msg, ap);
  printf(CLR);
}

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



const double mum_count =   7.16322e5;
const double mum_count_e = 0.05388e5;



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

const double c12nuc_cap = c_atomic_capture_prob*(1-f13)*mum_count*capprob12;
const double c13nuc_cap = c_atomic_capture_prob*f13    *mum_count*capprob13;



// Will subtract mean muon lifetime, 2028ns, and mean transit time for
// light from B-12, 12ns. Doesn't make a real difference.
const double offset = 2.040e-6;

const double lowtime = 1.0 - offset;
const double hightime = 100e3;
const double totaltime = hightime - lowtime;

const double b12energyeff = 0.8494;  // B-12 energy cut
const double b12energyeff_e = 0.0063;

const double b13energyeff = b12energyeff * 1.015; // estimate from MC
const double b13energyeff_e = 0.02; // BS

const double li8eff_energy = 0.71; // Not as careful as the boron ones.
const double li8eff_energy_e = 0.01;

// time until end of run
const double eor_eff = 1-(1-0.9709)*hightime/100e3;

const double mich_eff = 0.9996;

const double b12eff = mich_eff * eor_eff * b12energyeff;
const double b12ferr_energy = b12energyeff_e/b12energyeff;

const double b13eff = mich_eff * eor_eff * b13energyeff;
const double b13ferr_energy = b13energyeff_e/b13energyeff;

const double li8eff = mich_eff * eor_eff * li8eff_energy;

// Using no distance cut ("nlate") so as not to have to figure out
// what the efficiency for that is for the high purity sample.
const double neff = neff_dt_highpurity;

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

void fcn(int & npar, double * gin, double & like, double *par, int flag)
{
  const double p_b12 = par[0],
               p_b12n= par[1],
               p_b13 = par[2],
               p_li8 = par[3],
               p_li8n= par[4];

  const double n_b12 = p_b12 *c12nuc_cap*b12eff,
               n_b12n= p_b12n*c13nuc_cap*b12eff,
               n_b13 = p_b13 *c13nuc_cap*b13eff,
               n_li8 = p_li8 *c12nuc_cap*li8eff,
               n_li8n= p_li8n*c13nuc_cap*li8eff,
               n_n16 = par[5],
               acc   = par[6],
               accn  = par[7],
               b12t  = par[8],
               b13t  = par[9],
               li8t  = par[10],
               n16t  = par[11];

  like =    n_b12 *exp(-lowtime/b12t)
          + n_b12n*exp(-lowtime/b12t)
          + n_b13 *exp(-lowtime/b13t)
          + n_li8 *exp(-lowtime/li8t)
          + n_li8n*exp(-lowtime/li8t)
          + n_n16 *exp(-lowtime/n16t)
          + totaltime*(acc+accn);

  printf(".");  fflush(stdout);

  for(unsigned int i = 0; i < events.size(); i++){
    const double mt = -events[i].t;
    const double n = events[i].n;

    double f;
    if(n == 0)
      f =          n_b12 /b12t*exp(mt/b12t) +
                   n_b13 /b13t*exp(mt/b13t) +
                   n_li8 /li8t*exp(mt/li8t) +
         (1-neff)*(n_b12n/b12t*exp(mt/b12t) +
                   n_li8n/li8t*exp(mt/li8t)) +
                   n_n16 /n16t*exp(mt/n16t) +
          acc;
    else if(n == 1)
      f = neff*(n_b12n/b12t*exp(mt/b12t) +
                n_li8n/li8t*exp(mt/li8t)) +
          accn;
    else
      continue;

    if(f > 0) like -= log(f);
  }

  like *= 2;

  // pull terms
  like += pow((b12t - b12life)/b12life_err, 2)
        + pow((b13t - b13life)/b13life_err, 2)
        + pow((li8t - li8life)/li8life_err, 2)
        + pow((n16t - n16life)/n16life_err, 2);

  if(p_b12n + p_b13 + p_li8n > 1)
    like += pow((p_b12n + p_b13 + p_li8n - 1)/errcapprob13*capprob13, 2);
}

static double getpar(TMinuit * mn, int i)
{
  double answer, dum;
  mn->GetParameter(i, answer, dum);
  return answer;
}

void fullb12finalfit(const char * const cut =
"mx**2+my**2 < 1050**2 && mz > -1175 && "
"abs(fez + 62*ivdedx/2 - 8847.2) < 1000 && chi2 < 2 && "
"timeleft > %f && miche < 12 && !earlymich && "
"e > 4 && e < 15 && dt < %f")
{
  printf("%sB-12 selection efficiency is %.1f%%%s\n", RED, b12eff*100, CLR);

  TFile *_file0 = TFile::Open(rootfile3up);
  TTree * t = (TTree *)_file0->Get("t");
 
  const int npar = 12;
  mn = new TMinuit(npar);
  mn->SetPrintLevel(-1);
  mn->fGraphicsMode = false;
  //mn->Command("SET STRATEGY 2");
  mn->SetFCN(fcn);
  int err;
  mn->mnparm(0, "p_b12", 0.2,  1e2, 0, 1, err);
  mn->mnparm(1, "p_b12n",0.5,  1e2, 0, 1, err);
  mn->mnparm(2, "p_b13", 0.2,  1e2, 0, 1, err);
  mn->mnparm(3, "p_li8", 0.01,  1e1, 0, 1, err);
  mn->mnparm(4, "p_li8n",0.05,  1e1, 0, 1, err);
  mn->mnparm(5, "n_n16",   1,  1e1, 0, 1e3, err);
  mn->mnparm(6, "acc",   10 ,    1, 0, 1e3, err);
  mn->mnparm(7, "accn",  10 ,    1, 0, 1e3, err);
  mn->mnparm(8, "b12t", b12life,  b12life_err, 0, 0, err);
  mn->mnparm(9, "b13t", b13life,  b13life_err, 0, 0, err);
  mn->mnparm(10, "li8t", li8life,  li8life_err, 0, 0, err);
  mn->mnparm(11, "n16t", n16life,  n16life_err, 0, 0, err);

  // XXX kill N-16 for faster fitting.
  mn->Command("SET PAR 6 0");
  mn->Command("FIX 6");
  mn->Command("FIX 12");

  printf("Making cuts...\n");
  TFile tmpfile("/tmp/b12tmp.root", "recreate");
  selt = t->CopyTree(Form(cut, hightime+offset, hightime+offset));
  events.clear();

  float dt;
  int nn;
  selt->SetBranchAddress("dt", &dt);
  selt->SetBranchAddress("laten", &nn);

  printf("Filling in data array...\n");
  for(int i = 0; i < selt->GetEntries(); i++){
    selt->GetEntry(i);
    events.push_back(ev(dt-offset, nn));
  }

  printf("MIGRAD");
  mn->Command("MIGRAD");
  puts(""); mn->Command("show par");

  for(int i = 1; i <= 3; i++){
    printf(Form("MINOS %d", i));
    mn->Command(Form("MINOS 10000 %d", i));
    puts("");
    mn->Command("show min");
  }

  // I did the fit in terms of probability per capture so that it was
  // straightforward to impose a unitarity bound, but this was only
  // a constant factor shift, now shift back to doing it in terms of
  // counts so I can step through the uncertainties.

  const double fit = c12nuc_cap*b12eff*getpar(mn, 0);
  const double ferrorfitup = mn->fErp[0]/getpar(mn, 0);
  const double ferrorfitlo = mn->fErn[0]/getpar(mn, 0);

  printtwice("\nC-12 -> B-12 raw %f +%f %f\n", 0,
    fit, ferrorfitup*fit, ferrorfitlo*fit);

  const double b12like_central = fit/b12eff/mum_count * 100;
  const double staterrup = ferrorfitup*b12like_central;
  const double staterrlo = ferrorfitlo*b12like_central;
  const double muerr = mum_count_e/mum_count*b12like_central;
  const double b12err = b12ferr_energy*b12like_central;
  const double toterrup = sqrt(pow(ferrorfitup,2) +
                               pow(mum_count_e/mum_count,2) +
                               pow(b12ferr_energy, 2))*b12like_central;
  const double toterrlo = -sqrt(pow(ferrorfitlo,2) +
                               pow(mum_count_e/mum_count,2) +
                               pow(b12ferr_energy, 2))*b12like_central;
  printtwice("\nC-12 -> B-12, eff corrected, percent per mu- stop\n"
         "%f +%f%f(fit) +-%f(mu count) +-%f(B-12 eff), +%f%f(total)\n",
         3, b12like_central, staterrup, staterrlo, muerr, b12err,
         toterrup, toterrlo);

  const double b12like_central_percap = b12like_central/capprob12;
  const double staterr_percapup = staterrup/capprob12;
  const double staterr_percaplo = staterrlo/capprob12;
  const double muerr_percap = muerr/capprob12;
  const double b12err_percap = b12err/capprob12;
  const double capfracerr_percap =
    b12like_central_percap * err_capprob/capprob12;
  const double toterr_percapup = sqrt(pow(staterr_percapup,2)+
                                    pow(muerr_percap,2)+
                                    pow(b12err_percap,2)+
                                    pow(capfracerr_percap,2));
  const double toterr_percaplo = -sqrt(pow(staterr_percaplo,2)+
                                    pow(muerr_percap,2)+
                                    pow(b12err_percap,2)+
                                    pow(capfracerr_percap,2));

  printtwice("\nOr percent per C-12 mu- capture\n"
         "%f +%f%f(fit) +-%f(mu count) +-%f(B-12 eff) +-%f(cap frac), "
         "+%f%f(total)\n", 2, 
         b12like_central_percap, staterr_percapup, staterr_percaplo,
         muerr_percap, b12err_percap, capfracerr_percap,
         toterr_percapup, toterr_percaplo);

  const double b12like_central_rate = b12like_central/lifetime_c12/100;
  const double staterr_rateup = staterrup/lifetime_c12/100;
  const double staterr_ratelo = staterrlo/lifetime_c12/100;
  const double muerr_rate = muerr/lifetime_c12/100;
  const double b12err_rate = b12err/lifetime_c12/100;
  const double lifetimeerr_rate = b12like_central_rate
    * lifetime_c12_err/lifetime_c12;
  const double toterr_rateup = sqrt(pow(staterr_rateup,2)+
                                  pow(muerr_rate,2)+
                                  pow(b12err_rate,2)+
                                  pow(lifetimeerr_rate,2));
  const double toterr_ratelo = -sqrt(pow(staterr_ratelo,2)+
                                  pow(muerr_rate,2)+
                                  pow(b12err_rate,2)+
                                  pow(lifetimeerr_rate,2));

  printtwice("\nOr 10^3/s: %f +%f%f(fit) +-%f(mu count) +-%f(B-12 eff),\n"
         "+-%f(lifetime) +%f%f(total)\n", 2,
         b12like_central_rate, staterr_rateup, staterr_ratelo,
         muerr_rate, b12err_rate, lifetimeerr_rate, toterr_rateup,
         toterr_ratelo);
  puts("");
}
