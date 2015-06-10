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

static const double mum_count =   7.149e5;
static const double mum_count_e = 0.054e5;

static const double capprob12 = (37.9/(37.9 + 1/2197e-6));
static const double ferrcapprob12 = 5./379.;

static const double capprob13 = (35.0/(35.0 + 1/2197e-6));
static const double ferrcapprob13 = 2./35.0;

// Weighted average of T and GC measurements
static const double f13 = 0.010918;

// ratio of probability of capture on C-13 / C-12
static const double caprat = capprob13 / capprob12;
static const double ferrorp13op12 =
  sqrt(pow(ferrcapprob13,2)+pow(ferrcapprob12,2));

// Will subtract mean muon lifetime, 2028ns, and mean transit time for
// light from B-12, 12ns. Doesn't make a real difference.
const double offset = 2.040e-6;

const double lowtime = 1.0 - offset;
const double hightime = 100e3;
const double totaltime = hightime - lowtime;

static const double b12energyeff = 0.8494;  // B-12 energy cut
static const double b12energyeff_e = 0.0063;

static const double b13energyeff = b12energyeff * 1.015; // estimate from MC
static const double b13energyeff_e = 0.02; // BS

// time until end of run
static const double eor_eff = 1-(1-0.9709)*hightime/100e3;

static const double mich_eff = 0.9996;

static const double b12eff = mich_eff * eor_eff * b12energyeff;
static const double b12ferr_energy = b12energyeff_e/b12energyeff;

static const double b13eff = mich_eff * eor_eff * b13energyeff;
static const double b13ferr_energy = b13energyeff_e/b13energyeff;

// Using no distance cut ("nlate") so as not to have to figure out
// what the efficiency for that is for the high purity sample.
static const double neff = neff_dt_highpurity;

void fcn(int & npar, double * gin, double & like, double *par, int flag)
{
  const double n_b12 = par[0],
               n_b13 = par[1],
               n_li8 = par[2],
               n_n16 = par[3],
               acc   = par[4],
               b12t  = par[5],
               b13t  = par[6],
               li8t  = par[7],
               n16t  = par[8];

  like = 2*(n_b12*exp(-lowtime/b12t)
          + n_b13*exp(-lowtime/b13t)
          + n_li8*exp(-lowtime/li8t)
          + n_n16*exp(-lowtime/n16t)
          + totaltime*acc);

  printf(".");  fflush(stdout);

  for(unsigned int i = 0; i < events.size(); i++){
    const double mt = -events[i].t;
    const double f = n_b12/b12life*exp(mt/b12life) +
                     n_b13/b13life*exp(mt/b13life) +
                     n_li8/li8life*exp(mt/li8life) +
                     n_n16/n16life*exp(mt/n16life) +
                     acc;
    if(f > 0) like -= 2*log(f);
  }

  // pull terms
  like += pow((b12t - b12life)/b12life_err, 2);
  like += pow((b13t - b13life)/b13life_err, 2);
  like += pow((li8t - li8life)/li8life_err, 2);
  like += pow((n16t - n16life)/n16life_err, 2);
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
 
  const int npar = 9;
  mn = new TMinuit(npar);
  mn->SetPrintLevel(-1);
  mn->Command("SET STRATEGY 2");
  mn->SetFCN(fcn);
  int err;
  mn->mnparm(0, "n_b12", 1e4,  1e2, 0, 1e7, err);
  mn->mnparm(1, "n_b13",   1,  1e2, 0, 1e6, err);
  mn->mnparm(2, "n_li8", 1e3,  1e1, 0, 1e4, err);
  mn->mnparm(3, "n_n16",   1,  1e1, 0, 1e3, err);
  mn->mnparm(4, "acc",   10 ,    1, 0, 1e3, err);
  mn->mnparm(5, "b12t", b12life,  b12life_err, 0, 0, err);
  mn->mnparm(6, "b13t", b13life,  b13life_err, 0, 0, err);
  mn->mnparm(7, "li8t", li8life,  li8life_err, 0, 0, err);
  mn->mnparm(8, "n16t", n16life,  n16life_err, 0, 0, err);
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

  printf("MINOS B-12");
  mn->Command("MINOS 10000 1");
  puts(""); mn->Command("show min");

  printf("MINOS B-13");
  mn->Command("MINOS 10000 2");
  puts(""); mn->Command("show min");

  const double ferrorfit = (mn->fErp[0]-mn->fErn[0])/2/getpar(mn, 0);

  printf("\nStuff with b12 lifetime raw %f +- %f\n",
         getpar(mn, 0), ferrorfit * getpar(mn, 0));

  const double b12like_central = getpar(mn, 0)/b12eff/mum_count * 100;
  const double staterr = ferrorfit*b12like_central;
  const double muerr = mum_count_e/mum_count*b12like_central;
  const double b12err = b12ferr_energy*b12like_central;
  const double toterr = sqrt(pow(ferrorfit,2) +
                               pow(mum_count_e/mum_count,2) +
                               pow(b12ferr_energy, 2))*b12like_central;
  printf("\nStuff with b12 lifetime raw with overall eff "
       "corrected, %% per mu- stop\n"
       "%f +- %f (fit) +- %f (mu count) +- %f (B-12 eff), %f (total)\n",
       b12like_central, staterr, muerr, b12err, toterr);

  const double b12like_central_percap = b12like_central/(0.0762/0.998);
  const double staterr_percap = staterr/(0.0762/0.998);
  const double muerr_percap = muerr/(0.0762/0.998);
  const double b12err_percap = b12err/(0.0762/0.998);
  const double capfracerr_percap = b12like_central_percap
    * sqrt(pow(0.0092,2)+pow(0.003,2));
  const double toterr_percap = sqrt(pow(staterr_percap,2)+
                                    pow(muerr_percap,2)+
                                    pow(b12err_percap,2)+
                                    pow(capfracerr_percap,2));

  printf("\nOr %% per mu- capture\n"
         "%f +- %f (fit) +- %f (mu count) +- %f (B-12 eff),\n"
         "+- %f (capture fraction) %f (total)\n",
         b12like_central_percap, staterr_percap, muerr_percap,
         b12err_percap, capfracerr_percap, toterr_percap);

  const double b12like_central_rate = b12like_central/2.028e-6/1e5;
  const double staterr_rate = staterr/2.028e-6/1e5;
  const double muerr_rate = muerr/2.028e-6/1e5;
  const double b12err_rate = b12err/2.028e-6/1e5;
  const double lifetimeerr_rate = b12like_central_rate
    * 0.002e-6/2.028e-6;
  const double toterr_rate = sqrt(pow(staterr_rate,2)+
                                  pow(muerr_rate,2)+
                                  pow(b12err_rate,2)+
                                  pow(lifetimeerr_rate,2));

  printf("\nOr 10^3/s: %f +- %f (fit) +- %f (mu count) "
         "+- %f (B-12 eff),\n"
         "+- %f (lifetime) %f (total)\n",
         b12like_central_rate, staterr_rate, muerr_rate,
         b12err_rate, lifetimeerr_rate, toterr_rate);
}
