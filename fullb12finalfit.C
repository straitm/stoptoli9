#include <fstream>
#include "consts.h"
#include "TFile.h"
#include "TMinuit.h"
#include "TTree.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TROOT.h"
#include <stdio.h>
#include <vector>
#include <algorithm>
#include "deadtime.C" // <-- note inclusion of source

#define DISABLEN16

using std::vector;

const double b12life = 20.20/log(2.);
const double b13life = 17.33/log(2.); 

// This is the NNDC value. Could alternatively use 838.75+-0.32 from PRC
// 82, 027309
const double li8life = 839.9/log(2.);

const double n16life = 7130./log(2.);

const double b12life_err = 0.02/log(2.);
const double b13life_err = 0.17/log(2.);
const double li8life_err = 0.9/log(2.);
const double n16life_err = 20./log(2.);

struct ev{
  double t; // time
  int n; // number of neutrons
  double e; // nominal neutron efficiency
  ev(double t_, int n_, double e_){
    t = t_;
    n = n_;
    e = e_;
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
const double offset = 2.040e-3;

TH2D * hdisp = new TH2D("hdisp", "", 3, 0, 2, 250, 1-offset, 5001-offset);


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

// Constant 1% absolute error assumed on neutron efficiency.
// Not a very rigourous model, but it's something.
const double neff_err = 0.01;

// Measured probablity of getting one accidental neutron.  These
// are *detected* muons, so don't apply efficiency to them.
const double paccn = 1.1e-4;

/*
 * Prints the message once with the requested floating point precision
 * and in RED, then again with all digits in the default color, starting
 * with the first floating point number.
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

double clamp(double x, double lo, double hi)
{
  return x < lo? lo: x > hi? hi: x;
}

double zeron_f(
const double mt, const double acc,const double paccn, const double neff,
const double n_b12, const double n_b12n, const double b12t,
const double n_li8, const double n_li8n, const double li8t,
const double n_b13, const double b13t,
const double n_n16, const double n16t)
{
  return acc +
    // Be really careful. Can get zero neutrons from a real zero-neutron
    // event without an accidental or from a real one-neutron event with
    // an inefficiency and without an accidental
    ((1-paccn)*n_b12 + (1-neff)*(1-paccn)*n_b12n)/b12t*exp(mt/b12t) +
    ((1-paccn)*n_li8 + (1-neff)*(1-paccn)*n_li8n)/li8t*exp(mt/li8t) +
    // For B-13 and N-16, there are no real one-neutron events, so it
    // is easier. (Ok, N-16 might come from O-17, but it is only a
    // nuisance parameter here...)
    (1-paccn)*n_b13/b13t*exp(mt/b13t)
#ifndef DISABLEN16
    + (1-paccn)*n_n16/n16t*exp(mt/n16t)
#endif
    ;
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
               accn2 = par[8],
               b12t  = par[9],
               b13t  = par[10],
               li8t  = par[11],
               n16t  = par[12],
               neffdelta = par[13];

  const double norm =
      (n_b12 + n_b12n)*(exp(-lowtime/b12t) - exp(-hightime/b12t))
    + (n_li8 + n_li8n)*(exp(-lowtime/li8t) - exp(-hightime/li8t))
    + n_b13 *(exp(-lowtime/b13t) - exp(-hightime/b13t))
#ifndef DISABLEN16
    + n_n16 *(exp(-lowtime/n16t) - exp(-hightime/n16t))
#endif
    + totaltime*(acc+accn+accn2);

  like = norm;

  printf(".");  fflush(stdout);

  for(unsigned int i = 0; i < events.size(); i++){
    const double mt = -events[i].t;
    const double neff =
      clamp(events[i].e /* XXX neff_dt_highpurity  0.67 */ + neffdelta, 0.01, 0.95);

    double f = 0;
    switch(events[i].n){
      case 0:
        f = acc +
          // Be really careful. Can get zero neutrons from a real
          // zero-neutron event without an accidental or from a real
          // one-neutron event with an inefficiency and without an
          // accidental
          ((1-paccn)*n_b12 + (1-neff)*(1-paccn)*n_b12n)/b12t
            *exp(mt/b12t) +
          ((1-paccn)*n_li8 + (1-neff)*(1-paccn)*n_li8n)/li8t
            *exp(mt/li8t) +
           // For B-13 and N-16, there are no real one-neutron events,
           // so it is easier. (Ok, N-16 might come from O-17, but it 
           // is only a nuisance parameter here...)
           (1-paccn)*n_b13/b13t*exp(mt/b13t)
#ifndef DISABLEN16
           + (1-paccn)*n_n16/n16t*exp(mt/n16t)
#endif
          ;
        break;
      case 1:
        f = accn +
          // Can get one neutron from a real one-neutron event that
          // is efficient and doesn't have an accidental, or that is
          // inefficient and does have an accidental, or from a real
          // zero-neutron event that is inefficient.
          ((neff*(1-paccn)+(1-neff)*paccn)*n_b12n+paccn*n_b12)/b12t
            *exp(mt/b12t) +
          ((neff*(1-paccn)+(1-neff)*paccn)*n_li8n+paccn*n_li8)/li8t
            *exp(mt/li8t) +
           // Can only get B-13 with a neutron from an inefficiency,
           // assuming there's no O-16(mu,ppn)B-13 to speak of.
           paccn*n_b13/b13t*exp(mt/b13t)
#ifndef DISABLEN16
           + paccn*n_n16/n16t*exp(mt/n16t)
#endif
          ;
        break;
      case 2:
        f = accn2 +
          // Can get two neutrons from a real one-neutron event that
          // is efficient and has an accidental. Note that we *must*
          // handle this case for the above normalization component of
          // the likelihood to be correct.
          neff*paccn*n_b12n/b12t*exp(mt/b12t) +
          neff*paccn*n_li8n/li8t*exp(mt/li8t)
        ;
        break;
      default:
        printf("NOT REACHED with %d neutrons\n", events[i].n);
        break;
    }
    if(f > 0) like -= log(f);
  }

  like *= 2; // Convert to "chi2"

  // pull terms for lifetimes
  like += pow((b12t - b12life)/b12life_err, 2)
        + pow((b13t - b13life)/b13life_err, 2)
        + pow((li8t - li8life)/li8life_err, 2)
#ifndef DISABLEN16
        + pow((n16t - n16life)/n16life_err, 2)
#endif
          // and the neutron efficiency
        + pow(neffdelta/neff_err, 2);
  
  // pull term to impose unitarity bound on products of C-13.
  // Width is determined by the error on the number of captures
  //if(p_b12n + p_b13 + p_li8n > 1)
    //like += pow((p_b12n + p_b13 + p_li8n - 1)/errcapprob13*capprob13, 2);
}

double dispf(double * x, double * par)
{

}

static double getpar(int i)
{
  double answer, dum;
  mn->GetParameter(i, answer, dum);
  return answer;
}

double sani_minos(const double in)
{
  if(in == -54321) return 0;
  return in;
}

void results(const char * const iname, const int mni,
             const double ifrac, const double eff,
             const double ferr_energy, const double thiscapprob,
             const double thiscapprob_err, const double lifetime,
             const double lifetime_err, const double nuc_cap,
             const int prec1, const int prec2, const int prec3)
{
  // I did the fit in terms of probability per capture so that it was
  // straightforward to impose a unitarity bound, but this was only
  // a constant factor shift, now shift back to doing it in terms of
  // counts so I can step through the uncertainties.
  const double fit = nuc_cap*eff*getpar(mni);
  const double ferrorfitup = sani_minos(mn->fErp[mni])/getpar(mni);
  const double ferrorfitlo = sani_minos(mn->fErn[mni])/getpar(mni);

  printtwice("\n%s raw %f +%f %f\n", 0, iname,
    fit, ferrorfitup*fit, ferrorfitlo*fit);

  const double like_central = fit/eff/mum_count * 100/ifrac;
  const double staterrup = ferrorfitup*like_central;
  const double staterrlo = ferrorfitlo*like_central;
  const double muerr = mum_count_e/mum_count*like_central;
  const double err = ferr_energy*like_central;
  const double toterrup = sqrt(pow(ferrorfitup,2) +
                               pow(mum_count_e/mum_count,2) +
                               pow(ferr_energy, 2))*like_central;
  const double toterrlo = -sqrt(pow(ferrorfitlo,2) +
                               pow(mum_count_e/mum_count,2) +
                               pow(ferr_energy, 2))*like_central;
  printtwice("\n%s, eff corrected, percent per mu- stop\n"
    "%f +%f %f(fit) +-%f(mu count) +-%f(B-12 eff), +%f %f(total)\n",
    prec1, iname, like_central, staterrup, staterrlo, muerr, err,
    toterrup, toterrlo);

  const double like_central_percap = like_central/thiscapprob;
  const double staterr_percapup = staterrup/thiscapprob;
  const double staterr_percaplo = staterrlo/thiscapprob;
  const double muerr_percap = muerr/thiscapprob;
  const double err_percap = err/thiscapprob;
  const double capfracerr_percap =
    like_central_percap * thiscapprob_err/thiscapprob;
  const double toterr_percapup = sqrt(pow(staterr_percapup,2)+
                                    pow(muerr_percap,2)+
                                    pow(err_percap,2)+
                                    pow(capfracerr_percap,2));
  const double toterr_percaplo = -sqrt(pow(staterr_percaplo,2)+
                                    pow(muerr_percap,2)+
                                    pow(err_percap,2)+
                                    pow(capfracerr_percap,2));

  printtwice("\nOr percent per mu- capture\n"
         "%f +%f %f(fit) +-%f(mu count) +-%f(eff) +-%f(cap frac), "
         "+%f %f(total)\n", prec2, 
         like_central_percap, staterr_percapup, staterr_percaplo,
         muerr_percap, err_percap, capfracerr_percap,
         toterr_percapup, toterr_percaplo);

  const double like_central_rate = like_central/lifetime/100;
  const double staterr_rateup = staterrup/lifetime/100;
  const double staterr_ratelo = staterrlo/lifetime/100;
  const double muerr_rate = muerr/lifetime/100;
  const double err_rate = err/lifetime/100;
  const double lifetimeerr_rate = like_central_rate
    * lifetime_err/lifetime;
  const double toterr_rateup = sqrt(pow(staterr_rateup,2)+
                                  pow(muerr_rate,2)+
                                  pow(err_rate,2)+
                                  pow(lifetimeerr_rate,2));
  const double toterr_ratelo = -sqrt(pow(staterr_ratelo,2)+
                                  pow(muerr_rate,2)+
                                  pow(err_rate,2)+
                                  pow(lifetimeerr_rate,2));

  printtwice("\nOr 10^3/s: %f +%f %f(fit) +-%f(mu count) +-%f(eff),\n"
         "+-%f(lifetime) +%f %f(total)\n", prec3,
         like_central_rate, staterr_rateup, staterr_ratelo,
         muerr_rate, err_rate, lifetimeerr_rate, toterr_rateup,
         toterr_ratelo);
  puts("");
}

void printc12b12results()
{
  results("C-12 -> B-12", 0, 1-f13, b12eff, b12ferr_energy, capprob12,
    errcapprob12, lifetime_c12, lifetime_c12_err, c12nuc_cap, 3, 2, 2);
}

void printc13b12nresults()
{
  results("C-13 -> B-12+n", 1, f13, b12eff, b12ferr_energy, capprob13,
    errcapprob13, lifetime_c13, lifetime_c13_err, c13nuc_cap, 2, 1, 1);
}

void printc13b13results()
{
  results("C-13 -> B-13", 2, f13, b13eff, b13ferr_energy, capprob13,
    errcapprob13, lifetime_c13, lifetime_c13_err, c13nuc_cap, 2, 1, 1);
}


void fullb12finalfit(const char * const cut =
"mx**2+my**2 < 1050**2 && mz > -1175 && "
"abs(fez + 62*ivdedx/2 - 8847.2) < 1000 && chi2 < 2 && "
"timeleft > %f && miche < 12 && !earlymich && "
"e > 4 && e < 15 && dt < %f && laten <= 2")
{
  printtwice("B-12 selection efficiency is %.1f%%\n", 2, b12eff*100);

  TFile *_file0 = TFile::Open(rootfile3up);
  TTree * t = (TTree *)_file0->Get("t");
 
  const int npar = 14;
  mn = new TMinuit(npar);
  mn->SetPrintLevel(-1);
  mn->fGraphicsMode = false;
  //mn->Command("SET STRATEGY 2");
  mn->SetFCN(fcn);
  int err;
  mn->mnparm(0, "p_b12", 0.2,  0.01, 0, 1, err);
  mn->mnparm(1, "p_b12n",0.5,  0.01, 0, 1, err);
  mn->mnparm(2, "p_b13", 0.2,  0.01, 0, 1, err);
  mn->mnparm(3, "p_li8", 0.01,  0.01, 0, 1, err);
  mn->mnparm(4, "p_li8n",0.05,  0.01, 0, 1, err);
  mn->mnparm(5, "n_n16",   1,  1, 0, 1e3, err);
  mn->mnparm(6, "acc",   1,    0.1, 0, 100, err);
  mn->mnparm(7, "accn",  0.1,  0.01, 0, 10, err);
  mn->mnparm(8, "accn2", 0.01, 0.001, 0, 1, err);
  mn->mnparm(9, "b12t", b12life,  b12life_err, 0, 0, err);
  mn->mnparm(10, "b13t", b13life,  b13life_err, 0, 0, err);
  mn->mnparm(11, "li8t", li8life,  li8life_err, 0, 0, err);
  mn->mnparm(12, "n16t", n16life,  n16life_err, 0, 0, err);
  mn->mnparm(13, "neffdelta", 0,  neff_err, 0, 0, err);

#ifdef DISABLEN16
  mn->Command("SET PAR 6 0");
  mn->Command("FIX 6");
  mn->Command("FIX 13");
#endif

  // XXX fix Li-8 lifetime for faster fitting, since it really shouldn't
  // matter
  mn->Command("FIX 12");

  mn->Command("FIX 14"); // neutron efficiency

  mn->Command("FIX 10");  // more lifetimes
  mn->Command("FIX 11");
  
  printf("Making cuts...\n");
  TFile * tmpfile = new TFile("/tmp/b12tmp.root", "recreate");
  selt = t->CopyTree(Form(cut, hightime+offset, hightime+offset));
  selt->Write();
  events.clear();

  float dt, mx, my, mz, fq;
  int nn;
  selt->SetBranchAddress("dt", &dt);
  selt->SetBranchAddress("laten", &nn);
  selt->SetBranchAddress("mx", &mx);
  selt->SetBranchAddress("my", &my);
  selt->SetBranchAddress("mz", &mz);
  selt->SetBranchAddress("fq", &fq);

  printf("Filling in data array...");
  for(int i = 0; i < selt->GetEntries(); i++){
    selt->GetEntry(i);
    events.push_back(ev(dt-offset, nn, eff(fq, mx, my, mz)));
    hdisp->Fill(nn, dt-offset);
    if(i%10000 == 9999){ printf("."); fflush(stdout); }
  }

  {
    double emean = 0;
    for(unsigned int i = 0; i < events.size(); i++)
      emean += events[i].e/events.size();
    printf("Mean neutron efficiency: %f\n", emean);
  }

  printf("\nMIGRAD");
  mn->Command("MIGRAD");
  puts(""); mn->Command("show par");

/*
  mn->SetPrintLevel(2);
  mn->Command("MINOS 2000 1 2 3");
  mn->SetPrintLevel(-1);
  puts("");
  mn->Command("show min");

  printc12b12results();
  printc13b12nresults();
  printc13b13results();
*/
  mn->SetPrintLevel(1);
}
