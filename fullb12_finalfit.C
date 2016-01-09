#include <fstream>
#include <math.h>
#include "sub_muon_eff.out.h"
#include "totallivetime_finalfit.out.h"
#include "li8cutefficiency_finalfit.out.h"
#include "b12cutefficiency_finalfit.out.h"

/*** BEGIN For near detector pull terms ***/
#include "li9_finalfit_1.out.h"
#include "li9_finalfit_-1.out.h"
#include "li8_finalfit.out.h"

#include "noncarbondenominators_finalfit.out.h"

// XXX awkward since carbondenominators has mum_count in it.
const double n_c12cap     = 360.9878604679215073;
//#include "carbondenominators_finalfit.out.h"
/*** END For near detector pull terms ***/

#include "TFile.h"
#include "TMinuit.h"
#include "TTree.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TError.h"
#include "TROOT.h"
#include <stdio.h>
#include <vector>
using std::pair;
using std::string;
using std::cin;
#include <algorithm>
#include "neff.C" // <-- note inclusion of source

// Should only be changed by near_fullb12_finalfit()
static bool near = false;
static string NEUTRONDEF = "laten";

// Reset at the top of fullb12_finalfit
bool apply_dr_eff = true; 

using std::vector;

bool unitarity = true;

struct ev{
  double mt; // negative time
  double e; // nominal neutron efficiency
  int n; // number of neutrons
  ev(double t /* positive time goes in */, int n_, double e_){
    mt = -t; // negative time comes out
    n = n_;
    e = e_;
  }
};

const int npar = 18;
TMinuit * mn = NULL;
TTree * selt = NULL;

vector<ev> events;

// Mean neutron efficiency for events that actually have
// zero, one or two neutrons;
double meaneff[3] = {0};
int eventc[3] = {0}; // number of events with each n count

/*
 * Prints the message once with the requested floating point precision
 * and in RED, then again with all digits in the default color, starting
 * with the first floating point number.
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
  vprintf(pmsg, ap);

  va_start(ap, digits);
  vprintf(bmsg, ap);
}

/**********************************************************************/

#include "mucount_finalfit.C"

// set below
static string countcut = "asergd!";
static string CUT_PART_OK_FOR_FINDING_PACCN  = "skrree!";

// The number of mu- stopping, regardless of what they atomicly or
// nuclearly capture on, set below
static double mum_count   = 0;
static double mum_count_e = 0;

/**********************************************************************/

const double lifetime_c = lifetime_c12*(1-f13_HP)+lifetime_c13*f13_HP;
const double lifetime_c_err = sqrt(pow(lifetime_c12_err*(1-f13_HP),2)
                                  +pow(lifetime_c13_err*   f13_HP ,2));

const double capprob12 = 1-lifetime_c12/mulife;
const double errcapprob12 = (1-(lifetime_c12+lifetime_c12_err)/mulife)/2
                           -(1-(lifetime_c12-lifetime_c12_err)/mulife)/2;

const double capprob13 = 1-lifetime_c13/mulife;
const double errcapprob13 = (1-(lifetime_c13+lifetime_c13_err)/mulife)/2
                           -(1-(lifetime_c13-lifetime_c13_err)/mulife)/2;

const double capprob = c_atomic_capture_prob *
                      (capprob12*(1-f13_HP) + capprob13*f13_HP);

const double err_capprob = sqrt(pow(errcapprob12,2)*(1-f13_HP)
                              + pow(errcapprob13,2)*f13_HP +
  pow(c_atomic_capture_prob_err/c_atomic_capture_prob * capprob, 2));

// Set below
double c12nuc_cap = 0;
double c13nuc_cap = 0;


// Will subtract mean muon lifetime, 2028ns, and mean transit time for
// light from B-12, 12ns. Doesn't make a real difference.
const double offset = 2.040e-3;

const double lowtime = 1.0 - offset;
const double hightime = 100e3;
const double totaltime = hightime - lowtime;

const double b13energyeff = b12energyeff * 1.014; // estimate from my MC
const double b13energyeff_e = 0.02; // BS

const double li9eff_energy = 0.8069+0.05; // DOGS, with ad hoc
                                          // correction for the beta
                                          // branches being the
                                          // relevant ones rather
                                          // than the betan branches
const double li9eff_energy_e = 0.05; // made up!

// time until end of run
const double eor_eff = (livetime_s - num_runs*(hightime+offset)/1e3)/livetime_s;

/* const */ double b12eff = mich_eff * light_noise_eff * eor_eff * sub_muon_eff05 * b12energyeff;
const double b12ferr_energy = b12energyeff_e/b12energyeff;

/* const */ double b13eff = mich_eff * light_noise_eff * eor_eff * sub_muon_eff05 * b13energyeff;
const double b13ferr_energy = b13energyeff_e/b13energyeff;

/* const */ double li8eff = mich_eff * light_noise_eff * eor_eff * sub_muon_eff05 * li8energyeff4MeV;

/* const */ double li9eff = mich_eff * light_noise_eff * eor_eff * sub_muon_eff05 * li9eff_energy;

// Measured probablity of getting one accidental neutron.
static double nom_paccn = 0; // replaced by find_paccn()
static double paccn_e = 0; // ditto

bool isibd(const int in_run, const int in_trig)
{
  if(near) return false; // XXX until/unless we have an IBD list

  static bool firsttime = true;
  static vector< pair<int,int> > ibds;

  if(firsttime){
    TTree ibd;
    firsttime = false;
    ibd.ReadFile("/cp/s4/strait/li9ntuples/Hprompts20141119", "run:trig");
    ibd.ReadFile("/cp/s4/strait/li9ntuples/Gdprompts20140925");
    float run, trig; // Yes, float
    ibd.SetBranchAddress("run", &run);
    ibd.SetBranchAddress("trig", &trig);
    for(int i = 0; i < ibd.GetEntries(); i++){
      ibd.GetEntry(i);
      ibds.push_back(pair<int,int>(int(run), int(trig)));

      // Get the delayed event too.  ROUGH, since there could be
      // an intervening event
      ibds.push_back(pair<int,int>(int(run), int(trig+1)));
    }
  }

  // Really stupid linear search, but not quite as stupid
  // as asking the TTree how many matches it has.
  for(unsigned int i = 0; i < ibds.size(); i++)
    if(ibds[i] == pair<int,int>(in_run, in_trig))
      return true;
  return false;
}

static inline double clamp(double x, double lo, double hi)
{
  return x < lo? lo: x > hi? hi: x;
}

static inline double mexp(const double arg)
{
  return arg < -746? 0: exp(arg);
}

// Return mult*exp(arg) if adding this to prev could possibly make a difference.
// This is only OK for negative arg! (And it is always negative for us.)
//
// Get an overall improvement of 20% in speed from this!
static inline double expif(const double arg, const double mult, const double prev)
{
  // An upper bound for how big exp(arg) is is
  // 2^((int(arg)+1)*log_2(e)), e.g. for arg=-137.1, we know that
  // exp(arg) < 2**(-136/M_LN_2) So the 2**-136 place is an upper bound
  // on the bit this could affect. If the lowest bit in prev is higher
  // than this, exp(arg) is irrelevant. THe lowest bit is the exponent
  // of prev - 53, where 53 is the number of significant bits in the
  // mantissa.
  
  int multexponent;
  frexp(mult, &multexponent);

  const int upbound = (int(arg) + 1 + multexponent+1)/M_LN2;
  int prevexponent;
  frexp(prev, &prevexponent);
  const int lowdigit = prevexponent - 53;
  if(upbound < lowdigit) return 0;
  return mult*exp(arg); 
}

static inline double justexp(const double arg, const double mult,
__attribute__((unused)) double prev)
{
  return mult*exp(arg);
}

static inline double zeron_f(
const double mt,const double acc,const double unpaccn,const double neff,
const double unpaccn_x_n_b12_x_b12r, const double n_b12n_x_b12r,const double b12r,
const double unpaccn_x_n_li8_x_li8r, const double n_li8n_x_li8r,const double li8r,
const double unpaccn_x_n_li9_x_li9r, const double n_li9n_x_li9r,const double li9r,
const double unpaccn_x_n_b13_x_b13r, const double b13r,
const double unpaccn_x_n_n16_x_n16r, const double n16r)
{
  // Be really careful. Can get zero neutrons from a real zero-neutron
  // event without an accidental or from a real one-neutron event with
  // an inefficiency and without an accidental.
  const double uneffunpaccn = (1-neff)*unpaccn; // function of pars AND event
      
  // For B-13 and N-16, there are no real one-neutron events, so it is
  // easier. (Ok, N-16 might come from O-17, but it is only a nuisance
  // parameter here...)
      
  double answer = acc;
  answer += justexp(mt*n16r, unpaccn_x_n_n16_x_n16r, answer);

  answer += expif(mt*li8r, unpaccn_x_n_li8_x_li8r
                     + uneffunpaccn*n_li8n_x_li8r, answer);
  answer += expif(mt*li9r, unpaccn_x_n_li9_x_li9r
                     + uneffunpaccn*n_li9n_x_li9r, answer);
  answer += expif(mt*b12r, unpaccn_x_n_b12_x_b12r
                     + uneffunpaccn*n_b12n_x_b12r, answer);
  answer += expif(mt*b13r, unpaccn_x_n_b13_x_b13r, answer);
      
  return answer;
}

static inline double onen_f(
const double mt,const double accn,const double paccn, const double neff,
const double n_b12, const double n_b12n, const double b12r,
const double n_li8, const double n_li8n, const double li8r,
const double n_li9, const double n_li9n, const double li9r,
const double n_b13, const double b13r,
const double n_n16, const double n16r)
{
  return accn +
  // Can get one neutron from a real one-neutron event that is efficient
  // and doesn't have an accidental, or that is inefficient and does
  // have an accidental, or from a real zero-neutron event that is
  // inefficient.
  ((neff*(1-paccn)+(1-neff)*paccn)*n_b12n+paccn*n_b12)*b12r*exp(mt*b12r)+
  // Can only get B-13 with a neutron from an inefficiency, assuming
  // there's no O-16(mu,ppn)B-13 to speak of.
  paccn*n_b13*b13r*exp(mt*b13r)
 +((neff*(1-paccn)+(1-neff)*paccn)*n_li8n+paccn*n_li8)*li8r*exp(mt*li8r)
 +((neff*(1-paccn)+(1-neff)*paccn)*n_li9n+paccn*n_li9)*li9r*exp(mt*li9r)
  + paccn*n_n16*n16r*exp(mt*n16r)
  ;
}

static inline double twon_f(
const double mt,const double accn2,const double neff,const double paccn,
const double n_b12n, const double b12r,
const double n_li8n, const double li8r,
const double n_li9n, const double li9r)
{
  return accn2 +
  // Can get two neutrons from a real one-neutron event that is
  // efficient and has an accidental. Note that we *must* handle this
  // case for the above normalization component of the likelihood to be
  // correct.
  neff*paccn*(n_b12n*b12r*exp(mt*b12r)
            + n_li8n*li8r*exp(mt*li8r)
            + n_li9n*li9r*exp(mt*li9r)
  );
}

#define DECODEPARS \
  const double p_b12 = par[0], \
               p_b12n= par[1], \
               p_b13 = par[2], \
               p_li8 = par[3], \
               p_li8n= par[4], \
               p_li9 = par[5], \
               p_li9n= par[6]; \
  const double n_b12 = p_b12 *c12nuc_cap*b12eff, \
               n_b12n= p_b12n*c13nuc_cap*b12eff, \
               n_b13 = p_b13 *c13nuc_cap*b13eff, \
               n_li8 = p_li8 *c12nuc_cap*li8eff, \
               n_li8n= p_li8n*c13nuc_cap*li8eff, \
               n_li9 = p_li9 *c13nuc_cap*li9eff, \
               n_li9n= p_li9n*c12nuc_cap*li9eff, \
               n_n16 = par[7], \
               acc   = par[8], \
               accn  = par[9], \
               accn2 = par[10], \
               b12t  = par[11]*b12life_err + b12life, \
               b13t  = par[12]*b13life_err + b13life, \
               li8t  = par[13]*li8life_err + li8life, \
               li9t  = par[14]*li9life_err + li9life, \
               n16t  = par[15]*n16life_err + n16life, \
               b12r  = 1/b12t, \
               b13r  = 1/b13t, \
               li8r  = 1/li8t, \
               li9r  = 1/li9t, \
               n16r  = 1/n16t, \
               neffdelta = par[16], \
               paccn = par[17];

double priorerf(const double sig)
{
  return 0.5*(1 + erf(-sig/sqrt(2)));
}

double unit_penalty(const double x)
{
  static double unitwidth =
    sqrt(pow(errcapprob13/capprob13,2)+pow(mum_count_e/mum_count,2));

  const double sig = (x-1)/unitwidth;

  // When erf gets very close to 0, bad things happen, so switch to
  // an approximation past 5 sigma out. There may well be a better
  // way to handle this.
  static const double sigcut = 5;
  static const double corr = -2*log(priorerf(sigcut)) - pow(sigcut,2);
  const double mlogprior =
    sig < sigcut? -2*log(priorerf(sig)) : corr + pow(sig, 2);

  return mlogprior;
}

// Sums the input array, losing the minimum precision via the Kahan
// Summation Algorithm.
// See https://en.wikipedia.org/wiki/Kahan_summation_algorithm
//
// This may help MIGRAD converge by providing a smoother function. Note,
// however, that it requires storing all the intermediate numbers in
// memory, which in some cases means reducing from 80-bit floats to
// 64-bit ones. But at least it will be consistent, and I think that
// MINUIT figures out what the effective precision is.
//
// In theory, an optimizing compiler with options like GCC's -ffast-math
// set could analyze this function, see that c is algebraically always
// zero, and eliminate it, thereby destroying the algorithm. In
// practice, we use GCC and it doesn't happen.
static inline double kahan_sum(const double * const numbers,
                              const unsigned int n)
{
  double sum = 0, c = 0;
  for(unsigned int i = 0; i < n; i++){
    const double y = numbers[i] - c;
    const double t = sum + y;
    c = (t - sum) - y;
    sum = t;
  }
  return sum;
}

void fcn(__attribute__((unused)) int & X,
         __attribute__((unused)) double * gin,
                                 double & like,
                                 double *par,
         __attribute__((unused)) int flag)
{
  DECODEPARS;

  const double norm =
      (n_b12 + n_b12n)*(exp(-lowtime*b12r) - exp(-hightime*b12r))
    + (n_li8 + n_li8n)*(exp(-lowtime*li8r) - exp(-hightime*li8r))
    + (n_li9 + n_li9n)*(exp(-lowtime*li9r) - exp(-hightime*li9r))
    +      n_b13      *(exp(-lowtime*b13r) - exp(-hightime*b13r))
    +      n_n16      *(exp(-lowtime*n16r) - exp(-hightime*n16r))
    + totaltime*(acc+accn+accn2);

  like = norm;

  // These products are invariant over the set of events for a given 
  // set of parameters, so do them once out of the loop
  const double unpaccn = 1-paccn;
  const double unpaccn_x_n_b13_x_b13r = unpaccn*n_b13*b13r;
  const double unpaccn_x_n_n16_x_n16r = unpaccn*n_n16*n16r;
  const double unpaccn_x_n_b12_x_b12r = unpaccn*n_b12*b12r;
  const double unpaccn_x_n_li8_x_li8r = unpaccn*n_li8*li8r;
  const double unpaccn_x_n_li9_x_li9r = unpaccn*n_li9*li9r;
  const double n_b12n_x_b12r = n_b12n*b12r;
  const double n_li8n_x_li8r = n_li8n*li8r;
  const double n_li9n_x_li9r = n_li9n*li9r;

  vector<double> per_event_like(events.size());
  for(unsigned int i = 0; i < events.size(); i++){
    const double neff = clamp(events[i].e + neffdelta, 0.00, 1.00);

    double f = 0;
    switch(events[i].n){
      case 0:
        f = zeron_f(events[i].mt, acc, unpaccn, neff,
                    unpaccn_x_n_b12_x_b12r, n_b12n_x_b12r, b12r,
                    unpaccn_x_n_li8_x_li8r, n_li8n_x_li8r, li8r,
                    unpaccn_x_n_li9_x_li9r, n_li9n_x_li9r, li9r,
                    unpaccn_x_n_b13_x_b13r, b13r,
                    unpaccn_x_n_n16_x_n16r, n16r);
        break;
      case 1:
        f = onen_f(events[i].mt, accn, paccn, neff, n_b12, n_b12n, b12r,
                   n_li8, n_li8n, li8r, n_li9, n_li9n, li9r,
                   n_b13, b13r, n_n16, n16r);
        break;
      case 2:
        f = twon_f(events[i].mt, accn2, neff, paccn, n_b12n, b12r,
                   n_li8n, li8r, n_li9n, li9r);
        break;
      default:
        // Has turned out to be a useful check!
        printf("NOT REACHED with %d neutrons\n", events[i].n);
        break;
    }
    per_event_like[i] = f > 0? -log(f): 0;
  }

  like += kahan_sum(&(per_event_like[0]), events.size());

  like *= 2; // Convert to "chi2"

  // pull terms for lifetimes
  like += pow((b12t - b12life)/b12life_err, 2)
        + pow((b13t - b13life)/b13life_err, 2)
        + pow((li8t - li8life)/li8life_err, 2)
        + pow((li9t - li9life)/li9life_err, 2)
        + pow((n16t - n16life)/n16life_err, 2)
          // and the neutron efficiency
        + pow(neffdelta/sqrt(pow(f_neff_dt_error, 2) +
                apply_dr_eff?pow(f_neff_dr_800_avg_error, 2):0
                            ), 2)

          // and the accidental neutron probability
        + pow((paccn - nom_paccn)/paccn_e, 2);

  if(near){
    // Use FD constraints on Li-8 and Li-9 production
    like +=
    + pow((p_li8 -probEightLiFromTwelveC)/probEightLiFromTwelveC_statup, 2)
    + pow((p_li8n-probEightLiFromThirteenC)/probEightLiFromThirteenC_statup, 2)
    + pow((p_li9 -probNineLiFromTwelveC)/probNineLiFromTwelveC_uperr, 2)
    + pow((p_li9n-primaryresult/f13)/(primaryresult_uperr/f13), 2);

    // Use Measday's data for constraint on N-16 production
    const double o16nuc_cap = c12nuc_cap*n_o16cap_beta/n_c12cap;
    const double n16expected = o16nuc_cap * n16prob_measday;
    const double n16err = o16nuc_cap * n16prob_measday_err;
    like += pow((n_n16 - n16expected)/n16err, 2);
  }


  // pull terms for Li-9 from the betan analysis. Assume zero production
  // with a neutron so as not to double count. Since IBD candidates are
  // cut, this is only the non-betan rate. The multiper is to account
  // for the uncertainty in the Li-9 energy cut.
  static const double energymultiplier = 1 + li9eff_energy_e/li9eff_energy;

  //                      XXX no betan list for ND yet
  like += pow( p_li9n/(1-(!near)*0.508)/0.44e-4/energymultiplier, 2)
        + pow((p_li9 /(1-(!near)*0.508) - 2.4e-4)/0.9e-4/energymultiplier, 2);


  // Pull term to impose unitarity bound on products of C-13. Width is
  // determined by the error on the number of captures The concept here
  // is that there is a prior for the total number of captures with the
  // central value normalized to one, normally distributed. Therefore,
  // the prior for the sum of some isotopes is flat until it gets near
  // the total, and then drops off like the error function. Far away
  // from the central value, this looks a lot like adding a quadratic
  // if x>1 (I believe it approaches this asymptotically), but it is
  // gentler near 1.
  //
  // Do *not* inclue li-8 and li-9 in the unitarity penalty, because the
  // fit really likes adding lots of them in and being really sure about
  // it, which forces B-13 to a very low value. I don't think this is
  // physical. I suspect the presence of some other isotope that doesn't
  // come from C-13.  With this hypothesis, it makes sense to allow
  // them to float freely in the fit, since they are representing some
  // unknown background.  Could they be spallation reactions farther
  // up the muon track, here admitted since I don't use a distance cut?
  if(unitarity) like += unit_penalty(p_b12n +p_b13 /* +p_li8n +p_li9*/);

  if(mn->fCfrom != "CONtour   "){
    static int dotcount = 0;
    if(dotcount++%10 == 9) printf("%d %16f\n", dotcount, like);
    fflush(stdout);
  }
}

double getpar(int i)
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

  printtwice("\n%s%s, eff corrected, percent per C mu- stop "
    "%f +%f %f(fit) +-%f(mu count) +-%f(B-12 eff), "
    "+%f %f(total),  +%f -%f (non-fit)\n",
    prec1,
    !strcmp(iname, "C-12 -> B-12")?
      "TECHNOTE 4.3: ":
    !strcmp(iname, "C-13 -> B-12+n")?
      "TECHNOTE 4.3: ":
    !strcmp(iname, "C-13 -> B-13")?"":"?",
    iname, like_central, staterrup, staterrlo, muerr, err,
    toterrup, toterrlo,
    sqrt(pow(toterrup,2) - pow(staterrup,2)),
    sqrt(pow(toterrlo,2) - pow(staterrlo,2))
  );

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

  printtwice("\n%sPercent per nuclear mu- capture on this isotope "
         "%f +%f %f(fit) +-%f(mu count) +-%f(eff) +-%f(cap frac), "
         "+%f %f(total),  +%f -%f (non-fit)\n",
         prec2, 
         !strcmp(iname, "C-12 -> B-12")?
           "TECHNOTE 4.3 and results.tex probTwelveBFromTwelveC: ":
         !strcmp(iname, "C-13 -> B-12+n")?
           "TECHNOTE 4.3 and results.tex probTwelveBFromThirteenC: ":
         !strcmp(iname, "C-13 -> B-13")?"":"?",
         like_central_percap, staterr_percapup, staterr_percaplo,
         muerr_percap, err_percap, capfracerr_percap,
         toterr_percapup, toterr_percaplo,
         sqrt(pow(toterr_percapup,2) - pow(staterr_percapup,2)),
         sqrt(pow(toterr_percaplo,2) - pow(staterr_percaplo,2))
         );

  if(!strcmp(iname, "C-13 -> B-12+n")){
    printf("const double probTwelveBFromThirteenC = %.16f;\n",
           like_central_percap/100.);
    printf("const double probTwelveBFromThirteenC_statup = %.16f;\n",
           staterr_percapup/100.);
    printf("const double probTwelveBFromThirteenC_statlo = %.16f;\n",
           staterr_percaplo/100.);
  }
  else if(!strcmp(iname, "C-12 -> B-12")){
    printf("const double probTwelveBFromTwelveC = %.16f;\n",
           like_central_percap/100.);
    printf("const double probTwelveBFromTwelveC_statup = %.16f;\n",
           staterr_percapup/100.);
    printf("const double probTwelveBFromTwelveC_statlo = %.16f;\n",
           staterr_percaplo/100.);
  }

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

  printtwice("\n%s10^3/s: %f +%f %f(fit) +-%f(mu count) +-%f(eff), "
         "+-%f(lifetime) +%f %f(total)  +%f -%f\n",
         prec3,
         !strcmp("C-12 -> B-12", iname)?   "TECHNOTE 4.3: ":
         !strcmp("C-13 -> B-12+n", iname)? "TECHNOTE 4.3: ": "",
         like_central_rate, staterr_rateup, staterr_ratelo,
         muerr_rate, err_rate, lifetimeerr_rate, toterr_rateup,
         toterr_ratelo,
         sqrt(pow(toterr_rateup, 2) - pow(staterr_rateup,2)),
         sqrt(pow(toterr_ratelo, 2) - pow(staterr_ratelo,2))
         );
  if(!strcmp("C-12 -> B-12", iname)){
    printf("const double b12totalrate = %.16f;\n", like_central_rate);
    printf("const double b12totalrate_statup = %.16f;\n", staterr_rateup);
    printf("const double b12totalrate_statlo = %.16f;\n", staterr_ratelo);
    printf("const double b12totalrate_syst = %.16f;\n",
      sqrt(pow(muerr_rate,2)+pow(err_rate,2)+pow(lifetimeerr_rate,2)));
  }

  puts("");
}

void printc12b12results()
{
  results("C-12 -> B-12", 0, 1-f13_HP, b12eff, b12ferr_energy,
    capprob12*c_atomic_capture_prob, errcapprob12, lifetime_c12,
    lifetime_c12_err, c12nuc_cap, 3, 2, 2);
}

void printc13b12nresults()
{
  results("C-13 -> B-12+n", 1, f13_HP, b12eff, b12ferr_energy,
    capprob13*c_atomic_capture_prob, errcapprob13, lifetime_c13,
    lifetime_c13_err, c13nuc_cap, 2, 1, 1);
}

void printc13b13results()
{
  results("C-13 -> B-13", 2, f13_HP, b13eff, b13ferr_energy,
    capprob13*c_atomic_capture_prob, errcapprob13, lifetime_c13,
    lifetime_c13_err, c13nuc_cap, 2, 1, 1);
}

void mncommand()
{
  string command;
  while(true){
    printf("MINUIT> ");
    if(!getline(cin, command)) break;
    if(command == "exit") break;
    mn->Command(command.c_str());
  }
}

double b13limit()
{
  double sump = 0;

  const double increment = 0.1;
  const int N = 11;
  double ps[N];

  const double bestchi2 = mn->fAmin;

  mn->Command("fix 3");
  mn->Command("set strategy 0");
  for(int i = 0; i < N; i++){
    const double prob = i*increment;
    mn->Command(Form("set par 3 %f", prob));
    mn->Command("Migrad 2000 1");

    if(prob == 0 || fabs(prob-1) < 1e-5){
      printf("FIGURE const double pars_%sb13[%d] = { ",
             prob == 0?"no":"all", npar);
      for(int j = 0; j < npar; j++) printf("%.9f, ", getpar(j));
      printf("};\n");
    }

    const double p = exp(bestchi2-mn->fAmin);
    printf("\n%8.6f %8.3g %8.3g ", prob, mn->fAmin-bestchi2, p);
    for(int j = 0; j < p*10 - 1; j++) printf("#");
    if     (p*10 - int(p*10) > 0.67) printf("+");
    else if(p*10 - int(p*10) > 0.33) printf("|");
    printf("\n");
    sump += p;
    ps[i] = p;
  }

  printf("Norm: %f\n", sump);

  double sump2 = 0;
  double answer = 0;
  for(int i = 0; i < N; i++){
    const double prob = i*increment;
    sump2 += ps[i]/sump;
    if(sump2 > 0.9){
       answer = prob; // had been subtracting increment/2, but
                      // that makes the answer look funny.
                      // Just be conservative and take the
                      // next sample point after the crossing.
       printf("TECHNOTE 4.3: C-13 -> B-13 per nuclear capture "
              "90%% limit = %.0f%%\n", answer*100);
       break;
    }
  }

  mn->Command("release 3");
  mn->Command("Migrad");

  return answer;
}


// Measured probablity of getting one accidental neutron. These are 
// *detected* neutrons, so don't apply efficiency to them.
void find_paccn(TTree * t)
{
  // Assume that the only source of events 0-5.5mus after the muon
  // stop between 45 and 80 MeV is muon decay. Then all neutrons
  // following this are accidental, and we can directly count to find
  // the accidental probability for this sample.
  //
  // Complication: There is radiative capture to B-11 + n.  The photon
  // can, in principle, go to up ~88MeV, but I think it is very strongly
  // supressed at higher energies, so I am hoping that above 45MeV it is
  // negligible.
  //
  // Compliction: I am ignoring multiple neutrons.  These are rare
  // compared to single neutrons, but clearly come correlated to each
  // other (i.e. the probability of two is not P(one)^2).

  const double zero = t->GetEntries((CUT_PART_OK_FOR_FINDING_PACCN + " && "
    "miche > 45 && miche < 80 && " + NEUTRONDEF + " == 0 && ndecay == 0").c_str());
  const double one = t->GetEntries((CUT_PART_OK_FOR_FINDING_PACCN + " && "
    "miche > 45 && miche < 80 && " + NEUTRONDEF + " == 1 && ndecay == 0").c_str());

  // To cover the several complications above (hopefully), assume
  // that the rate of accidentals is overestimated by this amount,
  // and take a gaussian systematic equal to this amount (this
  // intentionally covers, to some extent, the possibility that we are
  // *underestimating* the rate, too).
  const double SYST = 0.05;

  nom_paccn = one/(zero + one) * (1-SYST);
  const double paccn_e_stat = sqrt(one)/(zero + one);
  const double paccn_e_syst = nom_paccn * SYST;
  paccn_e = sqrt(pow(paccn_e_stat, 2) + pow(paccn_e_syst, 2));
  
  printf("Accidental neutron prob: %.6g +- %.6g\n", nom_paccn, paccn_e);
}

static bool compare_events(const ev & a, const ev & b)
{
  // put the late time ones first on the theory that starting the
  // likelihood sum with smaller numbers makes it slightly more
  // accurate.  Might not be true, but shouldn't hurt anyway.
  if(a.n != b.n) return a.mt > b.mt;

  // Primary sort is by neutron number to aid branch prediction in fcn()
  return a.n < b.n;
}

void fullb12_finalfit()
{
  gErrorIgnoreLevel = kError;

  // Known by careful study to be good for FD (see tech note).  For ND, 
  // with a very rough study, mz looks pure down to -1300 or so and
  // sqrt(mx**2+my**2) out to 1400 given mz > -1300. So reusing the 
  // FD cuts should be quite conservative.
  string STD_POS_CUT = "mx**2+my**2 < 1050**2 && mz > -1175";

  // XXX dr efficiency
  if(near){
    const double dreff = 1; // /* for 600 vs. 1000 : */ 0.65
                         // /* * * for 1000: */ (0.1603 + 0.003)/0.1735;
    b12eff *= dreff;
    b13eff *= dreff;
    li8eff *= dreff;
    li9eff *= dreff;
  }

  // ND: Another very rough study suggests full purity up to about
  // this point. It might be ok to go to 6, but there was a mildly
  // significant dip after 4.
  //
  // XXX kludge in a minimum energy cut here for the ND, since I let
  // lots of non-muons in by accident due to confusion with EvisID(g).
  // This can be taken back out again once everything is reprocessed.
  string STD_CHI2_CUT = near?"rchi2 < 4 && fq > 1100e3":"rchi2 < 2";

  // ND: As far as I can tell, a cut of the form used at the FD does
  // nothing to improve the signal/bg at the ND. At the FD, there is a
  // drop-off in purity on the high end of the cut (but not the low end
  // -- that just cuts chimney muons, which may or may not be a good
  // thing). But at the ND, the purity is pretty much flat.
  string STD_IVDEDX_CUT =
    near?"1":"abs(fez + 62*ivdedx/2 - 8847.2) < 1000";

  countcut = 
    "ndecay == 0 && " +
    STD_POS_CUT + " && " +
    STD_IVDEDX_CUT + " && " +
    STD_CHI2_CUT;

  // For ND, not at all clear that this gives a pure sample, and surely
  // the contamination figures in consts.h are not justified for the ND
  // even if this is reasonably pure.
  CUT_PART_OK_FOR_FINDING_PACCN  = 
    STD_POS_CUT + " && " +
    STD_IVDEDX_CUT + " && " +
    STD_CHI2_CUT;

  TFile *_file0 = TFile::Open(near?rootfile3up_near:rootfile3up);
  TTree * t = (TTree *)_file0->Get("t");

  // The number of mu- stopping, regardless of what they atomicly or
  // nuclearly capture on
  const ve mum_count_ve=mucountfinalfit_cut(countcut.c_str(), !near, t);
  mum_count   = mum_count_ve.val;
  mum_count_e = mum_count_ve.err;

  c12nuc_cap = c_atomic_capture_prob*(1-f13_HP)*mum_count*capprob12;
  c13nuc_cap = c_atomic_capture_prob*f13_HP    *mum_count*capprob13;

  const string cut = 
  CUT_PART_OK_FOR_FINDING_PACCN
  + string(near?"":"&& !earlymich") +
  "&& miche < 12 && e > 4 && e < 15 "
  //" && dist < 600 " // XXX not standard, just testing
  "&& timeleft > %f && dt < %f && " + NEUTRONDEF + " <= 2";

  apply_dr_eff = strstr(NEUTRONDEF.c_str(), "near") != NULL;

  printtwice("The number of mu- stopping in the high-purity sample, "
             " regardless of what they atomicly or nuclearly capture "
             "on, is %f +- %f\n", 4, mum_count, mum_count_e);
  
  printtwice("B-12 selection efficiency is %f%%\n", 2, b12eff*100);
  printtwice("B-13 selection efficiency is %f%%\n", 2, b13eff*100);
  printtwice("Li-8 selection efficiency is %f%%\n", 2, li8eff*100);
  printtwice("Li-9 selection efficiency is %f%%\n", 2, li9eff*100);
  printtwice("Number of C-12 captures %g\n", 2, c12nuc_cap);
  printtwice("Number of C-13 captures %g\n", 2, c13nuc_cap);

  find_paccn(t);

  mn = new TMinuit(npar);
  mn->SetPrintLevel(1);
  mn->fGraphicsMode = false;
  mn->Command("SET STRATEGY 2");
  mn->SetFCN(fcn);
  int err;
  mn->mnparm(0, "p_b12" , 0.16,                 0.001, 0,   2, err);
  mn->mnparm(1, "p_b12n", 0.40,                  0.01, 0,   2, err);
  mn->mnparm(2, "p_b13" , 0.40,                  0.01, 0,
                                               unitarity?2:10, err);
  mn->mnparm(3, "p_li8" , 0.01,                  0.01, 0,   1, err);
  mn->mnparm(4, "p_li8n", 0.05,                  0.01, 0,   1, err);
  mn->mnparm(5, "p_li9" , probNineLiFromTwelveC, 0.01, 0,   1, err);
  mn->mnparm(6, "p_li9n", primaryresult/f13,     0.01, 0,   1, err);
  mn->mnparm(7, "n_n16" ,   10,                     1, 0, 1e3, err);

  mn->mnparm(8, "acc"   ,   10,                  0.01, 0,1000, err);
  mn->mnparm(9, "accn"  ,  0.1,                 0.001, 0,  10, err);
  mn->mnparm(10, "accn2", 0.01,                 0.001, 0,   1, err);

  // All constrained with pull terms
  mn->mnparm(11, "b12t",      0,  b12life_err/b12life, 0, 0, err);
  mn->mnparm(12, "b13t",      0,  b13life_err/b13life, 0, 0, err);
  mn->mnparm(13, "li8t",      0,  li8life_err/li8life, 0, 0, err);
  mn->mnparm(14, "li9t",      0,  li9life_err/li9life, 0, 0, err);
  mn->mnparm(15, "n16t",      0,  n16life_err/n16life, 0, 0, err);
  mn->mnparm(16, "neffdelta", 0,                 0.01, 0, 0, err);
  mn->mnparm(17, "paccn",     nom_paccn,      paccn_e, 0, 0, err);

  printf("Making cuts...\n");
  TFile * tmpfile = new TFile("/tmp/b12tmp.root", "recreate");
  tmpfile->cd();
  const string fullcut = Form(cut.c_str(), hightime+offset, hightime+offset);
  selt = t->CopyTree(fullcut.c_str());
  selt->Write();
  events.clear();

  float dt, mx, my, mz, fq, fqiv;
  int nn, run, trig;
  selt->SetBranchAddress("dt", &dt);
  selt->SetBranchAddress(NEUTRONDEF.c_str(), &nn);
  selt->SetBranchAddress("mx", &mx);
  selt->SetBranchAddress("my", &my);
  selt->SetBranchAddress("mz", &mz);
  selt->SetBranchAddress("fq", &fq);
  selt->SetBranchAddress("fqiv", &fqiv);
  selt->SetBranchAddress("run", &run);
  selt->SetBranchAddress("trig", &trig);

  printf("Filling in data array...");
  for(int i = 0; i < selt->GetEntries(); i++){
    selt->GetEntry(i);
    if(isibd(run, trig)) continue;
    events.push_back(ev(
      dt-offset,
      nn,
      neff_dt(near?0:fq, fqiv, mx, my, mz)

      // Apply the dr efficiency if we are using a "near" neutron variable
      *(apply_dr_eff?neff_dr_800(mx, my, mz):1)

      ));
    if(i%10000 == 9999){ printf("."); fflush(stdout); }
  }
  printf("\n%d events to be used in the fit\n", (int)events.size());
  std::sort(events.begin(), events.end(), compare_events);
  printf("Sorted by time and number of neutrons\n");

  {
    double emean = 0;
    for(unsigned int i = 0; i < events.size(); i++){
      emean += events[i].e/events.size();
      meaneff[events[i].n] += events[i].e;
      eventc[events[i].n]++;
    }
    for(unsigned int i = 0; i < 3; i++){
      if(eventc[i]) meaneff[i] /= eventc[i];
      printf("Mean neutron eff when there were %d neutrons: %f\n",
             i, meaneff[i]);
    }
    printf("Mean neutron efficiency: %f\n", emean);
  }

  {
    const int nbins = 1000;
    TH1D * n0 = new TH1D("n0","",nbins,1, 2001);
    TH1D * n1 = new TH1D("n1","",nbins,1, 2001);
    TH1D * n2 = new TH1D("n2","",nbins,1, 2001);

    selt->Draw("dt >> n0", Form("%s && %s == 0", fullcut.c_str(),
                                NEUTRONDEF.c_str()));
    selt->Draw("dt >> n1", Form("%s && %s == 1", fullcut.c_str(),
                                NEUTRONDEF.c_str()));
    selt->Draw("dt >> n2", Form("%s && %s == 2", fullcut.c_str(),
                                NEUTRONDEF.c_str()));

    printf("FIGURE const double n0contents[%d] = {", nbins+2);
    for(int i = 0; i <= nbins+1; i++)
      printf("%.9f, ", n0->GetBinContent(i));
    printf("};\n");

    printf("FIGURE const double n1contents[%d] = {", nbins+2);
    for(int i = 0; i <= nbins+1; i++)
      printf("%.9f, ", n1->GetBinContent(i));
    printf("};\n");

    printf("FIGURE const double n2contents[%d] = {", nbins+2);
    for(int i = 0; i <= nbins+1; i++)
      printf("%.9f, ", n2->GetBinContent(i));
    printf("};\n");

    printf("FIGURE const int nbin = %d;\n", n0->GetNbinsX());
    printf("FIGURE const double lowx = %.9f;\n", n0->GetBinLowEdge(1));
    printf("FIGURE const double highx = %.9f;\n",
           n0->GetBinLowEdge(n0->GetNbinsX()+1));
    printf("FIGURE const double b12eff = %.9f;\n", b12eff);
    printf("FIGURE const double b13eff = %.9f;\n", b13eff);
    printf("FIGURE const double li8eff = %.9f;\n", li8eff);
    printf("FIGURE const double li9eff = %.9f;\n", li9eff);
    printf("FIGURE const double c12nuc_cap = %.9f;\n", c12nuc_cap);
    printf("FIGURE const double c13nuc_cap = %.9f;\n", c13nuc_cap);
    printf("FIGURE const double nom_paccn = %.9f;\n", nom_paccn);
    printf("FIGURE const int npar = %d;\n", npar);
    printf("FIGURE const double meaneff[3] = { ");
    for(int i = 0; i < 3; i++) printf("%.9f, ", meaneff[i]);
    printf("};\n");

    printf("FIGURE const double mum_count = %.9f;\n", mum_count);
  }

  // Ease into the fit by fitting the major features first, then 
  // gradually relaxing constraints.
  mn->Command("FIX 1 2 3 4 5 6 7 8         12 13 14 15 16 17 18");
  mn->Command("MIGRAD");
  {
    printf("FIGURE const double pars_imed1[%d] = { ", npar);
    for(int i = 0; i < npar; i++) printf("%.9f, ", getpar(i));
    printf("};\n");
  }
  mn->Command("REL 1 2   4 5                                   ");
  mn->Command("MIGRAD");
  {
    printf("FIGURE const double pars_imed2[%d] = { ", npar);
    for(int i = 0; i < npar; i++) printf("%.9f, ", getpar(i));
    printf("};\n");
  }
  mn->Command("REL     3     6 7 8                             ");
  mn->Command("MIGRAD");
  {
    printf("FIGURE const double pars_imed3[%d] = { ", npar);
    for(int i = 0; i < npar; i++) printf("%.9f, ", getpar(i));
    printf("};\n");
  }
  mn->Command("REL                                        17 18");
  mn->Command("MIGRAD");
  {
    printf("FIGURE const double pars_imed4[%d] = { ", npar);
    for(int i = 0; i < npar; i++) printf("%.9f, ", getpar(i));
    printf("};\n");
  }
  mn->Command("REL                         12 13 14 15 16      ");

  const char * const commands[2] = { "MIGRAD", "HESSE" };
  for(int i = 0; i < 2; i++){
    printf("\n%s\n", commands[i]);
    int fails = 0;
    while(4 == mn->Command(commands[i])){
      puts(""); mn->Command("show par");
      if(++fails >= 3){
        printf("\nGiving up on %s\n", commands[i]);
        break;
      }
      else{
        printf("\nTrying %s again\n", commands[i]);
      }
    }
    puts(""); mn->Command("show par");
  }

  {
    printf("FIGURE const double pars[%d] = { ", npar);
    for(int i = 0; i < npar; i++) printf("%.9f, ", getpar(i));
    printf("};\n");
  }

  mn->SetPrintLevel(0);
  for(int i = 1; i <= 2; i++){
    int fails = 0;
    while(4 == mn->Command(Form("MINOS 10000 %d", i))){
      if(++fails >= 3){
        fprintf(stderr, "giving up on MINOS %d\n", i);
        break;
      }
      fprintf(stderr, "retrying MINOS %d\n", i);
    }
  }

  printc12b12results();
  printc13b12nresults();
  printc13b13results();

  mn->SetPrintLevel(-1);
  b13limit();
  mn->SetPrintLevel(0);
  puts("");
}

void near_fullb12_finalfit()
{
  near = true;
  NEUTRONDEF = "latennear";
  fullb12_finalfit();
}

int main()
{
  near_fullb12_finalfit();
}
