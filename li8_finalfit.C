#include <fstream>
#include "consts.h"
#include "b12cutefficiency_finalfit.out.h"
#include "carbondenominators_finalfit.out.h"
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
using std::pair;
using std::string;
using std::cin;
#include "deadtime.C" // <-- note inclusion of source

//#define HP

const double neff_dr_800_avg =
  (
  (n_c12cap - n_c12captarget*gd_fraction) * neff_dr_800_h
  +           n_c12captarget*gd_fraction  * neff_dr_800_gd
  )/n_c12cap;

using std::vector;

bool unitarity = true;

struct ev{
  double t; // time
  int n; // number of neutrons
  double e; // nominal neutron efficiency
  bool i; // is IBD?
  ev(double t_, int n_, double e_, bool i_){
    t = t_;
    n = n_;
    e = e_;
    i = i_;
  }
};

TMinuit * mn = NULL;
TTree * selt = NULL;
TF1 * zero= NULL;
TF1 * one = NULL;
TF1 * two = NULL;

vector<ev> events;

// Mean neutron efficiency for events that actually have
// zero, one or two neutrons;
double meaneff[3] = {0};
int eventc[3] = {0}; // number of events with each n count

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


#ifdef HP
const double mum_count_thisselection = mum_count;
const double mum_count_thisselection_e = mum_count_e;
const double c12nuc_cap = c_atomic_capture_prob*(1-f13)*mum_count*capprob12;
const double c13nuc_cap = c_atomic_capture_prob*f13    *mum_count*capprob13;
#else
const double c12nuc_cap = n_c12cap*livetime;
const double c13nuc_cap = n_c13cap*livetime;
const double mum_count_thiselection = n_c12cap*livetime/c_atomic_capture_prob/(1-f13)/capprob12;
const double mum_count_thiselection_e = n_c12cap_err*livetime/c_atomic_capture_prob/(1-f13)/capprob12;
#endif


// Will subtract mean muon lifetime, 2028ns, and mean transit time for
// light from B-12, 12ns. Doesn't make a real difference.
const double offset = 2.040e-3;

TH2D * hdisp = new TH2D("hdisp", "", 3, 0, 3, 40, 100, 20101);


const double lowtime = 1.0 - offset;
const double hightime = 100e3;
const double totaltime = hightime - lowtime;

const double li8eff_energy = 0.7165; // estimate from DOGS MC -- good 
const double li8eff_energy_e = 0.02; // made up!
const double li8ferr_energy = li8eff_energy_e/li8eff_energy; // made up!

const double b12energyeff_5MeV = 0.7484;  // B-12 energy cut for 5MeV
const double b12energyeff_5MeV_e = 0.0056;

const double b13energyeff = b12energyeff_5MeV * 1.014; // estimate from my MC
const double b13energyeff_e = 0.02; // BS

const double li9eff_energy = 0.6794 + 0.05; // from blessed li-9
                                            // spectrum but bumped
                                            // up since I only see
                                            // non-betan decays here.
const double li9eff_energy_e = 0.05;

// time until end of run
const double eor_eff = 1-(1-0.9709)*hightime/100e3;

const double b12likelihood_eff = 0.906; // For 0.02

const double b12eff = b12likelihood_eff * light_noise_eff *
  wholedet_dist400eff * mich_eff * eor_eff * sub_muon_eff05 * b12energyeff_5MeV;
const double b12ferr_energy = b12energyeff_5MeV_e/b12energyeff_5MeV;

const double b13eff = b12likelihood_eff * light_noise_eff *
  wholedet_dist400eff * mich_eff * eor_eff * sub_muon_eff05 * b13energyeff;
const double b13ferr_energy = b13energyeff_e/b13energyeff;

const double li8eff = b12likelihood_eff * light_noise_eff *
  wholedet_dist400eff * mich_eff * eor_eff * sub_muon_eff05 * li8eff_energy;

const double li9eff = b12likelihood_eff * light_noise_eff *
  wholedet_dist400eff * mich_eff * eor_eff * sub_muon_eff05 * li9eff_energy;

// Constant 1% absolute error assumed on neutron efficiency.
// Not a very rigourous model, but it's something.
const double neff_err = 0.01;

// Measured probablity of getting one accidental neutron.  These
// are *detected* muons, so don't apply efficiency to them.
const double paccn = 1.1e-4;

/*
 * Prints the message once with the requested floating point precision,
 * then again with all digits in the default color, starting
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


bool isibd(const int in_run, const int in_trig)
{
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

double clamp(double x, double lo, double hi)
{
  return x < lo? lo: x > hi? hi: x;
}

double zeron_f(
const double mt, const double acc,const double paccn, const double neff,
const double n_b12, const double n_b12n, const double b12t,
const double n_li8, const double n_li8n, const double li8t,
const double n_li9, const double n_li9n, const double li9t,
const double n_b13, const double b13t,
const double n_n16, const double n16t)
{
  // Be really careful. Can get zero neutrons from a real zero-neutron
  // event without an accidental or from a real one-neutron event with
  // an inefficiency and without an accidental. (Does ROOT optimize
  // these? No idea.)
  static double unpaccn = 1-paccn; // always the same
  double uneffunpaccn = (1-neff)*unpaccn; // function of pars
  return acc +
    (unpaccn*n_b12 + uneffunpaccn*n_b12n)/b12t*exp(mt/b12t) +
    // For B-13 and N-16, there are no real one-neutron events, so it
    // is easier. (Ok, N-16 might come from O-17, but it is only a
    // nuisance parameter here...)
    unpaccn*n_b13/b13t*exp(mt/b13t)
  + (unpaccn*n_li8 + uneffunpaccn*n_li8n)/li8t*exp(mt/li8t) +
    (unpaccn*n_li9 + uneffunpaccn*n_li9n)/li9t*exp(mt/li9t)
  + unpaccn*n_n16/n16t*exp(mt/n16t)
    ;
}

double onen_f(
const double mt,const double accn,const double paccn, const double neff,
const double n_b12, const double n_b12n, const double b12t,
const double n_li8, const double n_li8n, const double li8t,
const double n_li9, const double n_li9n, const double li9t,
const double n_b13, const double b13t,
const double n_n16, const double n16t)
{
  return accn +
  // Can get one neutron from a real one-neutron event that is efficient
  // and doesn't have an accidental, or that is inefficient and does
  // have an accidental, or from a real zero-neutron event that is
  // inefficient.
  ((neff*(1-paccn)+(1-neff)*paccn)*n_b12n+paccn*n_b12)/b12t*exp(mt/b12t)+
  // Can only get B-13 with a neutron from an inefficiency, assuming
  // there's no O-16(mu,ppn)B-13 to speak of.
  paccn*n_b13/b13t*exp(mt/b13t)
 +((neff*(1-paccn)+(1-neff)*paccn)*n_li8n+paccn*n_li8)/li8t*exp(mt/li8t)
 +((neff*(1-paccn)+(1-neff)*paccn)*n_li9n+paccn*n_li9)/li9t*exp(mt/li9t)
  + paccn*n_n16/n16t*exp(mt/n16t)
  ;
}

double twon_f(
const double mt,const double accn2,const double neff,const double paccn,
const double n_b12n, const double b12t,
const double n_li8n, const double li8t,
const double n_li9n, const double li9t)
{
  return accn2 +
  // Can get two neutrons from a real one-neutron event that is
  // efficient and has an accidental. Note that we *must* handle this
  // case for the above normalization component of the likelihood to be
  // correct.
  neff*paccn*(n_b12n/b12t*exp(mt/b12t)
            + n_li8n/li8t*exp(mt/li8t)
            + n_li9n/li9t*exp(mt/li9t)
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
               neffdelta = par[16];

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

void fcn(int & npar, double * gin, double & like, double *par, int flag)
{
  DECODEPARS;

  const double norm =
      (n_b12 + n_b12n)*(exp(-lowtime/b12t) - exp(-hightime/b12t))
    + (n_li8 + n_li8n)*(exp(-lowtime/li8t) - exp(-hightime/li8t))
    + (n_li9 + n_li9n)*(exp(-lowtime/li9t) - exp(-hightime/li9t))
    +      n_b13      *(exp(-lowtime/b13t) - exp(-hightime/b13t))
    +      n_n16      *(exp(-lowtime/n16t) - exp(-hightime/n16t))
    + totaltime*(acc+accn+accn2);

  like = norm;

  {
    static int dotcount = 0;
    if(mn->fCfrom != "CONtour   " && dotcount++%100 == 0){
      printf(".");  fflush(stdout); }
  }

  for(unsigned int i = 0; i < events.size(); i++){
    const double mt = -events[i].t;
    const double neff = clamp(events[i].e + neffdelta, 0.00, 1.00);

    double f = 0;
    switch(events[i].n){
      case 0:
        f = zeron_f(mt, acc, paccn, neff, n_b12, n_b12n, b12t,
                    n_li8, n_li8n, li8t, n_li9, n_li9n, li9t, n_b13, b13t, n_n16, n16t);
        break;
      case 1:
        f = onen_f(mt, accn, paccn, neff, n_b12, n_b12n, b12t,
                   n_li8, n_li8n, li8t, n_li9, n_li9n, li9t, n_b13, b13t, n_n16, n16t);
        break;
      case 2:
        f = twon_f(mt, accn2, neff, paccn, n_b12n, b12t, n_li8n, li8t, n_li9n, li9t);
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
        + pow((li9t - li9life)/li9life_err, 2)
        + pow((n16t - n16life)/n16life_err, 2)
          // and the neutron efficiency
        + pow(neffdelta/neff_err, 2);

  // pull terms for Li-9 from the betan analysis. Assume zero production
  // with a neutron so as not to double count. Since IBD candidates are
  // cut, this is only the non-betan rate. The multiper is to account
  // for the uncertainty in the Li-9 energy cut.
  static const double energymultiplier = 1 + li9eff_energy_e/li9eff_energy;
  like += pow( p_li9n/(1-0.508)/0.44e-4/energymultiplier, 2)
        + pow((p_li9 /(1-0.508) - 2.4e-4)/0.9e-4/energymultiplier, 2);
  
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

  static const double likeoffset = 252637;
  
  like += likeoffset;
}

double dispf0(double * x, double * par)
{
  DECODEPARS;
  return hdisp->GetYaxis()->GetBinWidth(1)*
    zeron_f(-x[0], acc, paccn, meaneff[0]+neffdelta, n_b12, n_b12n,
            b12t, n_li8, n_li8n, li8t, n_li9, n_li9n, li9t, n_b13, b13t, n_n16, n16t);
}

double dispf1(double * x, double * par)
{
  DECODEPARS;
  return hdisp->GetYaxis()->GetBinWidth(1)*
    onen_f(-x[0], accn, paccn, meaneff[1]+neffdelta, n_b12, n_b12n,
           b12t, n_li8, n_li8n, li8t, n_li9, n_li9n, li9t, n_b13, b13t, n_n16, n16t);
}

double dispf2(double * x, double * par)
{
  DECODEPARS;
  return hdisp->GetYaxis()->GetBinWidth(1)*
    twon_f(-x[0], accn2, paccn, meaneff[2]+neffdelta, n_b12, b12t,
           n_li8, li8t, n_li9, li9t);
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

  printtwice("\nTECHNOTE 10.1: %s raw %f +%f %f\n", 0, iname,
    fit, ferrorfitup*fit, ferrorfitlo*fit);

  const double like_central = fit/eff/mum_count_thiselection * 100/ifrac;
  const double staterrup = ferrorfitup*like_central;
  const double staterrlo = ferrorfitlo*like_central;
  const double muerr = mum_count_thiselection_e/mum_count_thiselection*like_central;
  const double err = ferr_energy*like_central;
  const double toterrup = sqrt(pow(ferrorfitup,2) +
                               pow(mum_count_thiselection_e/mum_count_thiselection,2) +
                               pow(ferr_energy, 2))*like_central;
  const double toterrlo = -sqrt(pow(ferrorfitlo,2) +
                               pow(mum_count_thiselection_e/mum_count_thiselection,2) +
                               pow(ferr_energy, 2))*like_central;
  printtwice("\n%s, eff corrected, percent per C mu- stop "
    "%f +%f %f(fit) +-%f(mu count) +-%f(Li-8 eff), "
    "+%f %f(total),  +%f -%f (non-fit)\n",
    prec1, iname, like_central, staterrup, staterrlo, muerr, err,
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

  printtwice("\nTECHNOTE results.tex %s: percent per nuclear mu- capture "
         "%f +%f %f(fit) +-%f(mu count) +-%f(eff) +-%f(cap frac), "
         "+%f %f(total),  +%f -%f (non-fit, a.k.a. syst)\n",
         prec2, 
         !strcmp(iname, "C-12 -> Li-8")?"probEightLiFromTwelveC":
         !strcmp(iname, "C-13 -> Li-8+n")?"probEightLiFromThirteenC":"???",
         like_central_percap, staterr_percapup, staterr_percaplo,
         muerr_percap, err_percap, capfracerr_percap,
         toterr_percapup, toterr_percaplo,
         sqrt(pow(toterr_percapup,2) - pow(staterr_percapup,2)),
         sqrt(pow(toterr_percaplo,2) - pow(staterr_percaplo,2))
         );

  if(!strcmp(iname, "C-13 -> Li-8+n")){
    printf("const double probEightLiFromThirteenC = %f;\n",
           like_central_percap/100.);
    printf("const double probEightLiFromThirteenC_statup = %f;\n",
           staterr_percapup/100.);
    printf("const double probEightLiFromThirteenC_statlo = %f;\n",
           staterr_percaplo/100.);
  }
  else if(!strcmp(iname, "C-12 -> Li-8")){
    printf("const double probEightLiFromTwelveC = %f;\n",
           like_central_percap/100.);
    printf("const double probEightLiFromTwelveC_statup = %f;\n",
           staterr_percapup/100.);
    printf("const double probEightLiFromTwelveC_statlo = %f;\n",
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

  printtwice("\n10^3/s: %f +%f %f(fit) +-%f(mu count) +-%f(eff), "
         "+-%f(lifetime) +%f %f(total)  +%f -%f\n",
         prec3,
         like_central_rate, staterr_rateup, staterr_ratelo,
         muerr_rate, err_rate, lifetimeerr_rate, toterr_rateup,
         toterr_ratelo,
         sqrt(pow(toterr_rateup, 2) - pow(staterr_rateup,2)),
         sqrt(pow(toterr_ratelo, 2) - pow(staterr_ratelo,2))
         );
  puts("");
}

void printli8results()
{
  results("C-12 -> Li-8", 3, 1-f13, li8eff, li8ferr_energy,
    capprob12*c_atomic_capture_prob, errcapprob12, lifetime_c12,
    lifetime_c12_err, c12nuc_cap, 3, 2, 2);
}

void printli8nresults()
{
  results("C-13 -> Li-8+n", 4, f13, li8eff, li8ferr_energy,
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


void li8_finalfit(const char * const cut =
#ifdef HP
  "mx**2+my**2 < 1050**2 && mz > -1175 && "
  "abs(fez + 62*ivdedx/2 - 8847.2) < 1000 && rchi2 < 2 && "
#endif
"b12like < 0.02 && timeleft > %f && miche < 12 && !earlymich && "
"e > 5 && e < 14 && dist < 400 && dt < %f && latennear <= 2")
{
  printtwice("B-12 selection efficiency is %f%%\n", 2, b12eff*100);
  printtwice("B-13 selection efficiency is %f%%\n", 2, b13eff*100);
  printtwice("Li-8 selection efficiency is %f%%\n", 2, li8eff*100);
  printtwice("Li-9 selection efficiency is %f%%\n", 2, li9eff*100);
  printtwice("Number of C-12 captures %g\n", 2, c12nuc_cap);
  printtwice("Number of C-13 captures %g\n", 2, c13nuc_cap);

  TFile *_file0 = TFile::Open(rootfile3up);
  TTree * t = (TTree *)_file0->Get("t");
 
  const int npar = 17;
  mn = new TMinuit(npar);
  mn->SetPrintLevel(-1);
  mn->fGraphicsMode = false;
  mn->Command("SET STRATEGY 2");
  mn->SetFCN(fcn);
  int err;
  mn->mnparm(0, "p_b12", 0.2,  0.001, 0, 2, err);
  mn->mnparm(1, "p_b12n",0.5,  0.01,  0, 2, err);
  mn->mnparm(2, "p_b13", 0.2,  0.01,  0, 2, err);
  mn->mnparm(3, "p_li8", 0.005,  0.001, 0, 2, err);
  mn->mnparm(4, "p_li8n",0.05,  0.01, 0, 2, err);
  mn->mnparm(5, "p_li9", 0.00,  0.01, 0, 2, err);
  mn->mnparm(6, "p_li9n",0.024,  0.01, 0, 2, err);
  mn->mnparm(7, "n_n16",   1,  1, 0, 1e3, err);
  mn->mnparm(8, "acc",   1,    0.01, 0, 100, err);
  mn->mnparm(9, "accn",  0.1,  0.001, 0, 10, err);
  mn->mnparm(10, "accn2", 0.01, 0.001, 0, 1, err);

  // Sometimes the B-12 lifetime spins out of control.  Not sure why,
  // but constrain it to be somewhat reasonable.
  mn->mnparm(11, "b12t", 0,  b12life_err/b12life, -5, +5, err);
  mn->mnparm(12, "b13t", 0,  b13life_err/b13life, -5, +5, err);
  mn->mnparm(13, "li8t", 0,  li8life_err/li8life, -5, +5, err);
  mn->mnparm(14, "li9t", 0,  li9life_err/li9life, -5, +5, err);
  mn->mnparm(15, "n16t", 0,  n16life_err/n16life, -5, +5, err);
  mn->mnparm(16, "neffdelta", 0,  neff_err, 0, 0, err);

  printf("Making cuts...\n");
  TFile * tmpfile = new TFile("/tmp/b12tmp.root", "recreate");
  tmpfile->cd();
  selt = t->CopyTree(Form(cut, hightime+offset, hightime+offset));
  selt->Write();
  events.clear();

  float dt, mx, my, mz, fq;
  int nn, run, trig;
  selt->SetBranchAddress("dt", &dt);
  selt->SetBranchAddress("latennear", &nn);
  selt->SetBranchAddress("mx", &mx);
  selt->SetBranchAddress("my", &my);
  selt->SetBranchAddress("mz", &mz);
  selt->SetBranchAddress("fq", &fq);
  selt->SetBranchAddress("run", &run);
  selt->SetBranchAddress("trig", &trig);

  printf("Filling in data array...");
  for(int i = 0; i < selt->GetEntries(); i++){
    selt->GetEntry(i);
    if(isibd(run, trig)) continue;
    events.push_back(ev(dt-offset, nn,
      neff_dr_800_avg*eff(fq, mx, my, mz), isibd(run, trig)));
    hdisp->Fill(nn, dt-offset);
    if(i%10000 == 9999){ printf("."); fflush(stdout); }
  }
  puts("");

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

  // I know that MIGRAD often fails, so instead of doing MINIMIZE, do
  // SIMPLEX in the first place to save a little time and get it to
  // converge. I think SIMPLEX is needed because the way I implement the
  // unitarity constraint is really harsh on MIGRAD.
  const char * const commands[3] = { "SIMPLEX", "MIGRAD", "HESSE" };
  for(int i = 0; i < 3; i++){
    printf("\n%s", commands[i]);
    mn->Command(commands[i]);
    puts(""); mn->Command("show par");
  }

  mn->SetPrintLevel(0);
  mn->Command("MINOS 2000 4 5");

  for(int cpar = 3; cpar <= 4; cpar++){
    int tries = 0;
    // Check for MINOS failure
    while(sani_minos(mn->fErn[cpar]) == 0 || 
          sani_minos(mn->fErp[cpar]) == 0){
      if(tries++ > 3) break;
      mn->Command(Form("MINOS 2000 %d", cpar+1));
    }
  }

  puts("");
  mn->Command("show min");

  zero= new TF1("zero",dispf0, 0, 100e3, npar);
  one = new TF1("one", dispf1, 0, 100e3, npar);
  two = new TF1("two", dispf2, 0, 100e3, npar);

  printli8results();
  printli8nresults();

  for(int i = 0; i < npar; i++){
    zero->SetParameter(i, getpar(i));
    one ->SetParameter(i, getpar(i));
    two ->SetParameter(i, getpar(i));
  }
}
