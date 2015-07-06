#include "unistd.h"
#include "consts.h"
#include <string>
#include "math.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TMinuit.h"
#include "TTree.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TROOT.h"
#include <stdio.h>

using std::string;

// number of neutrons.  Zero for the primary analysis.  One to determine 
// gamma background from C-13 -> B-12n
const int nn = 0;


//#define HP

TCanvas * c2 = new TCanvas("c2", "", 0, 0, 500, 400);

TH1D * bg     = new TH1D("bg"    , "", 600, 0.7, 60.7);
TH1D * corrbg = new TH1D("corrbg", "", 600, 0.7, 60.7);
TH1D * ehist  = new TH1D("ehist" , "", 600, 0.7, 60.7);

TF1 *gg = NULL;

const double neff = 0.858;
const double b12n_subfactor = (1-neff)/neff;

const double muClife = 2028.; // muon lifetime in C-12, ns
const double hn_e = 2.224573;
const double gn_e1 = 7.937; // Gd-157
const double gn_e2 = 8.536; // Gd-155

const double mulife = 2196.9811;
const double muClife_err = 2.;
const double capprob12 = 1-muClife/mulife;
const double errcapprob12 = (1-(muClife+muClife_err)/mulife)/2
                           -(1-(muClife-muClife_err)/mulife)/2;

#ifdef HP
  const double No16cap_betan = n_o16cap_betan_hp*livetime;
  const double Nc12captarget = n_c12captarget_hp*livetime;
  const double Nc12cap       = n_c12cap_hp*livetime;
  const double Nc13cap       = n_c13cap_hp*livetime;
#else
  const double No16cap_betan = n_o16cap_betan*livetime;
  const double Nc12captarget = n_c12captarget*livetime;
  const double Nc12cap       = n_c12cap*livetime;
  const double Nc13cap       = n_c13cap*livetime;
#endif

const double li8life = 839.9;

// Possible background from li-8 gammas, particularly at 980.8keV
const double li8lowt = 300, li8hight = 5*li8life;

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

double corrmiche(const double e, const double me, const double mt)
{
  const double early = e /(logis(me, 0.402/hn_e, 203., 25.0) + 1);
  const double late  = e /(logis(me, 0.559/hn_e, 286., 73.0) + 1);
  return early + (late - early) * (mt - 3625.)/1250.;
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

// eff: the efficiency
// energy: the gamma line energy
// n: the raw number of fitted events
// nelo: the lower error on that
// neup: the upper error on that
// sub: the raw number of B-12n events to subtract
// sublo: the lower error on that
// subup: the upper error on that
// The neutron efficiencies are applied here, so don't put in
// the B-12n number corrected for that.
// B-12n numbers do nothing unless nn == 0, so this can be used
// unmodified to *get* the B-12n numbers.
void print_results(const double eff, const double energy,
                   const double n, const double nelo, const double neup,
                   const double sub, const double sublo,
                                     const double subup)
{
  const double n_sub = n - (nn == 0)*sub*b12n_subfactor,
    nelo_sub=sqrt(nelo*nelo + (nn == 0)*pow(sublo*b12n_subfactor,2)),
    neup_sub=sqrt(neup*neup + (nn == 0)*pow(subup*b12n_subfactor,2));

  printtwice("\n%.0fkeV fitted number of events, raw: "
             "%f %f +%f\n", 1, energy, n, nelo, neup);

  printtwice("\n%.0fkeV fitted number of events, B12-n subbed: "
             "%f %f +%f\n", 1, energy, n_sub, nelo_sub, neup_sub);

  const double n_ec = n_sub/eff, nelo_ec = nelo_sub/eff,
                                 neup_ec = neup_sub/eff;

  /*printtwice("\n%.0fkeV fitted number of events, eff corrected: "
             "%f %f +%f\n", 1, energy, n_ec, nelo_ec, neup_ec); */

  /*printtwice("\n%.0fkeV per C-12 nuclear capture: "
             "(%f %f +%f)%%\n", 2, energy,
             100*n_ec   /Nc12cap,
             100*nelo_ec/Nc12cap,
             100*neup_ec/Nc12cap);*/

  const double p_b12_from_c12 = 0.1735;

  printtwice("\n%.0fkeV per bound B-12 production: "
             "(%f %f +%f)%%\n", 2, energy,
             100*n_ec   /Nc12cap/p_b12_from_c12,
             100*nelo_ec/Nc12cap/p_b12_from_c12,
             100*neup_ec/Nc12cap/p_b12_from_c12);

  const double ratemult = capprob12/(Nc12cap*muClife)*1e6;

  printtwice("\n%.0fkeV rate: %f %f +%f e-3\n", 3, energy,
             n_ec*ratemult,
             nelo_ec*ratemult,
             neup_ec*ratemult);
  puts("");
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

TF1 * drawpeak(TF1 * gg, const char * const peak)
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

  for(int i = 2; i < 4; i++)
    mypeak->SetParameter(i,
      gg->GetParameter(i+gg->GetParNumber(Form("n%s", peak))-2));

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
  gg->SetParameters(par);
  //gg->Draw("same"); c2->Update();  c2->Modified();
  for(int i = 0; i < ehist->GetNbinsX(); i++){
    const double bincenter = ehist->GetBinCenter(i);
    const double binlo = ehist->GetBinLowEdge(i);
    const double binhi = ehist->GetBinLowEdge(i+1);
    if(bincenter < 0.7 || bincenter > 15) continue;
    const double data = ehist->GetBinContent(i);

    const double theory = gg->Integral(binlo, binhi)/(binhi-binlo);

    chi2 += theory - data;
    if(data > 0 && theory > 0) chi2 += data*log(data/theory);
  }
  chi2 *= 2;
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

void b12gammafinalfit(const int region = 1)
{
  const double lowt   = region == 0? 4000 :region == 1?  3008: 2016;
  const double hight = 5024.; // ns
  const double highfq = 215; // MeV

#ifdef HP
  // From muon counting.  Exact.
  const double fq_eff = 481414./1628874;
#else
  // From B-12-like event counting, with a bit over 1% stat error. Since
  // the most precise output of this fit has 10% stat errors, I won't
  // obsess over that.
  const double fq_eff = 0.3857;
#endif

  const double b12eff = 1
    * 0.981  // Subsequent muon veto efficiency 
    * eff_eor_b12 // timeleft cut
    * (b12ecutlow == 3?0.9251:
       b12ecutlow == 4?0.8504:
       (exit(1),1)) // B-12 energy cut
    * (distcut == 400?wholedet_dist400eff:(exit(1),1))
    * (exp(-b12lowt *log(2)/b12hl)
      -exp(-b12hight*log(2)/b12hl)) // B-12 beta decay time
    * fq_eff
  ;

  const double gammatimecut_eff = 
    (exp(-lowt/muClife)-exp(-hight/muClife));
  const double eff = b12eff * gammatimecut_eff;

  printtwice("Gamma time efficiency: %f%%\n", 1, 100*gammatimecut_eff);
  printtwice("Efficiency: %f%%\n", 1, 100*eff);

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

    // [30] is the normalization of the Gd-n distribution, with a 
    // parameterization derived from the histogram in doc-5593-v3
    + " + ([30]*(1-[28]))*( "

    // normalization of the gaussian peak
    "0.725471/(sqrt([23]**2*8.19701+([24]*8.19701)**2+[25]**2 + 0.1091)*"
    "sqrt(2*TMath::Pi()))"

    // main gaussian.  The width is the quadrature sum of the intrinsic
    // width (sqrt(0.1091)) and the resolution
     "*exp(-(((x-8.19701*[1])/"
     "sqrt([23]**2* 8.19701+([24]*8.19701)**2+[25]**2 + 0.1091))**2)/2)"

    // Then we have a parametrization for the low-energy tail
    // First a double logisitic that goes up from zero at 5MeV, then
    // down in the middle of the gaussian peak
    "+0.0372863/(1+exp(-(x-5.0757)/0.332407))/(1+exp((x-8)/0.1))"

    // Then a logistic that handles the stuff below 5, and goes
    // away in the middle of the gaussian peak.
    "+0.0149869/(1+exp((x-8)/0.1))"

    ")"
     ).c_str(), 0,15);

  TFile * f = 
//#define FAST
#ifdef FAST
  new TFile(Form("b12gamma.sel%d.root", region), "read");
  if(!f || f->IsZombie()) f =
#endif
  new TFile(rootfile3up, "read");
  TTree * t = (TTree *)f->Get("t");

  bg->Sumw2();
  corrbg->Sumw2();
  ehist->Sumw2();

  TFile * tempfile = new
    TFile(Form("/tmp/b12gammafit.%d.root", getpid()), "recreate");

  printf("Pre-selecting...\n");
  TTree * seltree = t->CopyTree(Form(
         #ifdef HP
          "mx**2+my**2 < 1050**2 && mz > -1175 && "
          "abs(fez + 62*ivdedx/2 - 8847.2) < 1000 && chi2 < 2 && "
         #endif
         "!earlymich && "
         "latennear==%d && "
         "e > %f && e < 15 && "
         "dist < %f && "
         "fq < %f && "
         "micht >= %f && micht < %f"
         , nn, b12ecutlow, distcut, fq_per_mev*highfq, lowt, hight));
  seltree->Write();

  printf("%d events in t, %d in seltree\n",
         t->GetEntries(), seltree->GetEntries());
  if(seltree->GetEntries() == 0){ printf("No events in seltree\n"); exit(1);}

  int ndecay;
  float dt, timeleft, miche, micht, fq; // yes, all floats
  #define SBA(x) seltree->SetBranchAddress(#x, &x);
  SBA(ndecay); 
  SBA(dt); 
  SBA(timeleft); 
  SBA(miche); 
  SBA(micht); 
  SBA(fq); 

  printf("Drawing...\n");
  for(int i = 0; i < seltree->GetEntries(); i++){
    seltree->GetEntry(i);

    if(timeleft > b12hight && dt > b12lowt && dt < b12hight && ndecay == 0)
      ehist ->Fill(corrmiche(miche, fq/fq_per_mev, micht));

    // "ndecay == 0" how to handle this? Really need the first event in
    // lots of windows, which isn't convenient -- I think this is close
    // enough
    if(timeleft > acchight && dt > acclowt && dt < acchight)
      bg    ->Fill(corrmiche(miche, fq/fq_per_mev, micht));

    if(timeleft > li8hight && dt > li8lowt && dt < li8hight && ndecay == 0)
      corrbg->Fill(corrmiche(miche, fq/fq_per_mev, micht));
  }

  if(ehist->GetEntries() == 0){ printf("No signal events\n"); exit(1); }

  bg->SetLineColor(kRed);
  bg->SetMarkerColor(kRed);
  corrbg->SetLineColor(kViolet);
  corrbg->SetMarkerColor(kViolet);
  const double bgscaleb12 = (b12hight - b12lowt)/(acchight-acclowt)
           *eff_eor_b12/eff_eor_acc;
  const double bgscaleli8 = (li8hight - li8lowt)/(acchight-acclowt)
           *eff_eor_li8/eff_eor_acc;
  bg->Scale(bgscaleb12);

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
    (exp(-b12lowt/li8life) - exp(-b12hight/li8life))/
    (exp(-li8lowt/li8life) - exp(-li8hight/li8life))/eff_eor_li8;

  TF1 * corrbgfit = new TF1("corrbgfit", "gaus(0)", 0.7, 2);
  corrbgfit->SetLineColor(kViolet);
  const double li8gamma = 0.9808;
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
  corrbg->Fit("corrbgfit", "i");
  corrbgfit->SetNpx(300);

  const char * ggpars[31] = { "accidentals", "energyscale", "unused",
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

  TMinuit * mn = new TMinuit(gg->GetNpar());
  mn->Command("Set print 0");
  mn->fGraphicsMode = false;
  mn->SetFCN(fcn);

  for(int i = 0; i < gg->GetNpar(); i++) gg->SetParName(i, ggpars[i]);

  for(int mni = 0; mni < gg->GetNpar(); mni++){
    int err = 0;
    mn->mnparm(mni, gg->GetParName(mni), 0, 0.1, 0, 0, err);
  }

  gg->SetNpx(400);

  for(int i = 15; i <= 22; i++)
    fixat(mn, 1+i, bggg->GetParameter(i-15));

  fixat(mn, 1+23, 0.077);
  fixat(mn, 1+24, 0.018);
  fixat(mn, 1+25, 0.1);

  fixat(mn, 1+0, 0);
  fixat(mn, 1+1, 1);
  fixat(mn, 1+2, 0.22);

  fixat(mn, 1+4, 0.95314);
  fixat(mn, 1+6, 1.67365);
  fixat(mn, 1+8, 2.6208);

  fixat(mn, 1+9, 0);
  fixat(mn, 1+10,0);

  fixat(mn, 1+12,3.759);
  fixat(mn, 1+13, 0);
  fixat(mn, 1+14,9.040);

  fixat(mn, 1+26, max(0, corrbgfit->GetParameter(0)));
  fixat(mn, 1+27, corrbgfit->GetParameter(1));
  fixat(mn, 1+28, corrbgfit->GetParameter(2));

  // contamination by Hn and Gdn events.
  if(nn == 0){
    const double hn_t = 179.;

    // from my own measurements of B-12n and Li-8n from C-13
    const double Pn = 0.516 + 0.054 * b12hl/li8life;
    const double ntrueb12n = Nc13cap*Pn;
    const double nobsb12n = ntrueb12n*b12eff;

    printf("Getting muon stop positions...\n");
    const double targfrac = 
      double(t->GetEntries(Form(
          "mx**2+my**2 < 1154**2 && "
          "abs(mz) < 1229+4+0.03*(1154-sqrt(mx**2+my**2)) && "
        #ifdef HP
          "mx**2+my**2 < 1050**2 && mz > -1175 && "
          "abs(fez + 62*ivdedx/2 - 8847.2) < 1000 && chi2 < 2 && "
        #endif
         "ndecay == 0 && fq < %f", fq_per_mev*highfq)))/
      t->GetEntries(Form(
        #ifdef HP
          "mx**2+my**2 < 1050**2 && mz > -1175 && "
          "abs(fez + 62*ivdedx/2 - 8847.2) < 1000 && chi2 < 2 && "
        #endif
         "ndecay == 0 && fq < %f", fq_per_mev*highfq));

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
  else{
    fixat(mn, 1+gg->GetParNumber("n_totn"), 0);
    fixat(mn, 1+gg->GetParNumber("e_hn"), hn_e);
    fixat(mn, 1+gg->GetParNumber("frac_hn"), 0);
  }

  for(int i = 3; i <= 13; i+=2)
    mn->Command(Form("set limits %d 0. 100.", i+1));

  mn->Command("show par");

  ehist->Draw("e");
  ehist->GetXaxis()->SetRangeUser(0.7, 12);

  printf("Fitting signal, step 1/2...\n");
  mn->Command("MIGRAD");
  mn->Command("Release 26");
  mn->Command(Form("Set limits 26 0 %f", 1.92530e-01 + 1.81760e-02));
  if(ehist->GetEntries() < 50) fixat(mn, 26, 0.186);
  printf("Fitting signal, step 2/2...\n");
  mn->Command("MIGRAD");

  bg->Draw("histsame");
  corrbg->Draw("samee");
  bggg->Draw("same");
  corrbgfit->Draw("same");


  mn->Command("Set strategy 2");
  if(region == 0){
    mn->Command("MINOS 10000 4");
    mn->Command("MINOS 10000 6");
  }
  else{
    mn->Command("MINOS 10000 8");
    mn->Command("MINOS 10000 12");
  }

  for(int i = 0; i < gg->GetNpar(); i++){
    double p, dummy;
    mn->GetParameter(i, p, dummy);
    gg->SetParameter(i, p);
  }

  mn->Command("showmin");

  const double a = gg->GetParameter("er_a");
  const double b = gg->GetParameter("er_b");
  const double c = gg->GetParameter("er_c");

  TF1 * nh_peak = drawpeak(gg, "_hn");
  nh_peak->SetLineStyle(7);

  TF1 * gh_peak = drawnpeaks(gg);
  gh_peak->SetLineStyle(7);

  if(region == 0) {
    {
      drawpeak(gg, "1");


      // The number of events in the peak, corrected for efficiency.
      const double nev = gg->GetParameter("n1")/ehist->GetBinWidth(1);
      const double neveup = mn->fErp[0]/gg->GetParameter("n1")*nev;
      const double nevelo = mn->fErn[0]/gg->GetParameter("n1")*nev;

      print_results(eff, 953, nev, nevelo, neveup,
                     6.075891, -2.337429, +3.095766);
    }
    {
      drawpeak(gg, "2");

      const double nev = gg->GetParameter("n2")/ehist->GetBinWidth(1);
      const double neveup = mn->fErp[1]/gg->GetParameter("n2")*nev;
      const double nevelo = mn->fErn[1]/gg->GetParameter("n2")*nev;

      print_results(eff, 1674, nev, nevelo, neveup,
                    2.555937, -1.374821, +2.107725);
    }
  }
  else if(region == 1){
    {
      drawpeak(gg, "3");

      const double nev = gg->GetParameter("n3")/ehist->GetBinWidth(1);
      const double neveup = mn->fErp[2]/gg->GetParameter("n3")*nev;
      const double nevelo = mn->fErn[2]/gg->GetParameter("n3")*nev;

      print_results(eff, 2621, nev, nevelo, neveup, 
                     1.018417, -0.722168, +1.403013);
    }
    {
      drawpeak(gg, "5");

      const double nev = gg->GetParameter("n5")/ehist->GetBinWidth(1);
      const double neveup = mn->fErp[3]/gg->GetParameter("n5")*nev;
      const double nevelo = mn->fErn[3]/gg->GetParameter("n5")*nev;

      print_results(eff, 3759, nev, nevelo, neveup, 
                    1.001790, -0.722315, +1.401882);
    }
  }
  else if(region == 2 && gg->GetParameter("n6") != 0){
    drawpeak(gg, "6");

    const double nev = gg->GetParameter("n6")/ehist->GetBinWidth(1);
    const double neveup = mn->fErp[4]/gg->GetParameter("n6")*nev;
    const double nevelo = mn->fErn[4]/gg->GetParameter("n6")*nev;

    print_results(eff, 9040, nev, nevelo, neveup, 0, 0, 0);
  }
  gg->Draw("same");
}
