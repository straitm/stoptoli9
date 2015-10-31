#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include "consts.h"
#include "carbondenominators_finalfit.out.h"
#include "noncarbondenominators_finalfit.out.h"

#include "TCanvas.h"
#include "TColor.h"
#include "TF1.h"
#include "TGaxis.h"
#include "TGraph.h"
#include "TH1.h"
#include "TLatex.h"
#include "TMarker.h"
#include "TMath.h"
#include "TMinuit.h"
#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "TRandom3.h"

using std::vector;
using std::string;
using std::cin;
using std::cout;

//#define HP

// H efficiency is lower than GC efficiency by a little 
// because some H is in the Target.
const double neff_dt_h = 
  (
  (n_c12cap - n_c12captarget)      * neff_dt_gc
+ n_c12captarget * (1-gd_fraction) * neff_dt_targ
  )/
  (n_c12cap - gd_fraction*n_c12captarget);

static const double li9t = 0.257233,
                    he8t = 0.171825,
                    n17t = 6.020366,
                    c16t = 1.077693,
                    b13t = 0.025002,
                    li11t= 0.012624;

static TRandom3 superrand;

struct parres{
  double val, ep, en, lp, ln;
  bool fix;

  // !!! parameters are assumed to *always* have limits set !!!
  parres(double _val, double _ep, double _en, double _lp, double _ln,
         bool _fix)
  {
    if(_ep == -54321 || _ep == 0) _ep = _lp - _val;
    if(_en == -54321 || _en == 0) _en = _val - _ln;
    val = _val; ep = fabs(_ep), en = fabs(_en),
    lp = _lp, ln = _ln; fix = _fix;
  }

  double rand() const
  {
    if(fix) return val;
    if(val < ln || val > lp) return val;

    double ans = 0;
    if(superrand.Rndm() < 0.5){
      do{ 
        ans = val+superrand.Rndm()*ep*1.01;
      }while(ans > lp);
      return ans;
    }
    else{
      do{
        ans = val-superrand.Rndm()*en*1.01;
      }while(ans < ln);
      return ans;
    }
  }

  double max() const
  {
    if(fix) return val;
    if(val < ln || val > lp) return val;

    return val+ep < lp? val+ep: lp;
  }

  double min() const
  {
    if(fix) return val;
    if(val < ln || val > lp) return val;

    return val-en > ln? val-en: ln;
  }
};


struct ev{
  bool ish; // is it H-n?
  double t; // time
  int period; // RRM run period
  double nt_li9, nt_he8, nt_n17, nt_c16, nt_b13, nt_li11;

  ev(bool ish_, double t_, int period_){
    ish = ish_;
    t = t_;
    period = period_;
    nt_li9 = -t/li9t; 
    nt_he8 = -t/he8t; 
    nt_n17 = -t/n17t; 
    nt_c16 = -t/c16t; 
    nt_b13 = -t/b13t; 
    nt_li11 = -t/li11t; 
  }
};

vector<ev> events;
double Geff, Heff;

bool dopull = true;  // bad global variable modified as we go
bool ifdopulldob13pull = true;

const int npar = 12;

const float dist = 300;

#ifdef HP
  const double No16cap_betan = n_o16cap_betan_hp;
  const double Nc12captarget = n_c12captarget_hp;
  const double Nc12cap       = n_c12cap_hp;
  const double Nc13cap       = n_c13cap_hp;
#else
  const double No16cap_betan = n_o16cap_betan;
  const double Nc12captarget = n_c12captarget;
  const double Nc12cap       = n_c12cap;
  const double Nc13cap       = n_c13cap;
#endif

// bn decay probability multiplied by the number of captures relative
// to C-12
const double li9ebn = 0.5080,
             he8ebn = 0.1600,
             c16ebn = 0.99  * 88.0/102.5*0.00243*No16cap_betan/Nc12cap,
             n17ebn = 0.951 * 88.0/102.5*0.00243*No16cap_betan/Nc12cap,
             b13ebn = 0.00286 * Nc13cap/Nc12cap,
             li11ebn = 0.789 * Nc13cap/Nc12cap;
             
const double distcuteffgc = 0.7487,
             distcutefftarg = 0.9202;

// 100s begin-of-run requirement taken into account here
const double rrmlivetime=rrmlivetimes[0]+rrmlivetimes[1]+rrmlivetimes[2];
const double denominator = 0.9709*rrmlivetime*Nc12cap;

/* DC3rdPub product of muon, light noise, OV, multiplicity,
   neutron (E, t, R), FV and IV efficiencies */
const double Geff_sans_prompt_or_mun =
       distcutefftarg*
       (1-4.49/100.)*
       (1-0.01/100.)*
       (1-0.06/100.)*
       (1-1.06/100.)*
       0.9829*
       (1-0.66/100.)*
       (1-0.04/100.);

const double Heff_sans_prompt_or_mun =
  (
  (Nc12cap - Nc12captarget)      * distcuteffgc
+ Nc12captarget * (1-gd_fraction) * distcutefftarg
  )/  
  (Nc12cap - gd_fraction*Nc12captarget)* 
       (1-1.25*4.49/100.)* // muon - ok, straightforwards scaling
       (1-0.01/100.)* // ? LN - same cut, but not obviously same eff
                      // however, *very* small for Gd, so...
       (1-0.06/100.)* // OV - ok, same
       (1-((8.+9.)/(2.+6.))*1.06/100.)* // Mult - I think this is valid
       0.9512*       // neutron (ANN,E,t,R) -- plot on slide 6 of 
                    // doc-5863 -- I think this is right
       (1-0.806/100.)* // FV -- doc5480
       (1-0.025/100.)* // IV prompt -- doc5813
       (1-0.014/100.); // IV delayed -- doc5813


// But not the actual gd fraction because of geometrical effects
const double expectedgdfrac = gd_fraction*Nc12captarget/Nc12cap;


double getlimlo(TMinuit * mn, int i)
{
  double answer, dum;
  int idum;
  TString sdum;
  mn->mnpout(i, sdum, dum, dum, answer, dum, idum);
  return answer;
}

double getlimup(TMinuit * mn, int i)
{
  double answer, dum;
  int idum;
  TString sdum;
  mn->mnpout(i, sdum, dum, dum, dum, answer, idum);
  return answer;
}

static double getpar(TMinuit * mn, int i)
{
  double answer, dum;
  mn->GetParameter(i, answer, dum);
  return answer;
}

void fixat(TMinuit * mn, int i, float v)
{
  mn->Command(Form("REL %d", i));
  if(getlimup(mn, i-1)) mn->Command(Form("SET LIM %d", i));
  mn->Command(Form("SET PAR %d %g", i, v));
  mn->Command(Form("FIX %d", i));
}

void fixatzero(TMinuit * mn, int i)
{
  mn->Command(Form("SET LIM %d", i));
  mn->Command(Form("SET PAR %d 0", i));
  mn->Command(Form("FIX %d", i));
}

int reactorpowerbin(const int run)
{
  static bool inited = false;
  static vector<int> on_off, off_off;
  if(!inited){
    inited = true;
    // From doc-5095-v2. I see that most are included in the on-off
    ifstream offofffile("offoff.h");
    if(!offofffile.is_open()){
      fprintf(stderr, "Could not open offoff.h\n");
      exit(1);
    }
    // From doc-5341
    ifstream onofffile ("onoff.h");
    if(!onofffile.is_open()){
      fprintf(stderr, "Could not open onoff.h\n");
      exit(1);
    }
    int r;
    while(offofffile >> r) off_off.push_back(r);
    while(onofffile  >> r) on_off.push_back(r);
    if(off_off.empty() || on_off.empty()){
      fprintf(stderr, "Off-off has %d runs, on-off has %d: Bad(?)\n",
              (int)off_off.size(), (int)on_off.size());
    }
  }

  if(std::binary_search(off_off.begin(), off_off.end(), run)) return 0;
  if(std::binary_search( on_off.begin(),  on_off.end(), run)) return 1;
  return 2;
}


static bool fcnearlystop = false;
static float fcnstopat = 0;

void fcn(int & npar, double * gin, double & chi2, double *par, int flag)
{
  const double rrmbg0 = fabs(par[0]),
               gdfracacc = fabs(par[1]),
               par_li9 = fabs(par[2]),
               par_he8 = fabs(par[3]),
               par_n17 = fabs(par[4]),
               oxygen_gd_frac = fabs(par[5]),
               par_c16 = fabs(par[6]),
               par_b13 = fabs(par[7]),
               par_li11= fabs(par[8]),
               carbon_gd_frac = fabs(par[9]),
               rrmbg1 = fabs(par[10]),
               rrmbg2 = fabs(par[11]);

  const double carbon_h_frac = 1-carbon_gd_frac,
               oxygen_h_frac = 1-oxygen_gd_frac,
               hfracacc = 1-gdfracacc;

  chi2 = 2*(denominator*Heff*(
           li9ebn*par_li9*carbon_h_frac*exp(-1e-3/li9t)+// H Li-9
           he8ebn*par_he8*carbon_h_frac*exp(-1e-3/he8t)+// H He-8
           n17ebn*par_n17*carbon_h_frac*exp(-1e-3/n17t)+// H N-17
           c16ebn*par_c16*carbon_h_frac*exp(-1e-3/c16t)+// H C-16
           b13ebn*par_b13*carbon_h_frac*exp(-1e-3/b13t)+// H B-13
          li11ebn*par_li11*carbon_h_frac*exp(-1e-3/li11t)+// H Li-11
           99.999*(rrmbg0+rrmbg1+rrmbg2)*hfracacc) + // H bg
         denominator*Geff*(
           li9ebn*par_li9*carbon_gd_frac*exp(-1e-3/li9t)+ // Gd Li-9
           he8ebn*par_he8*carbon_gd_frac*exp(-1e-3/he8t)+ // Gd He-8
           n17ebn*par_n17*carbon_gd_frac*exp(-1e-3/n17t)+ // Gd N-17
           c16ebn*par_c16*carbon_gd_frac*exp(-1e-3/c16t)+ // Gd C-16
           b13ebn*par_b13*carbon_gd_frac*exp(-1e-3/b13t)+ // Gd B-13
          li11ebn*par_li11*carbon_gd_frac*exp(-1e-3/li11t)+ // Gd Li-11
           99.999*(rrmbg0+rrmbg1+rrmbg2)*gdfracacc));     // Gd bg

  // pull terms
  if(dopull){
    chi2 += pow((par_n17-0.5)/0.5, 2); // 50%+-70% for N-17
    chi2 += pow((par_c16-0.05)/0.05, 2); // 5%+-10%  for C-16
    
    // Approximation for the chi2 surface of B-13 from plain beta
    // decays, an independent measurement. Probably best to only use
    // this when specifically looking at B-13. Note also that this is
    // valid up to about 500% probability of B-13, and turns over at
    // about 800%, so you will have fit failures if you push it too
    // hard.
    if(ifdopulldob13pull){
      const double b13p = par_b13 < 0.6?
                       // portion mostly controlled by the data
                       9.12889e-05
                           +1.3982*par_b13
                          +15.7956*pow(par_b13, 2)
                          -246.523*pow(par_b13, 3)
                          +1732.19*pow(par_b13, 4)
                          -5859.71*pow(par_b13, 5)
                          +9179.73*pow(par_b13, 6)
                          -5134.55*pow(par_b13, 7)
                                       :
                       // portion mostly controlled by unitarity
                           41.7936-4.12173965962690403e-01
                          -150.075*par_b13
                          +118.477*pow(par_b13, 2)
                          +62.4278*pow(par_b13, 3)
                          -6.17429*pow(par_b13, 4)
                                       ;
      chi2 += b13p;
    }
  }

  if(fcnearlystop && chi2 > fcnstopat) return;

  static const double li9ebn_t =  li9ebn/ li9t,
                      he8ebn_t =  he8ebn/ he8t,
                      n17ebn_t =  n17ebn/ n17t,
                      c16ebn_t =  c16ebn/ c16t,
                      b13ebn_t =  b13ebn/ b13t,
                     li11ebn_t = li11ebn/li11t;


  const double li9_ebn_tp_gd =  li9ebn_t*par_li9*carbon_gd_frac,
               he8_ebn_tp_gd =  he8ebn_t*par_he8*carbon_gd_frac,
               n17_ebn_tp_gd =  n17ebn_t*par_n17*oxygen_gd_frac,
               c16_ebn_tp_gd =  c16ebn_t*par_c16*oxygen_gd_frac,
               b13_ebn_tp_gd =  b13ebn_t*par_b13*carbon_gd_frac,
              li11_ebn_tp_gd = li11ebn_t*par_li11*carbon_gd_frac;

  const double li9_ebn_tp_h =  li9ebn_t*par_li9*carbon_h_frac,
               he8_ebn_tp_h =  he8ebn_t*par_he8*carbon_h_frac,
               n17_ebn_tp_h =  n17ebn_t*par_n17*oxygen_h_frac,
               c16_ebn_tp_h =  c16ebn_t*par_c16*oxygen_h_frac,
               b13_ebn_tp_h =  b13ebn_t*par_b13*carbon_h_frac,
              li11_ebn_tp_h = li11ebn_t*par_li11*carbon_h_frac;

  static const double eff_gd_livefrac[3] = {
    rrmlivetimes[0]/rrmlivetime*Geff,
    rrmlivetimes[1]/rrmlivetime*Geff,
    rrmlivetimes[2]/rrmlivetime*Geff };

  static const double eff_h_livefrac[3] = {
    rrmlivetimes[0]/rrmlivetime*Heff,
    rrmlivetimes[1]/rrmlivetime*Heff,
    rrmlivetimes[2]/rrmlivetime*Heff };

  const double eff_h_bg [3]={ rrmbg0*Heff*hfracacc,
                              rrmbg1*Heff*hfracacc,
                              rrmbg2*Heff*hfracacc };
  const double eff_gd_bg[3]={ rrmbg0*Geff*gdfracacc,
                              rrmbg1*Geff*gdfracacc,
                              rrmbg2*Geff*gdfracacc };

  for(unsigned int i = 0; i < events.size(); i++){
    ev * e = &(events[i]);
    const int per = e->period;
    double f;
    if(e->ish)
      f = eff_h_livefrac[per]*(
          b13_ebn_tp_h*exp(e->nt_b13)+ li11_ebn_tp_h*exp(e->nt_li11)+
          li9_ebn_tp_h*exp(e->nt_li9)+ he8_ebn_tp_h*exp(e->nt_he8)+
          n17_ebn_tp_h*exp(e->nt_n17)+ c16_ebn_tp_h*exp(e->nt_c16)
        )+
        eff_h_bg[per]; // H bg
    else
      f = eff_gd_livefrac[per]*(
          b13_ebn_tp_gd*exp(e->nt_b13)+ li11_ebn_tp_gd*exp(e->nt_li11)+
          li9_ebn_tp_gd*exp(e->nt_li9)+ he8_ebn_tp_gd*exp(e->nt_he8)+
          n17_ebn_tp_gd*exp(e->nt_n17)+ c16_ebn_tp_gd*exp(e->nt_c16)
         )+
         eff_gd_bg[per]; // Gd bg

    if(f > 0){
      chi2 -= 2*log(f);
      if(fcnearlystop && chi2 > fcnstopat) return;
    }
  }
}

void setupmn(TMinuit * mn, const double expectedgdfrac)
{
  mn->SetPrintLevel(-1);
  mn->fGraphicsMode = false;
  dopull = true;
  for(int i = 0; i < npar; i++) mn->Command(Form("REL %d", i+1));
  mn->Command("SET LIM  1 0 1e-4");
  mn->Command("SET LIM 11 0 1e-3");
  mn->Command("SET LIM 12 0 1e-3");
  if(rrmlivetimes[0] == 0) fixatzero(mn, 1);
  if(rrmlivetimes[1] == 0) fixatzero(mn, 11);
  if(rrmlivetimes[2] == 0) fixatzero(mn, 12);

  mn->Command(Form("SET PAR 2 %g", expectedgdfrac));
  mn->Command("Set LIM 2 0 0.5");

  mn->Command("SET LIM 3 0 3e-3");
  mn->Command("SET LIM 4 0 3e-3");
  mn->Command("SET LIM 5 0 10"); // allow 100% production plus
                                 // a large factor for slop.  There's
                                 // a pull term, too.
  mn->Command("SET LIM 6 0 0.5");
  mn->Command("SET LIM 7 0 1");  // as with N-17
  mn->Command("SET LIM 8 0 1");
  mn->Command("SET LIM 9 0 1");

  mn->Command(Form("SET PAR 10 %g", expectedgdfrac));
  mn->Command("Set LIM 10 0 0.5");
}

static TMinuit * make_a_tminuit()
{
  TMinuit * mn = new TMinuit(npar);
  mn->SetPrintLevel(-1);
  mn->Command("SET STRATEGY 2");
  mn->SetFCN(fcn);
  int err;
  mn->mnparm(1 -1, "bgrrm0", 1e-7,  1e-7, 0, 0, err);
  mn->mnparm(2 -1, "gdfracacc",0.38,0.01, 0, 0, err);
  mn->mnparm(3 -1, "Li-9",   1e-4,  1e-6, 0, 0, err); // n: yes
  mn->mnparm(4 -1, "He-8",   1e-4,  1e-3, 0, 0, err); // n: yes
  mn->mnparm(5 -1, "N-17",    0.5,   0.5, 0, 0, err); // n: yes
  mn->mnparm(6 -1,"ocapgdfrac",0.01,1e-2, 0, 0, err);
  mn->mnparm(7 -1, "C-16",   0.05,   0.5, 0, 0, err); // n: yes
  mn->mnparm(8 -1, "B-13",   0.1,     1, 0, 0, err); // n: no
  mn->mnparm(9 -1, "Li-11",  0.1,  3e-3, 0, 0, err); // n: no
  mn->mnparm(10-1, "gdfracsig",0.38, 0.01,0, 0, err);
  mn->mnparm(11-1, "bgrrm1", 3e-5,  1e-6, 0, 0, err);
  mn->mnparm(12-1, "bgrrm2", 6e-5,  1e-6, 0, 0, err);

  setupmn(mn, expectedgdfrac);
  return mn;
}


TF1 * make_bounds_tf1(int i, int whichh, int runperiod, int t,
                      double * eedisppar, char * functionstring,
                      float low, float high2)
{
  TF1 * f = new TF1(
    Form("eedisp%d_%d_%d_%d", i, whichh, runperiod+10, t),
    functionstring, low, high2);

  f->SetLineWidth(1);
  f->SetNpx(400);
  if(runperiod == -1) f->SetParameter(0,
     eedisppar[0]+eedisppar[10]+eedisppar[11]);
  else if(runperiod == 0) f->SetParameter(0, eedisppar[0]);
  else if(runperiod == 1) f->SetParameter(0, eedisppar[10]);
  else/*(runperiod == 2)*/f->SetParameter(0, eedisppar[11]);

  for(int j = 1; j < f->GetNpar(); j++)
    f->SetParameter(j, eedisppar[j]);

  return f;
}

void mncommand()
{
  string command;
  TMinuit * mn = make_a_tminuit();
  prompt:
  printf("MINUIT> ");
  if(!getline(cin, command) || command == "exit") goto exit;
  mn->Command(command.c_str());
  goto prompt;
  exit:
  printf("MINUIT exiting on %s\n", command == "exit"?"request":"EOF");
}

int whichh = 0;
void drawhist(TTree * tgsel, TTree * thsel,
              const vector<parres> & parsave, 
              const int nbin, const double low, const double high,
              const int nbin2, const double high2)
{
  whichh++;

  TCanvas * c1 = new TCanvas("c1", "c1", 0, 0,800, 500);
  c1->cd();
  // XXX c1->Divide(2, 2);
  // XXX c1->cd(1);
  
  double bins[nbin+nbin2+1];
  for(int i = 0; i < nbin; i++)
    bins[i] = double(i)*(high-low)/nbin; 
  for(int i = 0; i <= nbin2; i++)
    bins[nbin+i] = high + double(i)*(high2-high)/nbin2; 

  for(int runperiod = -1; runperiod <= -1 /* XXX 2*/; runperiod++){
    // XXX c1->cd(runperiod+2);

    const double livefrac =
      runperiod < 0?1:rrmlivetimes[runperiod]/rrmlivetime;

    TH1D * htmp = new TH1D("htmp", "", nbin+nbin2, bins);

    thsel->Draw("dt/1000 >> htmp",
      runperiod>=0?Form("reactorpowerbin(run)==%d", runperiod):"");

    tgsel->Draw("dt/1000 >> +htmp",
      runperiod>=0?Form("reactorpowerbin(run)==%d", runperiod):"");

    TH1D * hdisp = new TH1D(Form("hdisp%d_%d", whichh, runperiod+10),"",
      nbin+nbin2, low, (nbin+nbin2)*htmp->GetBinWidth(1));

    for(int i = 1; i <= hdisp->GetNbinsX(); i++){
      hdisp->SetBinContent(i, htmp->GetBinContent(i) *
                              htmp->GetBinWidth(1)/
                              htmp->GetBinWidth(i));
      hdisp->SetBinError  (i, htmp->GetBinError(i) *
                              htmp->GetBinWidth(1)/
                              htmp->GetBinWidth(i));
    }


    char functionstring[1000];
    if(snprintf(functionstring, 999, "%f*("
       "%f*(%f*("
         "%f*[2]*(1-[9])/0.257233*exp(-x/0.257233)+"  // H Li-9
         "%f*[3]*(1-[9])/0.171825*exp(-x/0.171825)+"  // H He-8
         "%f*[4]*(1-[5])/6.020366*exp(-x/6.020366)+"  // H N-17
         "%f*[6]*(1-[5])/1.077693*exp(-x/1.077693)+"  // H C-16
         "%f*[7]*(1-[9])/0.025002*exp(-x/0.025002)+"  // H B-13
         "%f*[8]*(1-[9])/0.012624*exp(-x/0.012624))+" // H Li-11
         "[0]*(1-[1])) + "                            // H bg
       "%f*(%f*("
         "%f*[2]*[9]/0.257233*exp(-x/0.257233)+" // Gd Li-9
         "%f*[3]*[9]/0.171825*exp(-x/0.171825)+" // Gd He-8
         "%f*[4]*[5]/6.020366*exp(-x/6.020366)+" // Gd N-17
         "%f*[6]*[5]/1.077693*exp(-x/1.077693)+" // Gd C-16
         "%f*[7]*[9]/0.025002*exp(-x/0.025002)+" // Gd B-13
         "%f*[8]*[9]/0.012624*exp(-x/0.012624))+"// Gd Li-11
         "[0]*[1]))",                            // Gd bg
      denominator*(high-low)/nbin,
      Heff, livefrac, li9ebn, he8ebn, n17ebn, c16ebn, b13ebn, li11ebn,
      Geff, livefrac, li9ebn, he8ebn, n17ebn, c16ebn, b13ebn, li11ebn)
      >= 999){
      printf("Could not create function string\n");
      exit(1);
    };

    delete htmp;

    hdisp->GetYaxis()->SetRangeUser(0,
      (2*sqrt(hdisp->GetMaximum())+hdisp->GetMaximum())*2.5);
    hdisp->Draw("e");
    hdisp->SavePrimitive(cout);

    double bestpar[npar], bestchi2 = 0;
    for(int j = 0; j < npar; j++) bestpar[j] = parsave[j].val;
    int vnpar = npar;
    fcn(vnpar, NULL, bestchi2, bestpar, 0);

    const double chi2tolerance = 0.01;
    fcnstopat = bestchi2+1+chi2tolerance*2;

    vector<TF1 *> bounds;
    vector<int> resulttype;
    printf("Forming plot confidence band...\n");

    bounds.push_back(make_bounds_tf1(0, whichh, runperiod, 0, bestpar, functionstring, low, high2));
    resulttype.push_back(0);

    const int xpoints = 100;
    const double llow = low?low:0.001;
#ifdef BAND
    for(int t = 0; t < npar*2; t++){
      const bool lo = t%2;
      const int par = t/2;
      if(lo) printf("%d min/max\n", par);
      TMinuit * mn = make_a_tminuit();
      fixat(mn, par+1,
            lo?parsave[par].min():parsave[par].max());
      
      for(int y = 0; y < 2; y++) mn->Command("MIGRAD"); 
      const double dchi2 = mn->fAmin-bestchi2; 
      if(dchi2 > 1 + chi2tolerance || dchi2 < -1e-5){ 
        printf("MIGRAD gave a delta of %g :-(\n", dchi2); 
      }else{ 
        printf("OK!\n"); 
        double eedisppar[npar];
        for(int V = 0; V < npar; V++) eedisppar[V] = getpar(mn, V); 
        bounds.push_back(make_bounds_tf1(1, whichh, runperiod, t, eedisppar, functionstring, low, high2)); 
        resulttype.push_back(1);
        bounds[bounds.size()-1]->SetLineColor(kRed);
        bounds[bounds.size()-1]->Draw("same");
        c1->Modified(); c1->Update();
      }
      delete mn;
    }

    for(int t = 0; t < npar*(npar-1)/2; t++){
      int npar1 = 1, npar2 = 2;
      for(int moose = 0; moose < t; moose++){
        npar2++;
        if(npar2 > npar){
          npar1++;
          npar2 = npar1+1;
        }
      }
      if(parsave[npar1-1].fix || parsave[npar2-1].fix) continue;

      TMinuit * mn = make_a_tminuit();

      for(int X = 0; X < npar; X++){
        mn->Command(Form("SET PAR %d %g",X+1,parsave[X].val));
        if(parsave[X].fix) mn->Command(Form("FIX %d", X+1));
      }

      printf("%d %d contour\n", npar1, npar2);
      mn->Command("MIGRAD"); mn->Command("MIGRAD");
      mn->fGraphicsMode = true;
      mn->Command("Set print 0");
      mn->Command(Form("MNCONT %d %d 20", npar1, npar2));
      TGraph * gr = mn->GetPlot()?(TGraph*)((TGraph*)mn->GetPlot())->Clone():NULL;
      if(!gr){
        printf("Could not get contour!\n"); continue;
      }
      mn->Command("Set print -1");
      for(int squirrel = 0; squirrel < gr->GetN(); squirrel++){
        printf("  contour point %d: %f %f\n", squirrel, mn->fXpt[squirrel], mn->fYpt[squirrel]);
        fixat(mn, npar1, gr->GetX()[squirrel]);
        fixat(mn, npar2, gr->GetY()[squirrel]);

        for(int y = 0; y < 2; y++) mn->Command("MIGRAD"); 
        const double dchi2 = mn->fAmin-bestchi2; 
        if(dchi2 > 1 + chi2tolerance || dchi2 < -1e-5){ 
          printf("MIGRAD gave a delta of %g :-(\n", dchi2); 
        }else{ 
          printf("OK!\n"); 
          double eedisppar[npar];
          for(int V = 0; V < npar; V++) eedisppar[V] = getpar(mn, V); 
          bounds.push_back(make_bounds_tf1(2, whichh, runperiod, t, eedisppar, functionstring, low, high2)); 
          resulttype.push_back(2);
          bounds[bounds.size()-1]->SetLineColor(kBlack);
          bounds[bounds.size()-1]->Draw("same");
          c1->Modified(); c1->Update();
        }
      }
      delete mn;
      delete gr;
    }

    const unsigned int ncurves = 3000;
    for(unsigned int t = 0; t < ncurves; t++){
      double dchi2 = 0;
      if(t%10 == 0) printf("Random %d/%d\n", t, ncurves);
      double eedisppar[npar];
      do{
        for(int j=0; j<npar; j++) eedisppar[j]=parsave[j].rand();
        double chi2 = 0;
        fcnearlystop = true;
        fcn(vnpar, NULL, chi2, eedisppar, 0);
        fcnearlystop = false;
        dchi2 = chi2-bestchi2;
      }while(fabs(dchi2-1) > chi2tolerance);

      bounds.push_back(make_bounds_tf1(3, whichh, runperiod, t, eedisppar, functionstring, low, high2));
      resulttype.push_back(3);
      bounds[bounds.size()-1]->SetLineColor(kGreen+2);
      bounds[bounds.size()-1]->Draw("same");
      c1->Modified(); c1->Update();
    }

    printf("Making high and low graphs\n");
    TGraph ghigh, glow;
    for(int ix = -1; ix <= xpoints; ix++){
      const double x = ix == -1?0:
        exp(log(llow) + double(ix)*(log(high2)-log(llow))/xpoints);
      double lowest = 1e50, highest = -1e50;
      int besttl = -1, bestth = -1;
      for(unsigned int t = 0; t < bounds.size(); t++){
        const double y = bounds[t]->Eval(x);
        if(y < lowest) besttl=t, lowest = y;
        if(y > highest)bestth=t, highest= y;
      }

      const double dispx = x < high?x:
        high + (x - high)*(double(nbin2)/(high2-high))/
        (double(nbin )/(high -low ));
      ghigh.SetPoint(ghigh.GetN(), dispx, highest);
      glow .SetPoint( glow.GetN(), dispx,  lowest);
    }

#endif
    TGraph * gbest = new TGraph();
    for(int ix = -1; ix <= xpoints; ix++){
      const double x = ix == -1?0:
        exp(log(llow) + double(ix)*(log(high2)-log(llow))/xpoints);
      const double dispx = x < high?x:
        high + (x - high)*(double(nbin2)/(high2-high))/
        (double(nbin )/(high -low ));

      gbest->SetPoint(gbest->GetN(), dispx, bounds[0]->Eval(x));
    }
    gbest->SetLineColor(kBlack);
    gbest->SetLineWidth(3);
    gbest->SetLineStyle(kDashed);
#ifdef BAND

    printf("Making combined graph\n");
    TGraph * gall = new TGraph();
    for(int p = glow.GetN()-1; p >= 0; p--)
      gall->SetPoint(gall->GetN(),  glow.GetX()[p],  glow.GetY()[p]);
    for(int p = 0; p < ghigh.GetN(); p++)
      gall->SetPoint(gall->GetN(), ghigh.GetX()[p], ghigh.GetY()[p]);

    gall->SetFillColor(TColor::GetColor("#ccccff"));
    gall->SetFillStyle(1001);
    gall->SetNameTitle("errorband", "errorband");
    gall->SavePrimitive(cout);
    gall->Draw("f");
#endif
    gbest->SetNameTitle("bestfit", "bestfit");
    gbest->SavePrimitive(cout);
    gbest->Draw("l");
#ifdef BAND

    for(unsigned int soup = 0; soup < bounds.size(); soup++){
      if(resulttype[soup] == 1 || resulttype[soup] == 2) {
        bounds[soup]->SetLineWidth(1);
        bounds[soup]->SetNpx(400);
        bounds[soup]->Draw("same");
      }
    }

    hdisp->GetXaxis()->SetTickLength(0);
    hdisp->GetXaxis()->SetLabelSize(0);
    TGaxis * na1 = new TGaxis(low, 0, high, 0, low, high, 508, "+");
    TGaxis * na2 = new TGaxis(high, 0,
        hdisp->GetBinLowEdge(hdisp->GetNbinsX()+1),0,high,high2,503,"+");
    na1->SetLabelFont(42);
    na2->SetLabelFont(42);
    na1->SetLabelSize(0.05);
    na2->SetLabelSize(0.05);
    na1->Draw();
    na2->Draw();
#endif
    hdisp->Draw("samee");
  }
}

int whichc = -1;
TCanvas * cans[100]; // oh no
void contour(TMinuit * mn, const int par1, const int par2,
             const double xrange, const double yrange, const int points,
             const char * const comment,
             const char * const nintyname = "ninty",
             const char * const sigmaname = "sigma")
{
  mn->fGraphicsMode = true;
  whichc++;
  TCanvas * c = new TCanvas(Form("c%d_%d%d", whichc, par1, par2),
                            Form("c%d_%d%d", whichc, par1, par2),
                            600, 350);
  mn->Command("MIGRAD");
  const double minx = getpar(mn, par1-1);
  const double miny = getpar(mn, par2-1);

  const int oldprintlevel = mn->fISW[4];
  mn->Command("Set print 0");

  mn->fUp = 2.30258509299404590; // 68% in 2D
  mn->Command(Form("mncont %d %d %d", par1, par2, points));
  TGraph * sigma_2d =
    mn->GetPlot()?(TGraph*)((TGraph*)mn->GetPlot())->Clone():NULL;

  mn->Command("MIGRAD");
  mn->fUp = 4.61; // 90% in 2D
  mn->Command(Form("mncont %d %d %d", par1, par2, points));
  TGraph * ninty_2d =
    mn->GetPlot()?(TGraph*)((TGraph*)mn->GetPlot())->Clone():NULL;

  mn->Command("MIGRAD");
  mn->fUp = 2.71; // 90% in 1D
  mn->Command(Form("mncont %d %d %d", par1, par2, points));
  TGraph * ninty_1d =
    mn->GetPlot()?(TGraph*)((TGraph*)mn->GetPlot())->Clone():NULL;


  if(ninty_2d){
    ninty_2d->SetNameTitle(nintyname, nintyname);
    ninty_2d->SavePrimitive(cout);
    ninty_2d->SetFillColor(kViolet);
    ninty_2d->Draw("alf");
    ninty_2d->GetXaxis()->SetRangeUser(0, xrange);
    ninty_2d->GetYaxis()->SetRangeUser(0, yrange);
    ninty_2d->GetXaxis()->SetTitle(mn->fCpnam[par1-1]);
    ninty_2d->GetYaxis()->SetTitle(mn->fCpnam[par2-1]);
    ((TGaxis*)(ninty_2d->GetXaxis()))->SetMaxDigits(3);
    ((TGaxis*)(ninty_2d->GetYaxis()))->SetMaxDigits(3);
  }
  if(ninty_1d) ninty_1d->SetLineColor(kRed),   ninty_1d->Draw("l");
  else printf("ACK! Couldn't make 90%% 1D contour!\n");
  if(sigma_2d){
    sigma_2d->SetNameTitle(sigmaname, sigmaname);
    sigma_2d->SavePrimitive(cout);
    sigma_2d->SetLineColor(kBlack);
    sigma_2d->Draw("l");
  }
/*  mn->fUp = 11.83; // 99.73% in 2D
  mn->Command(Form("mncont %d %d %d", par1, par2, points));
  TGraph * threesigma_2d =
    mn->GetPlot()?(TGraph*)((TGraph*)mn->GetPlot())->Clone():NULL;

  if(threesigma_2d) threesigma_2d->SetLineStyle(kDashed),
                    threesigma_2d->Draw("l"); */

  TMarker * best = new TMarker(minx, miny, kStar);
  best->Draw();

  if(!strcmp(comment, "")){
    TLatex * t = new TLatex(0.5, 0.8, comment);
    t->SetTextSize(0.07);
    t->SetNDC(1);
    t->Draw();
  }

  cans[whichc] = c;
  for(int i = 0; i <= whichc; i++){
    cans[i]->Modified();
    cans[i]->Update();
  }

  mn->Command(Form("set print %d", oldprintlevel));
  mn->fGraphicsMode = false;
}


double lratsig(const double l1, const double l2)
{
  const double dll = (l1-l2)/2;
  const double rat = exp(dll);
  return TMath::NormQuantile(1-1/(2*rat));
}

float getredchi2(const int run, const int trig)
{
  static TTree * fidot = NULL;
  static TFile * fidof = NULL;
  static int lastrun = 0;
  static TBranch * br[3];
  static float ids_chi2;
  static int nidtubes, nivtubes;
  TDirectory * old = gDirectory;
  if(lastrun != run){
    if(fidot) delete fidot;
    if(fidof) delete fidof;
    fidof = new
      TFile(Form("/cp/s4/strait/fido_seq010/fido.%07d.root",run),"read");
    if(!fidof || fidof->IsZombie()) {
      lastrun = 0, fidot = NULL, fidof = NULL;
      return 0;
    }
    fidot = (TTree *) fidof->Get("RecoMuonFIDOInfoTree"); 
    if(!fidot){
      lastrun = 0, fidot = NULL, fidof = NULL;
      return 0;
    }
    fidot->SetMakeClass(1);
    fidot->SetBranchAddress("ids_chi2", &ids_chi2);
    fidot->SetBranchAddress("nidtubes", &nidtubes);
    fidot->SetBranchAddress("nivtubes", &nivtubes);
    fprintf(stderr, ".");
    br[0] = fidot->GetBranch("ids_chi2");
    br[1] = fidot->GetBranch("nidtubes");
    br[2] = fidot->GetBranch("nivtubes");
  }

  for(int i = 0; i < 3; i++) br[i]->GetEntry(trig);
  fidot->GetEntry(trig);

  lastrun = run;

  old->cd();

  if(nidtubes + nivtubes - 6 < 0) return 0;
  return ids_chi2/(nidtubes + nivtubes - 6);
}

double getprimaryresult()
{
  ifstream infile("li9_finalfit_-1.out.h");
  if(!infile.is_open()){
    printf("couldn't open li9_finalfit_-1.out.h\n");
    exit(1);
  }

  string scratch;
  double primaryresult;
  for(int i = 0; i < 11; i++){
    if(!(infile >> scratch)){
      printf("li9_finalfit_-1.out.h ended unexpectedly\n");
      exit(1);
    }
  }
  if(!(infile >> primaryresult)){
    printf("li9_finalfit_-1.out.h ended unexpectedly or had bad format\n");
    exit(1);
  }
   
  return primaryresult;
}

void li9_finalfit(int neutrons = -1, int contourmask = 0)
{
  if(neutrons < -1){
    puts("You gave me less than -1 neutrons, so I am just compiling");
    return;
  }

  // First factor takes into account the efficiency of selecting a
  // neutron after a muon, second is the prompt energy cut, third as
  // documented above
  // 
  // Neutron efficiencies are assuming any within the Michel window
  // are rejected.
  Geff=(neutrons > 0?pow(neff_dr_1000_gd*neff_dt_targ, neutrons):1)
    *0.996*Geff_sans_prompt_or_mun,
  Heff=(neutrons > 0?pow(neff_dr_1000_h*neff_dt_h, neutrons):1)
    *0.993*Heff_sans_prompt_or_mun;

  ///////////////////////////////////////////////////////////////////
 
  char nodistcut[1000];
  snprintf(nodistcut, 999,
#ifdef HP
  "mx**2+my**2 < 1050**2 && mz > -1175 && "
  "abs(dedxslant - 8847.2) < 1000 && " // cut on reduced chi2 below
                                       // since it isn't in the files
#endif
    "%smiche < 12 && !earlymich && prompttime > 100e9 && dt < 100e3",
           neutrons >= 0?Form("nlate==%d&&", neutrons):"");
  char cut[1000];
  snprintf(cut, 999, "%s && dist < %f %s %s %s", nodistcut, dist,
      rrmlivetimes[0]==0?"&& reactorpowerbin(run) != 0":"",
      rrmlivetimes[1]==0?"&& reactorpowerbin(run) != 1":"",
      rrmlivetimes[2]==0?"&& reactorpowerbin(run) != 2":"");

  printf("Cut is: %s\n", cut);
  ////////////////////////////////////////////////////////////////////

  TTree tg("t", "t");
  TTree th("t", "t");
  tg.ReadFile("/cp/s4/strait/li9ntuples/li9-20150219.Gd.ntuple");
  th.ReadFile("/cp/s4/strait/li9ntuples/li9-20150219.H.ntuple");

  TTree * tgsel = tg.CopyTree(cut);
  TTree * thsel = th.CopyTree(cut);
  if(!tgsel || !thsel){
    fprintf(stderr, "Failed to make cut trees\n");
    return;
  }

/*
  const double plotsigtime = 0.4;

  TCanvas * dxy = new TCanvas("dxy", "dxy", 200, 200);
  th.Draw("my-dy:mx-dx", Form("%s && dt/1000 > 5  ", cut), ".");
  tg.Draw("my-dy:mx-dx", Form("%s && dt/1000 > 5  ", cut), ".same");
  th.Draw("my-dy:mx-dx", Form("%s && dt/1000 < %f", cut, plotsigtime),"%same");
  lastgraphstyle(dxy, 0);
  tg.Draw("my-dy:mx-dx", Form("%s && dt/1000 < %f", cut, plotsigtime),"%same");
  lastgraphstyle(dxy, 0);

  TCanvas * dxz = new TCanvas("dxz", "dxz", 200, 200);
  th.Draw("mz-dz:mx-dx", Form("%s && dt/1000 > 5  ", cut), ".");
  tg.Draw("mz-dz:mx-dx", Form("%s && dt/1000 > 5  ", cut), ".same");
  th.Draw("mz-dz:mx-dx", Form("%s && dt/1000 < %f", cut, plotsigtime),"%same");
  lastgraphstyle(dxz, 0);
  tg.Draw("mz-dz:mx-dx", Form("%s && dt/1000 < %f", cut, plotsigtime),"%same");
  lastgraphstyle(dxz, 0);

  TCanvas * dyz = new TCanvas("dyz", "dyz", 200, 200);
  th.Draw("mz-dz:my-dy", Form("%s && dt/1000 > 5  ", cut), ".");
  tg.Draw("mz-dz:my-dy", Form("%s && dt/1000 > 5  ", cut), ".same");
  th.Draw("mz-dz:my-dy", Form("%s && dt/1000 < %f", cut, plotsigtime),"%same");
  lastgraphstyle(dyz, 0);
  tg.Draw("mz-dz:my-dy", Form("%s && dt/1000 < %f", cut, plotsigtime),"%same");
  lastgraphstyle(dyz, 0);

  TCanvas * muenergy = new TCanvas("muenergy", "muenergy", 200, 200);
  th.Draw("dedxslant >> muqivbg(25, 0, 12000)", Form("%s && dt/1000 > 5  ", cut), "hist");
  tg.Draw("dedxslant >> +muqivbg", Form("%s && dt/1000 > 5  ", cut), "esame");
  th.Draw("dedxslant >> muqivsig(25, 0, 12000)", Form("%s && dt/1000 < 0.4", cut),"same");
  tg.Draw("dedxslant >> +muqivsig", Form("%s && dt/1000 < 0.4", cut),"same");
*/

  // yes, all of them are floats
  float tim, run;
  #ifdef HP
    float mutrig;
  #endif

  TTree * tsels[2] = { tgsel, thsel };

  for(int tree = 0; tree < 2; tree++){
    printf("Reading cut tree %d with %lld entries\n",
           tree, tsels[tree]->GetEntries());
    tsels[tree]->SetBranchAddress("dt", &tim);
    tsels[tree]->SetBranchAddress("run", &run);
#ifdef HP
    tsels[tree]->SetBranchAddress("mutrig", &mutrig);
#endif
    for(int i = 0; i < tsels[tree]->GetEntries(); i++){
      tsels[tree]->GetEntry(i);
#ifdef HP
      if(getredchi2(int(run), int(mutrig)) >= 2){
        printf("dropping bad mu fit\n");
        continue;
      }
#endif
      const int rpb = reactorpowerbin(int(run));
      events.push_back(ev(tree == 1, (tim - 0.002)/1000, rpb));
    }
  }

  //////////////////

  TMinuit * mn = make_a_tminuit();

  vector< vector<parres> > parsaves;
  {
    const unsigned int nfix = 8;
    const int fix[nfix] = { 3, 4, 5, 6, 7, 8, 9, 10 };
    for(unsigned int i = 0; i < nfix; i++) fixatzero(mn, fix[i]);
    mn->Command("MIGRAD");
    printf("%sNo Li-9 or other bn isotopes (%.2f)%s\n",
           RED, mn->fAmin, CLR);
  }
  const double chi2nothing = mn->fAmin;

  {
    setupmn(mn, expectedgdfrac);
    dopull = false;
    const unsigned int nfix = 2;
    const int fix[nfix] = { 3, 4 };
    for(unsigned int i = 0; i < nfix; i++) fixatzero(mn, fix[i]);
    if(neutrons > 0){ fixatzero(mn, 8); fixatzero(mn, 9); }
    mn->Command("MIGRAD");
  }
  const double chi2all_exceptli9he8_nopull = mn->fAmin;

  {
    printf("%sRead the no-pull values of N-17 and C-16 here:%s\n", RED, CLR);
    setupmn(mn, expectedgdfrac);
    dopull = false;
    if(neutrons > 0){ fixatzero(mn, 8); fixatzero(mn, 9); }
    mn->Command("SET LIM 7 0 10");
    mn->Command("MIGRAD");
    mn->Command("MIGRAD");
    mn->Command("MINOS 10000 5");
    mn->Command("MINOS 10000 7");
    mn->Command("SHOW MIN");
  }
  const double chi2all_nopull = mn->fAmin;

  printf("%s%sSignificance of any betan (%f vs. %f), no cheaty pulls: %.1f%s\n",
         RED,
         neutrons == -1?"TECHNOTE results.tex betansignificance: ":
         neutrons ==  1?"TECHNOTE 8.4.1: ":"",
         chi2nothing, chi2all_nopull, 
         lratsig(chi2nothing, chi2all_nopull), CLR);

  printf("%sSignificance of li9/he8 over other bn&accidental without pull: %.2f%s\n",
         RED, lratsig(chi2all_exceptli9he8_nopull, chi2all_nopull), CLR);

  {
    setupmn(mn, expectedgdfrac);
    const unsigned int nfix = 2;
    const int fix[nfix] = { 3, 4 };
    for(unsigned int i = 0; i < nfix; i++) fixatzero(mn, fix[i]);
    if(neutrons > 0){ fixatzero(mn, 8); fixatzero(mn, 9); }
    mn->Command("MIGRAD");
  }
  const double chi2all_exceptli9he8_withpull = mn->fAmin;

  {
    setupmn(mn, expectedgdfrac);
    dopull = true;
    if(neutrons > 0){ fixatzero(mn, 8); fixatzero(mn, 9); }
    mn->Command("MIGRAD");
  }
  const double chi2all_withpull = mn->fAmin;

  printf("%s%sSignificance of li9/he8 over other bn&accidental with pull: %.1f%s\n",
         RED,
         neutrons == -1?"TECHNOTE results.tex linineheeightsignificance: ":
         neutrons ==  1?"TECHNOTE 8.4.1: ":"",
         lratsig(chi2all_exceptli9he8_withpull, chi2all_withpull), CLR);

  {
    setupmn(mn, expectedgdfrac);
    const unsigned int nfix = 6;
    const int fix[nfix] = { 4, 5, 6, 7, 8, 9 };
    for(unsigned int i = 0; i < nfix; i++) fixatzero(mn, fix[i]);
    mn->Command("MIGRAD");
    mn->Command("MINOS 10000 3");
    mn->Command("SHOW min");
    printf("%sLi-9 prob without other isotopes (%.2f): %g %g +%g%s\n",
           RED, mn->fAmin, getpar(mn, 2), mn->fErn[2], mn->fErp[2], CLR);
  }
  //const double chi2justli9 = mn->fAmin;

  {
    setupmn(mn, expectedgdfrac);
    const unsigned int nfix = 5;
    const int fix[nfix] = { 4, 6, 5, 8, 9 };
    for(unsigned int i = 0; i < nfix; i++) fixatzero(mn, fix[i]);
    mn->Command("MIGRAD");
    mn->Command("MINOS 10000 3");
    mn->Command("SHOW min");
    printf("%sLi-9 with C-16 only (%.2f): %f %f +%f%s\n", RED,
      mn->fAmin, getpar(mn, 2), mn->fErn[2], mn->fErp[2], CLR);
  }
  //const double chi2_li9_c16 = mn->fAmin;

  {
    setupmn(mn, expectedgdfrac);
    const unsigned int nfix = 5;
    const int fix[nfix] = { 4, 6, 7, 8, 9 };
    for(unsigned int i = 0; i < nfix; i++) fixatzero(mn, fix[i]);
    mn->Command("MIGRAD");
    mn->Command("MINOS 10000 3");
    mn->Command("SHOW min");
    printf("%sLi-9 with N-17 only (%.2f): %f %f +%f%s\n", RED,
      mn->fAmin, getpar(mn, 2), mn->fErn[2], mn->fErp[2], CLR);
  }
  //const double chi2_li9_n17 = mn->fAmin;

  {
    setupmn(mn, expectedgdfrac);
    const unsigned int nfix = 5;
    const int fix[nfix] = { 5, 6, 7, 8, 9 };
    for(unsigned int i = 0; i < nfix; i++) fixatzero(mn, fix[i]);
    mn->Command("MIGRAD");
    mn->Command("MINOS 10000 3");
    mn->Command("SHOW min");
    printf("%sLi-9 prob with only He-8 (%.2f): %f %f +%f%s\n", RED,
      mn->fAmin, getpar(mn, 2), mn->fErn[2], mn->fErp[2], CLR);
  }
  //const double chi2_li9_he8 = mn->fAmin;

  {
    setupmn(mn, expectedgdfrac);
    const unsigned int nfix = 1;
    const int fix[nfix] = { 4 };
    for(unsigned int i = 0; i < nfix; i++) fixatzero(mn, fix[i]);
    if(neutrons > 0){
      const unsigned int nfixn = 2;
      const int fixn[nfixn] = { 8, 9 };
      for(unsigned int i = 0; i < nfixn; i++) fixatzero(mn, fixn[i]);
    }
    mn->Command("MIGRAD");
    mn->Command("MIGRAD");
    mn->Command("MINOS 10000 1 3 5 11 12");
    mn->Command("SHOW min");
    printf("%sLi-9 prob/C-12 capture with everything except "
           "He-8 (%.2f): %g %g +%g%s\n",
           RED, mn->fAmin, getpar(mn, 2), mn->fErn[2], mn->fErp[2], CLR);
    printf("%s%sLi-9 prob/C-13 capture with everything except "
           "He-8 (%.2f): (%.1f %.1f +%.1f)%%%s\n",
           RED,
           neutrons == -1?"TECHNOTE results.tex probNineLiFromThirteenC: ":"",
           mn->fAmin,
           getpar(mn, 2)*Nc12cap/Nc13cap*100,
           mn->fErn[2]*Nc12cap/Nc13cap*100,
           mn->fErp[2]*Nc12cap/Nc13cap*100,
           CLR);
    printf("%s%sLi-9 prob/C-nat capture w/ everything except "
           "He-8 (%.2f): (%.2f %.2f +%.2f)*10**-4%s\n",
           RED,
           neutrons == -1? "TECHNOTE results.tex primaryresult: ":
           neutrons ==  1? "TECHNOTE results.tex probNineLiFromTwelveC: ":"",
           mn->fAmin,
           getpar(mn, 2)*Nc12cap/(Nc12cap+Nc13cap)*10000,
           mn->fErn[2]*Nc12cap/(Nc12cap+Nc13cap)*10000,
           mn->fErp[2]*Nc12cap/(Nc12cap+Nc13cap)*10000,
           CLR);
    if(neutrons == -1)
      // NOTE-driftcheck: horrible and fragile passing of this result to
      // the 1-neutron code. Note the space after the number, which is
      // necessary!
      printf("const double primaryresult = %f ;\n",
             getpar(mn, 2)*Nc12cap/(Nc12cap+Nc13cap));
  }
  //const double chi2_allbut_he8 = mn->fAmin;

  {
    setupmn(mn, expectedgdfrac);
    if(neutrons > 0){
      const unsigned int nfixn = 2;
      const int fixn[nfixn] = { 8, 9 };
      for(unsigned int i = 0; i < nfixn; i++) fixatzero(mn, fixn[i]);
    }
    mn->Command("MIGRAD");
    mn->Command("MINOS 10000");
    mn->Command("SHOW min");
    printf("%sLi-9 prob with all other isotopes (%.2f): %f %f +%f%s\n",
      RED, mn->fAmin, getpar(mn, 2), mn->fErn[2], mn->fErp[2], CLR);
    vector<parres> parsave;
    for(int i = 0; i < npar; i++)
      parsave.push_back(parres(getpar(mn, i), mn->fErp[i], mn->fErn[i],
                        getlimup(mn, i), getlimlo(mn, i),
      /* XXX fragile */ neutrons > 0 && (i == 7 || i == 8)));
    parsaves.push_back(parsave);
  }
  const double chi2_all = mn->fAmin;

  const double maxb13 = getpar(mn, 7) + mn->fErp[7] * sqrt(2.30258509299404590);
  printf("%sApprox 90%% upper limit B-13: %f\n%s", RED, maxb13, CLR);
  printf("%sBut look at the contour!\n%s", RED, CLR);

  /* // These are wrong if we are using any pulls
    printf("%s", RED);
    printf("Li-9 preferred over nothing by %f\n",
           chi2nothing - chi2justli9 < 0?0:lratsig(chi2nothing, chi2justli9));

    printf("Li-9 + C-16 preferred over just Li-9 by %f\n",
           chi2justli9 - chi2_li9_c16<0?0:lratsig(chi2justli9, chi2_li9_c16));

    printf("Li-9 + N-17 preferred over just Li-9 by %f\n",
           chi2justli9 - chi2_li9_n17<0?0:lratsig(chi2justli9, chi2_li9_n17));

    printf("All preferred over Li-9 + N-17 by %f\n",
           chi2_li9_n17 - chi2_all<0?0:lratsig(chi2_li9_n17, chi2_all));
    printf("All preferred over nothing by %f\n",
           chi2nothing - chi2_all<0?0:lratsig(chi2nothing, chi2_all));
    printf("%s", CLR);
  }
  */

  const int npoint = 100;

  // Li-9 vs. N-17 with nothing else
  if(contourmask & 0x01){
    setupmn(mn, expectedgdfrac);
    fixatzero(mn, 4);
    fixatzero(mn, 7);
    fixatzero(mn, 8);
    fixatzero(mn, 9);
    contour(mn, 3, 5,  0.0079, 2, npoint, "Nothing else");
  }

  // Li-9 vs. C-16 with no He-8
  if(contourmask & 0x02){
    setupmn(mn, expectedgdfrac);
    fixatzero(mn, 4);
    contour(mn, 3, 7, 0.0079, 2, npoint, "No ^{8}He");
  }

  // Li-9 vs. B-13 with no Li-11 or He-8
  if(contourmask & 0x04){
    setupmn(mn, expectedgdfrac);
    fixatzero(mn, 4);
    fixatzero(mn, 9);
    contour(mn, 3, 8,  0.0079, 2, npoint, "No ^{8}He or ^{11}Li");
  }

  // Li-9 vs. He-8
  if(contourmask & 0x08){
    setupmn(mn, expectedgdfrac);
    if(neutrons == -1) mn->Command("SET LIM 3 -3e-5 1e-3");
    contour(mn, 3, 4, 0.0079, 0.00149, npoint, "",
      neutrons == -1?"ignoring90":neutrons==0?"without90":"with90",
      neutrons == -1?"ignoring1s":neutrons==0?"without1s":"with1s");
  }

  // Li-9 vs. He-8 with nothing else
  if(contourmask & 0x10){
    setupmn(mn, expectedgdfrac);
    fixatzero(mn, 5);
    fixatzero(mn, 6);
    fixatzero(mn, 7);
    fixatzero(mn, 8);
    fixatzero(mn, 9);
    contour(mn, 3, 4, 0.0079, 0.00149, npoint, "Nothing else");
  }

  // B-13 vs. Li-11
  if(contourmask & 0x20){
    setupmn(mn, expectedgdfrac);
    mn->Command("SET LIM 8 0 3");
    contour(mn, 8, 9, 3, 0.008, npoint, "");
  }

  if(neutrons == 1){
    setupmn(mn, expectedgdfrac);
    const double primaryresult = getprimaryresult(); // NOTE-driftcheck
    fixat(mn, 3, primaryresult);
    mn->Command("MIGRAD");
    printf("%sTECHNOTE 8.4.1: Sigmas between here & Li-9 prob = %f: %.1f%s\n",
           RED, primaryresult, lratsig(mn->fAmin, chi2_all), CLR);
  }
  
  //////////////////////////////////////////////////////////////////////
  //drawhist(tgsel, thsel, parsaves[0], 48, 1, 97);
  if(neutrons == -1)  drawhist(tgsel, thsel, parsaves[0], 15, 0,  3,10, 100);
  else if(neutrons==1)drawhist(tgsel, thsel, parsaves[0],  3, 0,  3, 2, 100);
  else                drawhist(tgsel, thsel, parsaves[0], 12, 0,0.6, 1, 100);


  /* setupmn(mn, expectedgdfrac);
  fixatzero(mn, 4);
  const bool withall = true;
  if(!withall){
    fixatzero(mn, 5);
    fixatzero(mn, 7);
    fixatzero(mn, 8);
    fixatzero(mn, 9);
  }
  mn->SetPrintLevel(-1);
  for(int i = 0; i <= 150; i++){
    const double val = double(i)*4e-6;
    mn->Command("REL 3");
    mn->Command(Form("SET PAR 3 %E", val));
    mn->Command("FIX 3");
    mn->Command("MIGRAD");
    printf("%s%s->SetPoint(%d, %f, %f);\n",
           neutrons == -1?"ignoring":neutrons == 0?"withoutn":"withn",
           withall?"all":"",
           i, val*1e3, mn->fAmin-chi2_all ); 
  }
*/

/*
  TCanvas * xy = new TCanvas("xy", "xy", 200, 200);
  th.Draw("dy/1000:dx/1000", Form("%s && dt/1000 > 10 && dt/1000 < 30", cut), ".");
  tg.Draw("dy/1000:dx/1000", Form("%s && dt/1000 > 10 && dt/1000 < 30", cut), ".same");
  th.Draw("dy/1000:dx/1000", Form("%s && dt/1000 < 0.4", cut),"%same");
  lastgraphstyle(xy, 0);
  tg.Draw("dy/1000:dx/1000", Form("%s && dt/1000 < 0.4", cut),"%same");
  lastgraphstyle(xy, 0);

  TCanvas * rz = new TCanvas("rz", "rz", 200, 200);
  th.Draw("dz:(dx**2+dy**2)", Form("%s && dt/1000 > 10 && dt/1000 < 30", cut), ".");
  tg.Draw("dz:(dx**2+dy**2)", Form("%s && dt/1000 > 10 && dt/1000 < 30", cut), ".same");
  th.Draw("dz:(dx**2+dy**2)", Form("%s && dt/1000 < 0.4", cut),"%same");
  lastgraphstyle(rz, 0);
  tg.Draw("dz:(dx**2+dy**2)", Form("%s && dt/1000 < 0.4", cut),"%same");
  lastgraphstyle(rz, 0);

  TCanvas * dr = new TCanvas("dr", "dr", 200, 200);
  th.Draw("(dist**3)/1e9 >>  distbg(25,0,0.5)",Form("%s && dt/1000 > 5",  nodistcut), "hist");
  th.Draw("(dist**3)/1e9 >> distsig(25,0,0.5)",Form("%s && dt/1000 < 0.4",nodistcut),"esame");
  tg.Draw("(dist**3)/1e9 >>  +distbg",           Form("%s && dt/1000 > 5",  nodistcut), "histsame");
  tg.Draw("(dist**3)/1e9 >> +distsig",           Form("%s && dt/1000 < 0.4",nodistcut),"esame");
  distsig=distsig;
  distbg=distbg;
  distbg->Sumw2();
  distbg->Scale(distsig->Integral(5, 100)/distbg->Integral(5, 100));
  distsig->SetLineColor(kRed);
  distsig->SetMarkerColor(kRed);
  distbg->GetYaxis()->SetRangeUser(0, 25);
*/

/*
  TCanvas * rzh = new TCanvas("rzh", "rzh", 200, 200);
  th.Draw("dz >> dzall(10, -1800, 1800)", cut, "hist");
  th.Draw("dz >> dzn16", Form("%s && dt/1000 > 0.75 && dt/1000 < 5", cut),"samee");
  tg.Draw("dz >> +dzall", cut, "hist");
  tg.Draw("dz >> +dzn16", Form("%s && dt/1000 > 0.75 && dt/1000 < 5", cut),"samee");
  dzall=dzall;
  dzn16=dzn16;
  dzn16->Scale(dzall->Integral()/dzn16->Integral());
*/

/*
  printf("Candidates:\n");
  th.SetScanField(0);
  tg.SetScanField(0);
  th.Scan("run:trig:dt:sqrt(dx*dx+dy*dy):dx:dy:dz", cut);
  tg.Scan("run:trig:dt:sqrt(dx*dx+dy*dy):dx:dy:dz", cut);
*/
}
