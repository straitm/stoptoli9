#include <math.h>

struct ve{
  double val, err;
};

// Disable coloring now that I am marking important output with "TECHNOTE"
// and "^const "
const char * const RED     = ""; // bold red
const char * const CLR     = ""; // clear

const char * const rootfile3up_extended = // XXX now out of date
  "/cp/s4/strait/fullfido-300s-3-25MeV-20150924-3rdpub+post3rdpub.root";
const char * const rootfile3up =
  "/cp/s4/strait/fullfido-300s-3-25MeV-20151023.root";
const char * const rootfile0up =
  "/cp/s4/strait/fullfido-100s-0-25MeV-20151023.root";
const char * const rootfile_be12 =
  "/cp/s4/strait/be12-20151102.root";
const char * const rootfilemuinfo =
  "/cp/s4/strait/fullfido-muinfo-20151105.root";
const char * const rootfilethru = // XXX now out of date -- but doesn't matter, since it is only for energy calibration
  "/cp/s4/strait/fullfido-stopandthru-1ms-3-25MeV-20150325.root";
const char * const rootfilemusic =
  "/cp/s4/strait/DC_FD_Muon_Surface_Dist_Ecut_both.root";

const double mum_contamination = 0.0028,
             mum_contamination_err = 0.0019;

const double mass_n14targ = 4.763769, mass_n14gc = 2.831346;
const double mass_o16targ = 20.865041,
             mass_o16targves = 138.24*85./(85.+58.),
             mass_o16targbits = 138.24*58./(85.+58.),
             mass_o16gc = 3.235824,
             mass_o16gcves_effective = 200.727 * 0.4 * 0.8 * 0.45;

const double gd_fraction = 0.851;

const double mich_eff = 0.9996;

// Light noise efficiency, based on DC3rdPub's 0.0124+-0.0008% for
// delayed coincidences.
const double light_noise_eff = 0.999938;

// I made this up, with some studies to back it.
const double f_neff_dt_error = 0.01;

// Taking the plot from doc-4807, slide 3, corrected by MC correction
// factor from doc-4450, slide 8.
const double neff_dr_800_gd =  0.991*0.993;
const double neff_dr_1000_gd = 0.997*0.993;

// Direct from doc-4450, slide 17.
const double neff_dr_800_h = 0.933297;
const double neff_dr_1000_h = 0.972960;

const double neff_dr_800_targ = neff_dr_800_gd*gd_fraction +
                                neff_dr_800_h*(1-gd_fraction);
const double neff_dr_1000_targ = neff_dr_1000_gd*gd_fraction +
                                neff_dr_1000_h*(1-gd_fraction);

// Using the method of CDF memo 5928, particularly section 3.5
// TF1 hey("hey", "ROOT::Math::inc_gamma(1, x*[0])", 0, 10)
// TGraph g
// a = 6.13/9.64 (lower and mean o16beta rates)
// b = 13.15/9.64 (upper and mean o16beta rates)
// for(int i = 0; i < 1000; i++){
//   double val = i/100.;
//   hey.SetParameter(0, val);
//   g.SetPoint(i, val, hey.Integral(a,b)/(b-a));
// }
// Where this crosses 0.9 is the limit for zero signal, zero background
// Divide by the case of perfectly known target mass, 2.3026 to get this
// number.
const double lim_inflation_for_obeta = 1.055239;

// Systematic error on B-12 gamma lines due to uncertainty in the energy
// correction function.  Technote section 11.7
const double b12lineEsyst[4] = { 0.00170994 / 0.315747,
                                 0.00174438/0.0617006,
                                 0.00107529/0.464471,
                                 0.000441855/0.0258321};

// Systematic error on B-12 gamma lines due to uncertainly on neutron 
// background.  Technote section 11.7
const double b12lineNsyst[4] = { 0,
                                 0.03,
                                 0.01,
                                 0 };

// Given a mu- stop, the fraction that capture atomicly on carbon
const double c_atomic_capture_prob = 0.998;
const double c_atomic_capture_prob_err = 0.001;


// All in milliseconds
const double b12life = 20.20/log(2.);
const double b12life_err = 0.02/log(2.); // 0.1%

const double be12life = 21.49/log(2.);
const double be12life_err = 0.03/log(2.);

const double n16life = 7130./log(2.);
const double n16life_err = 20./log(2.); // 0.3%

const double b13life = 17.33/log(2.);
const double b13life_err = 0.17/log(2.);

const double be11life=13.81  /log(2.);
const double be11life_err= 0.08  /log(2.); // +- 0.08 (0.6%)

const double c15life = 2.449 /log(2.);
const double c15life_err = 0.005 /log(2.); // +- 0.005 (0.2%)

// This is the NNDC value. Could alternatively use 838.75+-0.32 from PRC
// 82, 027309
const double li8life = 839.9/log(2.);
const double li8life_err = 0.9/log(2.); // 0.1%

const double li9life = 178.3/log(2.);
const double li9life_err = 0.4/log(2.);

const double he8life =   119.1/log(2.);
const double he8life_err = 1.2/log(2.);

const double n17life =  4173./log(2.);
const double n17life_err = 4./log(2.);

const double c16life =   747./log(2.);
const double c16life_err = 8./log(2.);

const double li11life =     8.75/log(2.);
const double li11life_err = 0.14/log(2.);

const double c9life =   126.5/log(2.);
const double c9life_err = 0.9/log(2.);

const double n12life =    11.000/log(2.);
const double n12life_err = 0.016/log(2.);

// In milliseconds
const double mulife = 2196.9811e-6;

const double lifetime_c12 = 2028.e-6;
const double lifetime_c12_err = 2.e-6;

const double lifetime_c13 = 2037.e-6;
const double lifetime_c13_err = 8.e-6;

const double lifetime_n14   = 1919.e-6;
const double lifetime_n14_err = 15.e-6;

const double lifetime_o16   = 1796.e-6;
const double lifetime_o16_err =  3.e-6;

const double lifetime_o18   = 1844.e-6;
const double lifetime_o18_err =  5.e-6;

const double lifetime_gd   = (80.1+81.8)/2 * 1e-6;
const double lifetime_gd_err =  2.e-6;

// Isotopic fraction of C-13. Weighted average of T and GC measurements,
// which are only 0.2% different relative to each other and so really
// this figure is fine even for samples of events in other the T or only
// the GC, too.
const double f13 = 0.010919;

// But just to be obsessive, here's the weighted average of T and GC
// measurements for the HP region
const double f13_HP = 0.010921;


// isotopic fraction of N-15
const double f15 = 0.00364; // N-15 NNDC value

// isotopic fraction of O-17 and O-18
const double f17 = 0.00038; // O-17 NNDC value
const double f18 = 0.00205; // O-18 NNDC value

