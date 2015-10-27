const char * const RED     = "\033[31;1m"; // bold red
const char * const CLR     = "\033[m"    ; // clear

const char * const rootfile3up_extended = // XXX now out of date
  "/cp/s4/strait/fullfido-300s-3-25MeV-20150924-3rdpub+post3rdpub.root";
const char * const rootfile3up =
  "/cp/s4/strait/fullfido-300s-3-25MeV-20151023.root";
const char * const rootfile0up =
  "/cp/s4/strait/fullfido-100s-0-25MeV-20151023.root";
const char * const rootfilethru = // XXX now out of date -- but doesn't matter, since it is only for energy calibration
  "/cp/s4/strait/fullfido-stopandthru-1ms-3-25MeV-20150325.root";
const char * const rootfilemusic =
  "/cp/s4/strait/DC_FD_Muon_Surface_Dist_Ecut_both.root";

const double mum_frac = 0.4410,
             mum_frac_err = 0.0032,
             mum_contamination = 0.0028,
             mum_contamination_err = 0.0019;

const double livetime = 489.509;
const double rrmlivetimes[3] = { 7.570, 201.505, 280.434 };

// c12 and c13 captures per day
#include "loosecaptures_finalfit_out.h"

// From dcfluids.ods
const double n_o16cap_beta  = 9.8; // after acrylic eff - central value
const double n_o16cap_betan = 9.0; // after acrylic eff - central value

const double n_o16cap_beta_hp  = 5.0;
const double n_o16cap_betan_hp = 5.0; // yes, the same as beta

const double n_n14cap = (0.19+0.31)/2;

const double n_n14cap_hp = (0.11+0.18)/2;

const double n_c12captarget_hp = 131.8504 * 6743.417186/7882.377282;

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

// Subsequent muon veto efficiency (an efficiency on the isotope decay,
// NOT on the muon), for the hard cut imposed on events in order to get
// into the ntuples of 0.5ms.
// 
// This is seasonally dependent! But the variation is only
// 0.9773-0.9780. XXX. The blessed value of Sept 2015 for this is 0.981, but
// a more careful evaluation gets that range.  Hopefully I can update it
// for the paper.  The DC3rdpub period average is 0.977614.
const double sub_muon_eff05 = 0.981;
//const double sub_muon_eff05 = 0.977614;

// XXX. The more careful (but unblessed) value is 0.955744 for the
// DC3rdpub period.
const double sub_muon_eff10 = 0.962;
//const double sub_muon_eff10 = 0.955744;

// I made this up, with some studies to back it.
const double f_neff_dt_error = 0.01;

const double neff_dt_targ = 0.5724; // same as Gd eff
const double neff_dt_gc   = 0.8968;

// Might be useful to break this into Target and GC, but it is mostly
// Target, so maybe not.
const double neff_dt_highpurity = 0.6193; // only late neutrons

// Efficiency for four neutrons (not the same as the efficiency for one
// neutron**4)
const double n4eff_dt_targ        = 0.1697;
const double n4eff_dt_targ_wearly = 0.2514;
const double n4eff_dt_gc          = 0.6616;

// Efficiency for three out of four neutrons
const double n3of4eff_dt_targ        = 0.2801;
const double n3of4eff_dt_targ_wearly = 0.3021;
const double n3of4eff_dt_gc          = 0.2719;

// Efficiency for two out of four neutrons
const double n2of4eff_dt_targ        = 0.2870;
const double n2of4eff_dt_targ_wearly = 0.2606;
const double n2of4eff_dt_gc          = 0.0590;

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

const double wholedet_dist400eff = 0.9148;
const double wholedet_dist300eff = 0.8121;
const double wholedet_dist200eff = 0.5343;
