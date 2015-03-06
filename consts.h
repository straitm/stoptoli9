const char * const RED     = "\033[31;1m"; // bold red
const char * const CLR     = "\033[m"    ; // clear

const double livetime = 489.509, n_c12cap = 350.079246, n_c13cap=3.517328;
const double n_o16cap_betan = 8.886; // after acrylic efficiencies - central value
const double n_o16cap_beta  = 9.639; // after acrylic efficiencies - central value
const double n_n14cap = (0.298+0.180)/2;
const double n_c12captarget = 128.307666;
const char * const rootfile3up = "/cp/s4/strait/fullfido-300s-3-25MeV-20150219.root";
const char * const rootfile0up = "/cp/s4/strait/fullfido-100s-0-25MeV-20150219.root";

const double gd_fraction = 0.851;

const double neff_dt_targ = 0.5862; // same as Gd eff
const double neff_dt_gc   = 0.9060;

// H efficiency is lower than GC efficiency by a little 
// because some H is in the Target.
const double neff_dt_h = 
  (
  (n_c12cap - n_c12captarget)      * neff_dt_gc
+ n_c12captarget * (1-gd_fraction) * neff_dt_targ
  )/
  (n_c12cap - gd_fraction*n_c12captarget);

const double neff_dt_avg = ((n_c12cap - n_c12captarget)*neff_dt_gc
                            + n_c12captarget * neff_dt_targ)/n_c12cap;

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

const double neff_dr_800_avg =
  (
  (n_c12cap - n_c12captarget*gd_fraction) * neff_dr_800_h
  +           n_c12captarget*gd_fraction  * neff_dr_800_gd
  )/n_c12cap;

const double neff_dr_1000_avg =
  (
  (n_c12cap - n_c12captarget*gd_fraction) * neff_dr_1000_h
  +           n_c12captarget*gd_fraction  * neff_dr_1000_gd
  )/n_c12cap;

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
