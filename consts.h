const char * const RED     = "\033[31;1m"; // bold red
const char * const CLR     = "\033[m"    ; // clear

const double livetime = 489.509, n_c12cap = 350.149336, n_c13cap=3.518032;
const double n_o16cap_betan = 9.1; // after acrylic efficiencies - central value
const double n_o16cap_beta  = 9.8; // after acrylic efficiencies - central value
const double n_n14cap = (0.30+0.18)/2;
const double n_c12captarget = 128.850030;
const char * const rootfile3up = "/cp/s4/strait/fullfido-300s-3-25MeV-20150219.root";
const char * const rootfile0up = "/cp/s4/strait/fullfido-100s-0-25MeV-20150219.root";

const double neff_dt_targ = 0.5862;
const double neff_dt_gc   = 0.9060;

const double neff_dt_avg = ((n_c12cap - n_c12captarget)*neff_dt_gc
                            + n_c12captarget * neff_dt_targ)/n_c12cap;


