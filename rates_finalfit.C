void rates_finalfit()
{
  const double days = 489.509;

  const double c_atomic_cap_hp = 7.149e5;
  const double c_atomic_cap_hp_err = 0.054e5;
  
  const double c_atomic_cap_loose = 2.309e6;
  const double c_atomic_cap_loose_err = 0.032e6;

  const double mulife = 2196.9811e-6;

  // Weighted average of T and GC measurements, which are only 0.2%
  // different relative to each other and so really this figure is fine
  // even for samples of events in other the T or only the GC, too.
  const double f13 = 0.010919;

  const double f15 = 0.00364; // N-15 NNDC value

  const double f17 = 0.00038; // O-17 NNDC value
  const double f18 = 0.00205; // O-18 NNDC value

  const double c12_atomic_cap_hp = c_atomic_cap_hp*(1-f13);
  const double c12_atomic_cap_loose = c_atomic_cap_loose*(1-f13);
  const double c13_atomic_cap_hp = c_atomic_cap_hp*f13;
  const double c13_atomic_cap_loose = c_atomic_cap_loose*f13;

  const double c12_atomic_cap_hp_err = c_atomic_cap_hp_err*(1-f13);
  const double c12_atomic_cap_loose_err = c_atomic_cap_loose_err*(1-f13);
  const double c13_atomic_cap_hp_err = c_atomic_cap_hp_err*f13;
  const double c13_atomic_cap_loose_err = c_atomic_cap_loose_err*f13;

  const double lifetime_c12  = 2028.e-6;
  const double lifetime_c12_err = 2.e-6;

  const double lifetime_c13  = 2037.e-6;
  const double lifetime_c13_err = 8.e-6;

  const double lifetime_n14   = 1919.e-6;
  const double lifetime_n14_err = 15.e-6;

  const double lifetime_o16   = 1796.e-6;
  const double lifetime_o16_err =  3.e-6;

  const double lifetime_o18   = 1844.e-6;
  const double lifetime_o18_err =  5.e-6;

  const double lifetime_gd   = (80.1+81.8)/2 * 1e-6;
  const double lifetime_gd_err =  2.e-6;

  const double capprob_c12 = 1-lifetime_c12/mulife;
  const double capprob_c12_err =-(1-(lifetime_c12+lifetime_c12_err)/mulife)/2
                              +(1-(lifetime_c12-lifetime_c12_err)/mulife)/2;

  const double capprob_c13 = 1-lifetime_c13/mulife;
  const double capprob_c13_err =-(1-(lifetime_c13+lifetime_c13_err)/mulife)/2
                              +(1-(lifetime_c13-lifetime_c13_err)/mulife)/2;

  const double capprob_n14 = 1-lifetime_n14/mulife;
  const double capprob_n14_err =-(1-(lifetime_n14+lifetime_n14_err)/mulife)/2
                              +(1-(lifetime_n14-lifetime_n14_err)/mulife)/2;

  const double capprob_o16 = 1-lifetime_o16/mulife;
  const double capprob_o16_err =-(1-(lifetime_o16+lifetime_o16_err)/mulife)/2
                              +(1-(lifetime_o16-lifetime_o16_err)/mulife)/2;

  const double capprob_o18 = 1-lifetime_o18/mulife;
  const double capprob_o18_err =-(1-(lifetime_o18+lifetime_o18_err)/mulife)/2
                              +(1-(lifetime_o18-lifetime_o18_err)/mulife)/2;

  const double capprob_gd = 1-lifetime_gd/mulife;
  const double capprob_gd_err =-(1-(lifetime_gd+lifetime_gd_err)/mulife)/2
                              +(1-(lifetime_gd-lifetime_gd_err)/mulife)/2;

  const double c12_nuclear_cap_hp = c12_atomic_cap_hp*capprob_c12;
  const double c12_nuclear_cap_loose = c12_atomic_cap_loose*capprob_c12;
  const double c13_nuclear_cap_hp = c13_atomic_cap_hp*capprob_c13;
  const double c13_nuclear_cap_loose = c13_atomic_cap_loose*capprob_c13;

  const double c12_nuclear_cap_hp_err = c12_nuclear_cap_hp*sqrt(
    pow(c12_atomic_cap_hp_err/c12_atomic_cap_hp, 2) +
    pow(capprob_c12_err/capprob_c12, 2));

  const double c12_nuclear_cap_loose_err = c12_nuclear_cap_loose*sqrt(
    pow(c12_atomic_cap_loose_err/c12_atomic_cap_loose, 2) +
    pow(capprob_c12_err/capprob_c12, 2));

  const double c13_nuclear_cap_hp_err = c13_nuclear_cap_hp*sqrt(
    pow(c13_atomic_cap_hp_err/c13_atomic_cap_hp, 2) +
    pow(capprob_c13_err/capprob_c13, 2));

  const double c13_nuclear_cap_loose_err = c13_nuclear_cap_loose*sqrt(
    pow(c13_atomic_cap_loose_err/c13_atomic_cap_loose, 2) +
    pow(capprob_c13_err/capprob_c13, 2));

#define ABE(x) x, x##_err, x##_err/x*100

  printf("Muons in the high-purity sample: %.4f +- %.4f (%.1f%%)\n",
         ABE(c_atomic_cap_hp));
  printf("Muons in the loose sample: %.4f +- %.4f (%.1f%%)\n",
         ABE(c_atomic_cap_loose));

  puts("");

  printf("C-12 nuclear capture prob: %.4f +- %.4f (%.1f%%)\n",
         ABE(capprob_c12));
  printf("C-13 nuclear capture prob: %.4f +- %.4f (%.1f%%)\n", 
         ABE(capprob_c13));
  printf("N-14 nuclear capture prob: %.4f +- %.4f (%.1f%%)\n", 
         ABE(capprob_n14));
  printf("O-16 nuclear capture prob: %.4f +- %.4f (%.1f%%)\n", 
         ABE(capprob_o16));
  printf("O-18 nuclear capture prob: %.4f +- %.4f (%.1f%%)\n", 
         ABE(capprob_o18));
  printf("Gd   nuclear capture prob: %.4f +- %.4f (%.1f%%)\n", 
         ABE(capprob_gd));

  puts("");
 
  printf("C-12 mu- lifetime: %.4g +- %.4g (%.2f%%)\n",
         ABE(lifetime_c12));
  printf("C-13 mu- lifetime: %.4g +- %.4g (%.2f%%)\n", 
         ABE(lifetime_c13));
  printf("N-14 mu- lifetime: %.4g +- %.4g (%.1f%%)\n", 
         ABE(lifetime_n14));
  printf("O-16 mu- lifetime: %.4g +- %.4g (%.2f%%)\n", 
         ABE(lifetime_o16));
  printf("O-18 mu- lifetime: %.4g +- %.4g (%.2f%%)\n", 
         ABE(lifetime_o18));
  printf("Gd   mu- lifetime: %.4g +- %.4g (%.1f%%)\n", 
         ABE(lifetime_gd));

  puts("");

  const double c12_atomic_cap_hp_days = c12_atomic_cap_hp/days;
  const double c12_atomic_cap_hp_days_err = c12_atomic_cap_hp_err/days;
  const double c13_atomic_cap_hp_days = c13_atomic_cap_hp/days;
  const double c13_atomic_cap_hp_days_err = c13_atomic_cap_hp_err/days;
 
  const double c12_atomic_cap_loose_days = c12_atomic_cap_loose/days;
  const double c12_atomic_cap_loose_days_err=c12_atomic_cap_loose_err/days;
  const double c13_atomic_cap_loose_days = c13_atomic_cap_loose/days;
  const double c13_atomic_cap_loose_days_err=c13_atomic_cap_loose_err/days;
 
  const double c12_nuclear_cap_hp_days = c12_nuclear_cap_hp/days;
  const double c12_nuclear_cap_hp_days_err=c12_nuclear_cap_hp_err/days;
  const double c13_nuclear_cap_hp_days = c13_nuclear_cap_hp/days;
  const double c13_nuclear_cap_hp_days_err=c13_nuclear_cap_hp_err/days;
 
  const double c12_nuclear_cap_loose_days=c12_nuclear_cap_loose/days;
  const double c12_nuclear_cap_loose_days_err=c12_nuclear_cap_loose_err/days;
  const double c13_nuclear_cap_loose_days=c13_nuclear_cap_loose/days;
  const double c13_nuclear_cap_loose_days_err=c13_nuclear_cap_loose_err/days;
 
  printf("Atomic C-12 captures, high purity/day: %.0f +- %.0f (%.2f%%)\n",
         ABE(c12_atomic_cap_hp_days));
  printf("Nuclear C-12 captures, high purity/day: %g +- %.2g (%.1f%%)\n",
         ABE(c12_nuclear_cap_hp_days));
  puts("");
  printf("Atomic C-13 captures, high purity/day: %.2f +- %.2f (%.2f%%)\n",
         ABE(c13_atomic_cap_hp_days));
  printf("Nuclear C-13 captures, high purity/day: %.3f +- %.2g (%.1f%%)\n",
         ABE(c13_nuclear_cap_hp_days));

  puts("");
  puts("");
 
  printf("Atomic C-12 captures, loose/day: %.4g +- %.2g (%.1f%%)\n",
         ABE(c12_atomic_cap_loose_days));
  printf("Nuclear C-12 captures, loose/day: %.4g +- %.2g (%.1f%%)\n",
         ABE(c12_nuclear_cap_loose_days));
  puts("");
  printf("Atomic C-13 captures, loose/day: %.1f +- %.2g (%.1f%%)\n",
         ABE(c13_atomic_cap_loose_days));
  printf("Nuclear C-13 captures, loose/day: %.2f +- %.2g (%.1f%%)\n",
         ABE(c13_nuclear_cap_loose_days));



}
