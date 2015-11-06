#include "b12like_finalfit.C"
#include "distcuteff_include.C"

const string targvesin_cut  = "(dx**2+dy**2 < 950**2 && abs(dz) < 1029)";
const string targvesout_cut = "(dx**2+dy**2 < 1350**2 && abs(dz) < 1429)";

const string gcves_cut = "(dx**2+dy**2 > 1508**2 || abs(dz) > 1586)";

void distcuteff_targetgc_finalfit()
{
  do_distcuteff((string(othercuts) + " && " + string(target_cut)).c_str(),
         "dist < 300",
         "TECHNOTE 8.1: Efficiency for 300mm, loose sample, target",
         "targ_dist300eff");
  do_distcuteff((string(othercuts)+"&& !("+string(target_cut)+")").c_str(),
         "dist < 300",
         "TECHNOTE 8.1: Efficiency for 300mm, loose sample, GC",
         "gc_dist300eff");

  do_distcuteff((string(othercuts) + " && " + string(target_cut)).c_str(),
         "dist < 400",
         "TECHNOTE 7: Efficiency for 400mm, loose sample, target",
         "targ_dist400eff");
  do_distcuteff((string(othercuts)+"&& !("+string(target_cut)+")").c_str(),
         "dist < 400",
         "TECHNOTE 7: Efficiency for 400mm, loose sample, GC",
         "gc_dist400eff");

  do_distcuteff((targvesout_cut+"&& !("+targvesin_cut+")").c_str(),
         "dist < 400",
         "TECHNOTE 7: Efficiency for 400mm, loose sample, near the target vessel",
         "targves_dist400eff");
 
  do_distcuteff(gcves_cut.c_str(),
         "dist < 400",
         "TECHNOTE 7: Efficiency for 400mm, loose sample, near the GC vessel",
         "gcves_dist400eff");
}
