#include "b12like_finalfit.C"
#include "distcuteff_include.C"

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
         "TECHNOTE 7: Efficiency for 400mm, loose sample, target");
  do_distcuteff((string(othercuts)+"&& !("+string(target_cut)+")").c_str(),
         "dist < 400",
         "TECHNOTE 7: Efficiency for 400mm, loose sample, GC");
}
