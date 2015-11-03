#include "b12like_finalfit.C"
#include "distcuteff_include.C"
#include <string>
using std::string;

const string t0_cut = "(dx**2+dy**2 < 740**2 && abs(dz) < 740)";
const string t1_cut = "(dx**2+dy**2 < 933**2 && abs(dz) < 933)";
const string t2_cut = "(dx**2+dy**2 < 1068**2 && abs(dz) < 1068)";

void distcuteff_he6_finalfit()
{
  do_distcuteff((string(othercuts)+"&&"+t0_cut).c_str(),
         "dist < 200",
         "TECHNOTE 10.4: Efficiency for 200mm, loose sample, He-6 T0",
         "t0_dist200eff");
  do_distcuteff((string(othercuts) + "&&"+t1_cut+"&&!"+t0_cut).c_str(),
         "dist < 200",
         "TECHNOTE 10.4: Efficiency for 200mm, loose sample, He-6 T1",
         "t1_dist200eff");
  do_distcuteff((string(othercuts) + " &&"+t2_cut+"&&!"+t1_cut).c_str(),
         "dist < 200",
         "TECHNOTE 10.4: Efficiency for 200mm, loose sample, He-6 T2",
         "t2_dist200eff");
  do_distcuteff((string(othercuts) + " &&"+string(target_cut)+"&&!"+t2_cut).c_str(),
         "dist < 200",
         "TECHNOTE 10.4: Efficiency for 200mm, loose sample, He-6 T3",
         "t3_dist200eff");
  do_distcuteff((string(othercuts)+"&&!("+string(target_cut)+")").c_str(),
         "dist < 200",
         "TECHNOTE 10.4: Efficiency for 200mm, loose sample, GC",
         "gc_dist200eff");
}
