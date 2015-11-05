#include "b12like_finalfit.C"
#include "distcuteff_include.C"
#include <string>
using std::string;

const string t0_cut = "(dx**2+dy**2 < 740**2 && abs(dz) < 740)";
const string t1_cut = "(dx**2+dy**2 < 933**2 && abs(dz) < 933)";
const string t2_cut = "(dx**2+dy**2 < 1068**2 && abs(dz) < 1068)";

// Partially a misnomer.  The first thing we do is find the number of
// stopping muons in each region, which has nothing to do with a 
// distance cut.  Then we evaluate distance cut efficiencies for
// each region.
void distcuteff_he6_finalfit()
{
  const double allregions = b12like_finalfit("eff",othercuts,0,0).val;
  const double t0 = b12like_finalfit("eff",
    (string(othercuts)+"&&"+t0_cut).c_str(),0,0).val;
  const double t1 = b12like_finalfit("eff",
    (string(othercuts)+"&&!("+t0_cut+")&&"+t1_cut).c_str(),0,0).val;
  const double t2 = b12like_finalfit("eff",
    (string(othercuts)+"&&!("+t1_cut+")&&"+t2_cut).c_str(),0,0).val;
  const double t3 = b12like_finalfit("eff",
    (string(othercuts)+"&&!("+t2_cut+")&&"+target_cut).c_str(),0,0).val;
  const double gc = allregions - t3 - t2 - t1 - t0;

  printf("const double he6mus[5] = { %f, %f, %f, %f, %f }; "
    "// The number of stopping muons in each region, from the inside "
    "out, with an arbitrary overall normalization\n",
  t0/allregions,t1/allregions,t2/allregions,t3/allregions,gc/allregions);

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
