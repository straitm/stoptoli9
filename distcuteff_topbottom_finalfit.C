#include "b12like_finalfit.C"
#include "distcuteff_include.C"

void distcuteff_topbottom_finalfit()
{
  do_distcuteff((string(othercuts) + "&& dz > 0").c_str(), "dist < 400",
         "TECHNOTE 7: Efficiency for 400mm, loose sample, top half");
  do_distcuteff((string(othercuts) + "&& dz < 0").c_str(), "dist < 400",
         "TECHNOTE 7: Efficiency for 400mm, loose sample, bottom half");
}
