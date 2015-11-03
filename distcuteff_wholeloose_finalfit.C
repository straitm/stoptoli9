#include "b12like_finalfit.C"
#include "distcuteff_include.C"

void distcuteff_wholeloose_finalfit()
{
  do_distcuteff(othercuts, "dist < 400",
         "TECHNOTE 7: Efficiency for 400mm cut, loose sample",
         "wholedet_dist400");
  do_distcuteff(othercuts, "dist < 300",
         "TECHNOTE 7: Efficiency for 300mm cut, loose sample",
         "wholedet_dist300");
  do_distcuteff(othercuts, "dist < 200",
         "TECHNOTE 7: Efficiency for 200mm cut, loose sample",
         "wholedet_dist200");
}
