#include "b12like_finalfit.C"
#include "distcuteff_include.C"

// Misnomer. I'm using the distcuteff framework to measure the b12like
// cut efficiency. Should I rename it everywhere?
void distcuteff_b12like_finalfit()
{
  do_distcuteff(othercuts, "b12like<0.02",
         "TECHNOTE 7: Efficiency for b12like < 0.02 cut, loose sample",
         "b12like_eff_002");
  do_distcuteff(othercuts, "b12like<0.06",
         "TECHNOTE 7: Efficiency for b12like < 0.06 cut, loose sample",
         "b12like_eff_006");
  do_distcuteff(othercuts, "b12like<0.4",
         "TECHNOTE 7: Efficiency for b12like < 0.4 cut, loose sample",
         "b12like_eff_040");
}
