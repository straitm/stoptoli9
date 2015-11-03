#include "b12like_finalfit.C"
#include "distcuteff_include.C"

void distcuteff_byenergy_finalfit()
{
  do_distcuteff((string(othercuts) + "&& e > 8").c_str(), "dist < 400",
         "TECHNOTE 7: Efficiency for 400mm cut, loose sample, over 8MeV");
  do_distcuteff((string(othercuts) + "&& e < 8").c_str(), "dist < 400",
         "TECHNOTE 7: Efficiency for 400mm cut, loose sample, under 8MeV");
}
