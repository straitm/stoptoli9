#include "b12like_finalfit.C"
#include "distcuteff_include.C"

// Misnomer. I'm using the distcuteff framework to measure the b12like
// cut efficiency. Should I rename it everywhere?
void distcuteff_b12like_finalfit()
{
  // These efficiencies are functions of (at least) the distance cut.
  // The efficiency is higher for beta decays that reconstruct farther
  // from their parent stopping muon. I conjecture that this is because
  // these tend to be near the edge, which puts them mostly near
  // subsequent through-going muon tracks that have a short path length
  // and therefore don't lead to many detected neutrons, even if they
  // shower, reducing their likelihood. I haven't made the necessary
  // plots to check this.

  // Wasn't I clever to use a different set of cuts, including three
  // different values for the b12like cut, for all these analyses?

  do_distcuteff((string(othercuts) + "&&dist < 400").c_str(), "b12like<0.06",
    "TECHNOTE somewhere: B-14 efficiency for b12like < 0.06 cut, "
      "loose sample with dist < 400",
    "b12like006_dist400_eff");

  do_distcuteff((string(othercuts) + "&&dist < 400").c_str(), "b12like<0.02",
    "TECHNOTE somewhere: Li-8 efficiency for b12like < 0.02 cut, "
      "loose sample with dist < 400",
    "b12like002_dist400_eff");

  do_distcuteff((string(othercuts) +
    "&&dist < 200 && ttlastvalid>0.1 && ttlastmuon>1").c_str(), "b12like<0.02",
    "TECHNOTE somewhere: N-16 efficiency for b12like < 0.02 cut, "
       "loose sample with dist < 200, ttlastvalid > 0.1 and ttlastmuon > 1",
   "b12like002_dist200_ttlv01_ttlm1_eff");

  do_distcuteff((string(othercuts) +
    "&&dist < 200 && ttlastvalid>0.1 && ttlastmuon>1").c_str(), "b12like<0.4",
    "TECHNOTE somewhere: He-6 efficiency for b12like < 0.4 cut, "
       "loose sample with dist < 200, ttlastvalid > 0.1 and ttlastmuon > 1",
   "b12like040_dist200_ttlv01_ttlm1_eff");
}
