#include "TF1.h"
#include <math.h>
#include "consts.h"

TF1 *efffitg = new TF1("efffitg","min(1-0.0726, [0]*exp(-x/[1])+gaus(2))", 0, 800);
TF1 *efffith = new TF1("efffith","min(exp(-5.5/179.) - exp(-800./179.), [0]*exp(-x/[1])+gaus(2))",0, 800);

double neff_dr_800(const double x, const double y, const double z)
{
  if(x*x+y*y < 1154*1154 &&
     fabs(z) < 1229 + 0.03*(1154 - sqrt(x*x+y*y)))
    return neff_dr_800_targ;
  else
    return neff_dr_800_h;
}

/* fidoqid == 0 is taken as a code that this is near detector data with
DDR, because in this case fidoqid has no effect (I think). Gd captures
have no dt inefficiency in this case, but H captures have some IV-energy
dependent inefficiency starting at approximately fqiv == 200e3 for t
< 50mus. Why the efficiency is not exactly zero is a mystery to me,
although could presumably be learned from Manu. For now, I'm going to
cut events in this range and therefore be able to set their efficiency
to zero. */
double neff_dt(const float fidoqid, const float fidoqiv, const double x,
               const double y, const double z, const bool early = false)
{
// With the logisitic function as deadtime cutoff
   static bool firsttime = true;
   if(firsttime){
     firsttime = false;
     efffitg->SetParameter(0,1.28922);
     efffitg->SetParameter(1,348.103);
     efffitg->SetParameter(2,0.0757327);
     efffitg->SetParameter(3,194.879);
     efffitg->SetParameter(4,22.29);

     efffith->SetParameter(0,1.00176);
     efffith->SetParameter(1,2419.96);
     efffith->SetParameter(2,0.0121678);
     efffith->SetParameter(3,195.013);
     efffith->SetParameter(4,15.0475);
   }

// With the error function as deadtime cutoff
/*
   efffith->SetParameter(0,1.00578);
   efffith->SetParameter(1,2366.49);
   efffith->SetParameter(2,0.00463461);
   efffith->SetParameter(3,218.553);
   efffith->SetParameter(4,49.9893);
   efffitg->SetParameter(0,1.28925);
   efffitg->SetParameter(1,348.029);
   efffitg->SetParameter(2,0.0756832);
   efffitg->SetParameter(3,194.095);
   efffitg->SetParameter(4,23.1781);
*/

  const double r2 = x*x + y*y;
  const double r = sqrt(r2);

  static const double htau = 179.;

  // This is the FD value, but if ND, fidoqid is always passed in as
  // zero, so it doesn't matter.
  static const double qid_per_mev = 8300;

  // Can't possibly see a capture until the next accepted trigger
  // window, and then can see "early" captures up until 5.5mus, by
  // definition. VERY roughly handle the other sort of deadtime in the
  // FD "early" window.
  static const double hearlyprobbase = (exp(-0.512/htau)-exp(-5.5/htau));
  double hearlyprob = fidoqid > 0?hearlyprobbase/2.:hearlyprobbase;

  // Assume that no Gd captures are lost to being in the same
  // trigger window as the previous event, since the probability
  // of very early capture is very low.
  static double gdearlyprob =  1-efffitg->Eval(0);

  // If this is the near detector (fidoqid == 0) and over the approximate
  // threshold for the DDR to reduce the efficiency of Hn events, reduce
  // the efficiency accordingly.  Otherwise, this has no effect.
  static const double nearhmult = 1 - (exp(-5.5/htau) - exp(-50./htau));
  const double hmult = fidoqid == 0 && fidoqiv > 200e3? nearhmult: 1;
  const double earlyhmult = fidoqid == 0 && fidoqiv > 200e3? 0: 1;

  const double targeff =
    /* For the far detector, only time matters, but at the near detector,
    the DDR makes the capture type matter too. */
    (efffitg->Eval(fidoqid/qid_per_mev) + early*gdearlyprob)
     *gd_fraction
                                        /* XXX is this close enough?  */
                                        /* sic on h and gd mix, but   */
                                        /* does the h capture follow  */
                                        /* gd timing at early time?   */
    +(hmult*efffitg->Eval(fidoqid/qid_per_mev) + earlyhmult*early*gdearlyprob)
     *(1-gd_fraction);

  const double gceff =
      hmult*efffith->Eval(fidoqid/qid_per_mev) + earlyhmult*early*hearlyprob;

  // In the target
  if(fabs(z) < 1229 + 0.03*(1150-r) && r < 1150)
    return targeff;

  // In the target acrylic -- 0.5449 chance of capturing in the target
  static const double acry_nt_prob = 0.5449;
  if(fabs(z) < 1237 + 0.03*(1158-r) && r < 1158)
    return gceff*(1-acry_nt_prob) + targeff*acry_nt_prob;

  // In the GC
  return gceff;
}
