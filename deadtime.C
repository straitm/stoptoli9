#include "TF1.h"

TF1 *efffitg = new TF1("efffitg","min(1-0.0726, [0]*exp(-x/[1])+gaus(2))", 0, 800);
TF1 *efffith = new TF1("efffith","min(exp(-5.5/179) - exp(-800/179.), [0]*exp(-x/[1])+gaus(2))",0, 800);

double eff(const float fidoqid, const double x, const double y,
           const double z, const bool early = false, const bool reallyearlyh = true)
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

  static double earlyprob = 1-exp(-5.5/179);

  // In the target
  if(fabs(z) < 1229 + 0.03*(1150-r) && r < 1150) 
    return efffitg->Eval(fidoqid/8300) + early*0.0726;

  // In the target acrylic -- 0.5449 chance of capturing in the target
  if(fabs(z) < 1237 + 0.03*(1158-r) && r < 1158) 
    return (efffith->Eval(fidoqid/8300) + reallyearlyh*early*earlyprob)*
      (1-0.5449)+
    (efffitg->Eval(fidoqid/8300) + early*0.0726*(0.84+reallyearlyh*0.16))*
      (0.5449);

  // In the GC
  return efffith->Eval(fidoqid/8300) + reallyearlyh*early*earlyprob;
}
