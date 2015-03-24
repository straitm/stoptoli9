TF1 *efffitg = new TF1("efffitg","min(1-0.0726, [0]*exp(-x/[1])+gaus(2))", 0, 800);
TF1 *efffith = new TF1("efffith","min(exp(-5.5/179) - exp(-800/179.), [0]*exp(-x/[1])+gaus(2))",0, 800);

double eff(const float fidoqid, const double x, const double y,
           const double z, const bool early = false)
{
// With the logisitic function as deadtime cutoff
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


  return fabs(z) < 1233 + 0.03*(1154-sqrt(x*x+y*y)) && x*x+y*y < 1154*1154? 
         efffitg->Eval(fidoqid/8300) + early*0.0726:
         efffith->Eval(fidoqid/8300) + early*(1-exp(-5.5/179));
}
