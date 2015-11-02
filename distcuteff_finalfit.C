#include "b12like_finalfit.C"

const char * const target_cut =
  "dx**2+dy**2 < 1154**2 && "
  "abs(dz) < 1229 + 0.03*(1154 - sqrt(dx**2+dy**2))";

const char * const HPcut =
  "mx**2+my**2 < 1050**2 && mz > -1175 && "
  "abs(fez + 62*ivdedx/2 - 8847.2) < 1000 && rchi2 < 2";

const char * const othercuts =
  "timeleft > %f && miche < 12 && !earlymich && "
  "e > 4 && e < 15 && dt < %f";

// He-6 detector regions, copied from he6_finalfit.C because I couldn't
// think of a better way that I wanted to do.
int classi(const double x, const double y, const double z)
{
  const double r2 = x*x+y*y;
  const double r = sqrt(r2);
  const double az = abs(z);
  if(r2 > 1154*1154 || az > 1233 + 0.03*(1154-r)) return 4;
  if(r2 > 1068.*1068. || az > 1068.) return 3;
  if(r2 > 933.*933. || az > 933.) return 2;
  if(r2 > 740.*740. || az > 740) return 1;
  return 0;
}

ve geteff(const ve pass, const ve fail)
{
  ve answer;
  answer.val = pass.val/(pass.val+fail.val);
  const double N2 = pow(pass.val + fail.val,2);
  answer.err = sqrt(pow(pass.val*fail.err/N2, 2)
                  + pow(fail.val*pass.err/N2, 2));
  return answer;
}

void doit_wholedet_dist(const double distcut)
{
  const ve deltar_result = b12like_finalfit("eff",
      string(Form("%s && dist < %f", othercuts, distcut)).c_str(), false, false);
  const ve inverse_deltar_result = b12like_finalfit("eff",
      string(Form("%s && dist >= %f", othercuts, distcut)).c_str(), false, false);

  const ve result = geteff(deltar_result, inverse_deltar_result);

  printf("TECHNOTE 7: Efficiency for Dr = %.0fmm, loose cuts, whole det: (%.2f +- %.2f)%%\n",
         distcut, 100*result.val, 100*result.err);
  printf("const double wholedet_dist%.0feff = %f;\n", distcut, result.val);
  printf("const double wholedet_dist%.0feff_err = %f;\n", distcut, result.err);
}

void doit_other(const char * const basecut, const char * const effcut, 
                const char * const verbiage, const char * const headername = NULL)
{
  const ve cut =
    b12like_finalfit("eff", string(Form("%s &&   %s ", basecut, effcut)).c_str(), false, false);
  const ve anticut =
    b12like_finalfit("eff", string(Form("%s && !(%s)", basecut, effcut)).c_str(), false, false);

  printf("%s: (%f +- %f)%%\n", verbiage, 100*geteff(cut, anticut).val,
                                         100*geteff(cut, anticut).err);

  if(headername != NULL){
    printf("const double %s = %f;\n", headername, geteff(cut, anticut).val);
    printf("const double %s_err = %f;\n", headername, geteff(cut, anticut).val);
  }
}

void distcuteff_finalfit()
{
  doit_wholedet_dist(200);
  doit_wholedet_dist(300);
  doit_wholedet_dist(400);

  doit_other((string(othercuts) + "&& e > 8").c_str(), "dist < 400",
         "TECHNOTE 7: Efficiency for 400mm cut, loose sample, over 8MeV");
  doit_other((string(othercuts) + "&& e < 8").c_str(), "dist < 400",
         "TECHNOTE 7: Efficiency for 400mm cut, loose sample, under 8MeV");

  doit_other((string(othercuts) + "&& dz > 0").c_str(), "dist < 400",
         "TECHNOTE 7: Efficiency for 400mm, loose sample, top half");
  doit_other((string(othercuts) + "&& dz < 0").c_str(), "dist < 400",
         "TECHNOTE 7: Efficiency for 400mm, loose sample, bottom half");

  doit_other((string(othercuts) + " && " + string(target_cut)).c_str(),
         "dist < 400",
         "TECHNOTE 7: Efficiency for 400mm, loose sample, target");
  doit_other((string(othercuts)+"&& !("+string(target_cut)+")").c_str(),
         "dist < 400",
         "TECHNOTE 7: Efficiency for 400mm, loose sample, GC");

  doit_other((string(othercuts)+"&&!("+string(target_cut)+")").c_str(),
         "dist < 300",
         "TECHNOTE 7: Efficiency for 300mm, loose sample, GC"
         "gc_dist300eff");
  doit_other((string(othercuts) + " && " + string(target_cut)).c_str(),
         "dist < 300",
         "TECHNOTE 7: Efficiency for 300mm, loose sample, target",
         "targ_dist300eff");

  doit_other((string(othercuts) + " && classi(dx, dy, dz) == 0").c_str(),
         "dist < 200",
         "TECHNOTE 10.4: Efficiency for 200mm, loose sample, He-6 T0"
         "t0_dist200eff");
  doit_other((string(othercuts) + " && classi(dx, dy, dz) == 1").c_str(),
         "dist < 200",
         "TECHNOTE 10.4: Efficiency for 200mm, loose sample, He-6 T1"
         "t1_dist200eff");
  doit_other((string(othercuts) + " && classi(dx, dy, dz) == 2").c_str(),
         "dist < 200",
         "TECHNOTE 10.4: Efficiency for 200mm, loose sample, He-6 T2"
         "t2_dist200eff");
  doit_other((string(othercuts) + " && classi(dx, dy, dz) == 3").c_str(),
         "dist < 200",
         "TECHNOTE 10.4: Efficiency for 200mm, loose sample, He-6 T3"
         "t3_dist200eff");
  doit_other((string(othercuts) + " && classi(dx, dy, dz) == 4").c_str(),
         "dist < 200",
         "TECHNOTE 10.4: Efficiency for 200mm, loose sample, GC"
         "gc_dist200eff");
}
