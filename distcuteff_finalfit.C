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
      Form("%s && dist < %f", othercuts, distcut), false, false);
  const ve inverse_deltar_result = b12like_finalfit("eff",
      Form("%s && dist >= %f", othercuts, distcut), false, false);

  const ve result = geteff(deltar_result, inverse_deltar_result);

  printf("TECHNOTE 7: Efficiency for Dr = %.0fmm, loose cuts, whole det: (%.2f +- %.2f)%%\n",
         distcut, 100*result.val, 100*result.err);
  printf("const double wholedet_dist%.0feff = %f;\n", distcut, result.val);
  printf("const double wholedet_dist%.0feff_err = %f;\n", distcut, result.err);
}

void doit_other(const char * const basecut, const char * const effcut, 
                const char * const verbiage)
{
  const ve cut =
    b12like_finalfit("eff", Form("%s &&   %s ", basecut, effcut), false, false);
  const ve anticut =
    b12like_finalfit("eff", Form("%s && !(%s)", basecut, effcut), false, false);

  printf("%s: (%f +- %f)%%\n", verbiage, 100*geteff(cut, anticut).val,
                                         100*geteff(cut, anticut).err);
}

void distcuteff_finalfit()
{
  doit_wholedet_dist(200);
  doit_wholedet_dist(300);
  doit_wholedet_dist(400);

  doit_other((string(othercuts) + "e > 8").c_str(), "dist < 400",
         "TECHNOTE 7: Efficiency for 400mm cut over 8MeV");
  doit_other((string(othercuts) + "e < 8").c_str(), "dist < 400",
         "TECHNOTE 7: Efficiency for 400mm cut under 8MeV");

  doit_other((string(othercuts) + "dz > 0").c_str(), "dist < 400",
         "TECHNOTE 7: Efficiency for 400mm, top half");
  doit_other((string(othercuts) + "dz < 0").c_str(), "dist < 400",
         "TECHNOTE 7: Efficiency for 400mm, bottom half");

  doit_other((string(othercuts) + target_cut).c_str(), "dist < 400",
         "TECHNOTE 7: Efficiency for 400mm, target");
  doit_other((string(othercuts) + "!(" + string(target_cut) + ")").c_str(), "dist < 400",
         "TECHNOTE 7: Efficiency for 400mm, GC");
}
