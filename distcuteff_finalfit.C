#include "b12like_finalfit.C"

void distcuteff_finalfit()
{
  const char * const target_cut =
    "dx**2+dy**2 < 1154**2 && "
    "abs(dz) < 1229 + 0.03*(1154 - sqrt(dx**2+dy**2))";

  const char * const HPcut =
    "mx**2+my**2 < 1050**2 && mz > -1175 && "
    "abs(fez + 62*ivdedx/2 - 8847.2) < 1000 && rchi2 < 2";

  const char * const othercuts =
    "timeleft > %f && miche < 12 && !earlymich && "
    "e > 4 && e < 15 && dt < %f";

  const ve loose_result =
    b12like_finalfit("loose for delta r", othercuts, false, false);
  const ve deltar400_result =
    b12like_finalfit("loose for delta r", Form("%s && dist < 400", othercuts), false, false);

  const double eff = deltar400_result.val/loose_result.val;

  const double err = sqrt(deltar400_result.val*(loose_result.val-deltar400_result.val))/
                     loose_result.val;

  printf("TECHNOTE: Efficiency for Dr = 400m, loose cuts, whole det: (%.2f +- %.2f)%%\n",
         eff, err);
  printf("const double eff_deltar_loose_wholedet_400mm = %f;\n", eff);
  printf("const double eff_err_deltar_loose_wholedet_400mm = %f;\n", err);
}
