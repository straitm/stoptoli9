#include <string>
using std::string;

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

void do_distcuteff(const char * const basecut, const char * const effcut, 
                const char * const verbiage, const char * const headername = NULL)
{
  runningforefficiency = true;

  const ve cut =
    b12like_finalfit("eff", string(Form("%s &&   %s ", basecut, effcut)).c_str(), false, false);
  const ve anticut =
    b12like_finalfit("eff", string(Form("%s && !(%s)", basecut, effcut)).c_str(), false, false);

  printf("%s: (%f +- %f)%%\n", verbiage, 100*geteff(cut, anticut).val,
                                         100*geteff(cut, anticut).err);

  if(headername != NULL){
    printf("const double %s = %f;\n", headername, geteff(cut, anticut).val);
    printf("const double %s_err = %f;\n", headername, geteff(cut, anticut).err);
  }
}

