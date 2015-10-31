#include <stdio.h>
#include <math.h>
#include "consts.h"
#include "b12cutefficiency_finalfit.out.h"
#include "carbondenominators_finalfit.out.h"
#include "fullb12_finalfit.out.h"
#include "b12gamma_finalfit_0.out.h"
#include "b12gamma_finalfit_1.out.h"

void b12groundstate_finalfit()
{
  const double rate1 = probGammaOneRate;
  const double rate2 = probGammaTwoRate;
  const double rate3 = probGammaThreeRate;

  // Take the rate of all other states that cascade down
  // to the ground state of B-12 as twice the measured rate of 
  // the 3759 level.
  const double rate_rest = 2*probGammaFourRate;
 
  const double rerr1 =
    (fabs(probGammaOneRate_statup)+fabs(probGammaOneRate_statlo))/2;
  const double rerr2 =
    (fabs(probGammaTwoRate_statup)+fabs(probGammaTwoRate_statlo))/2;
  const double rerr3 =
    (fabs(probGammaThreeRate_statup)+fabs(probGammaThreeRate_statlo))/2;

  // Take the error on the rest as the error on a flat distribution
  // between zero and twice the nominal (twice the nominal -> 4 times
  // the measured 3759 central value) added in quadrature with the error
  // on the 3759 rate, also scaled up by a factor of two.
  const double rerr_rest_up = sqrt(pow(2*probGammaFourRate_statup,2) 
                                 + pow(2*probGammaFourRate/sqrt(12),2));
  const double rerr_rest_lo = sqrt(pow(2*probGammaFourRate_statlo,2) 
                                 + pow(2*probGammaFourRate/sqrt(12),2));
  const double rerr_rest = (rerr_rest_up+rerr_rest_lo)/2;


  const double cov12 = b12cc_12*rerr1*rerr2;
  const double cov13 = b12cc_13*rerr1*rerr3;
  const double cov23 = b12cc_23*rerr2*rerr3;

  // Conservatively (I think) take the covariance between each lower
  // line and the rest from the correlation coefficient between them and
  // 3759.
  const double cov1rest = b12cc_14*rerr1*rerr_rest;
  const double cov2rest = b12cc_24*rerr2*rerr_rest;
  const double cov3rest = b12cc_34*rerr3*rerr_rest;

  const double answercentral = b12totalrate
                                - rate1 - rate2 - rate3 - rate_rest;

  // normalization error due to mum counting, B-12 efficiency and mum
  // lifetime in C-12
  const double fsystmumcount = mum_count_e/mum_count;
  const double fsystb12eff = b12energyeff_e/b12energyeff;
  const double fmumlife = lifetime_c12_err/lifetime_c12;
  const double foverallsyst = sqrt(pow(fsystmumcount,2)
                                  +pow(fsystb12eff  ,2)
                                  +pow(fmumlife     ,2));
  const double overallsyst = answercentral*foverallsyst;

  
  const double staterror_lo = sqrt(pow(probGammaOneRate_statup,2)
                                 + pow(probGammaTwoRate_statup,2)
                                 + pow(probGammaThreeRate_statup,2)
                                 + pow(rerr_rest_up,2)
            + pow(b12totalrate_statlo, 2)
   + 2 * (cov12 + cov13 + cov23 + cov1rest + cov2rest + cov3rest));

  const double staterror_up = sqrt(pow(probGammaOneRate_statlo,2)
                                 + pow(probGammaTwoRate_statlo,2)
                                 + pow(probGammaThreeRate_statlo,2)
                                 + pow(rerr_rest_lo,2)
            + pow(b12totalrate_statup, 2)
   + 2 * (cov12 + cov13 + cov23 + cov1rest + cov2rest + cov3rest));

  const double fsharedgammasyst = n_c12cap_forb12gamma_additional_ferr;
  const double add_abs_syst = sqrt(
    pow(b12lineEsyst[0]*rate1,2)     + pow(b12lineNsyst[0]*rate1,2)
   +pow(b12lineEsyst[1]*rate2,2)     + pow(b12lineNsyst[1]*rate2,2)
   +pow(b12lineEsyst[2]*rate3,2)     + pow(b12lineNsyst[2]*rate3,2)
   +pow(b12lineEsyst[3]*rate_rest,2) + pow(b12lineNsyst[3]*rate_rest,2)
   +pow(fsharedgammasyst*(rate1 + rate2 + rate3 + rate_rest),2));
                                
  const double total_syst = sqrt(
    pow(add_abs_syst, 2)
   +pow(overallsyst,  2));

  const double total_err_up = sqrt(pow(total_syst,2) + pow(staterror_up,2));
  const double total_err_lo = sqrt(pow(total_syst,2) + pow(staterror_lo,2));

  printf("TECHNOTE results.tex probGammaZeroRateCent: B-12 ground "
         "state rate %.2f\n", answercentral);
  printf("%f\n", answercentral);
  printf("TECHNOTE results.tex probGammaZeroRate errors: B-12 ground "
         "state rate stat: +%.2f -%.2f syst: +-%.2f total: +%.2f -%.2f\n",
         staterror_up, staterror_lo, 
         total_syst,
         total_err_up, total_err_lo);
  printf("stat: +%f -%f syst: +-%f total: +%f -%f\n",
         staterror_up, staterror_lo, 
         total_syst,
         total_err_up, total_err_lo);
}
