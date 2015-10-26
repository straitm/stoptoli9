/*
 * Terrible name, dcfluids_finalfit is, but what it does is finds
 * the number of nitrogen and oxygen captures like the spreadsheet
 * dcfluids.ods used to. It depends a lot on the fluid compositions, so
 * maybe it isn't that bad of a name. Well, it really depends more on
 * the acrylic composition and wild guesses, but whatever.
 *
 * This code is a direct translation from the spreadsheet. Please
 * consider it machine-generated even though I did it by hand for lack
 * of knowing how to do an automatic conversion.
 */

#include "consts.h"

const double f13 = 0.010921;

const double hydrogenmass = 1.008;
const double carbonmass = 12.011;
const double nitrogenmass = 14.007;
const double oxygenmass = 15.999;

const double mulife = 2196.9811e-6;

const double lifetime_c12 = 2028.e-6;
const double lifetime_c13 = 2037.e-6;

const double lifetime_c = lifetime_c12*(1-f13)+lifetime_c13*f13;
const double capprob_c = 1-lifetime_c/mulife;

const double lifetime_o16 = 1796.e-6;
const double capprob_o16 = 1-lifetime_o16/mulife;

const double probrat = capprob_o16/capprob_c;


double looseE28();
double looseH50();
double looseD28();
double looseI50();
double looseG50();
double looseE26();
double looseD23();
double looseC16();
double looseC17();
double looseC18();
double looseD20();
double looseD21();
double looseD22();
double looseD24();
double looseD25();
double looseD26();
double looseB20();
double looseB21();
double looseB22();
double looseB24();
double looseB25();
double looseI18();
double looseD8();
double looseC8();
double looseB11();
double looseC11();
double looseE3();
double looseI3();
double looseI16();
double looseI17();
double looseB23();
double looseI13();
double looseC10();
double looseC4();
double looseD6();
double looseD5();
double looseD7();
double looseF4();
double looseH4();
double looseC3();
double looseH3();
double looseB10();
double looseF2();
double looseH2();
double looseH5();
double looseI10();
double looseI12();
double looseI14();
double looseI19();
double looseJ14();
double looseJ19();
double looseB26();
double looseD27();
double looseC26();
double looseE27();
double looseG49();
double looseH49();
double looseI49();





const double gc_inner_r = 1708;
const double gc_inner_h = 1786;
const double nt_inner_r = 1150;
const double nt_inner_h = 1229;
const double lidslope = 0.03;
const double mass_immersed_acrylic = 432.; // in kg
const double kg_oxygen_ppo_nt = 75.25;
const double mass_halfimmersed_acrylic = 814; // in kg

// XXX needs to be generated in loosecaptures_finalfit.C
const double n_c12cap_target = 129.979327;


double looseE26(){ return looseD26()*2.*0.85; }

double looseD23(){ return looseB23(); }

double looseC16(){ return TMath::Pi() * gc_inner_r*gc_inner_r; }

double looseC17(){ return TMath::Pi()*2.*gc_inner_r*gc_inner_h*2.; }

double looseC18(){ return (0.9*looseC17()+looseC16())/(2.*looseC16()+looseC17()); }

double looseD20(){ return looseC18()*mass_halfimmersed_acrylic*(oxygenmass*2./(oxygenmass*2.+carbonmass*5.+8.*hydrogenmass)); }

double looseD21(){ return looseD20()/((looseC10()+looseC11())*12./14.1); }

double looseD22(){ return looseD21()*probrat; }

double looseD24(){ return looseD22()*looseD23(); }

double looseD25(){ return looseD24()*0.4; }

double looseD26(){ return looseD25()*0.75; }

double looseB20(){ return mass_immersed_acrylic*(oxygenmass*2./(oxygenmass*2+carbonmass*5.+8.*hydrogenmass)); } 

double looseB21(){ return looseB20()/((looseC10()+looseC11())*12./14.1); }

double looseB22(){ return looseB21()*probrat; }

double looseB24(){ return looseB22()*looseB23(); }

double looseB25(){ return looseB24(); }

double looseI18(){ return looseB23()-looseI13(); }

double looseD8(){  return gc_inner_h-nt_inner_h; }

double looseC8(){  return gc_inner_r-nt_inner_r; }

double looseB11() { return ((nt_inner_r+looseC8())*(nt_inner_r+looseC8())*(nt_inner_h+looseD8())*2*TMath::Pi()
             + 2./3. * (nt_inner_r+looseC8())*lidslope*(nt_inner_r+looseC8())*(nt_inner_r+looseC8())*TMath::Pi())
             /1000000.-looseB10(); }

double looseC11() { return looseD5()*looseB11(); }

double looseE3()  { return 120.-kg_oxygen_ppo_nt; }

double looseI3()  { return looseE3()*looseC3()*1000.; }

double looseI16() { return looseI3()/(looseC11() * 12./14.1)/1000. * carbonmass/oxygenmass; }

double looseI17() { return looseI16()*probrat; }

double looseB23() { return 358.8; }

double looseI13() { return looseB23()*n_c12cap_target/n_c12cap; }

double looseC10() { return looseB10()*looseD5(); }

double looseC4() { return oxygenmass/(oxygenmass + 8*hydrogenmass + 4.*carbonmass); }

double looseD6() { return 0.5; }

double looseD5() { return 0.804; }

double looseD7() { return looseD5()*looseD6()/100.*1000.; }

double looseF4() { return looseD7()*looseC4(); }

double looseH4() { return looseF4()*looseB10(); }

double looseC3() { return oxygenmass/(oxygenmass + nitrogenmass + 11.*hydrogenmass+15*carbonmass); }

double looseH3() { return kg_oxygen_ppo_nt*looseC3()*1000.; }

double looseB10(){ return (nt_inner_r*nt_inner_r*nt_inner_h*2*TMath::Pi()
         + 2./3. * nt_inner_r*0.03*nt_inner_r*nt_inner_r*TMath::Pi())/1000000; }

double looseF2(){ return 0.99 * 6.*oxygenmass/157.25;}

double looseH2(){ return looseF2()*looseB10(); }

double looseH5(){ return looseH2() + looseH3() + looseH4(); }

double looseI10(){ return looseH5()/(looseC10() * 12./14.1) * carbonmass/oxygenmass / 1000.; }

double looseI12(){ return looseI10()*probrat; }

double looseI14(){ return looseI12() * looseI13(); }

double looseI19(){ return looseI17() * looseI18(); }

double looseJ14(){ return looseI14() * 2.; }

double looseJ19(){ return looseI19()*2.; }

double looseB26(){ return looseB25()*0.85; }

double looseD27(){ return looseD26()*0.43;}

double looseC26(){ return looseB25()*2*0.95;}

double looseE27(){ return looseE26()*0.47;}

double looseG49(){ return looseI14() + looseI19() + looseB26() + looseD27(); }

double looseH49(){ return looseJ14() + looseJ19() + looseC26() + looseE27(); }

double looseI49(){ return (looseG49() + looseH49())/2.; }

double looseI50(){ return (looseG50()+looseH50())/2; }

double looseE28(){ return looseE27()*0.47; }

double looseH50(){ return looseJ14() + looseJ19() + looseC26() + looseE28(); }

double looseD28(){ return looseD27()*0.43; }

double looseG50() { return looseI14() + looseI19() + looseB26() + looseD28(); }

void dcfluids_finalfit()
{
  const double n_o16cap_beta = looseI49();
  const double n_o16cap_betan = looseI50();

  printf("TECHNOTE 5.3: Gaussian central value of number of effective "
    "beta O-16 captures per day: %.1f\n", n_o16cap_beta);
  printf("TECHNOTE 5.3: Gaussian central value of number of effective "
    "beta-n O-16 captures per day: %.1f\n", n_o16cap_betan);

  printf("const double n_o16cap_beta  = %f;\n", n_o16cap_beta);
  printf("const double n_o16cap_betan = %f;\n", n_o16cap_betan);

  
}
