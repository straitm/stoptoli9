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
const double gdmass = 157.25;

const double mulife = 2196.9811e-6;

const double lifetime_c12 = 2028.e-6;
const double lifetime_c13 = 2037.e-6;

const double lifetime_c = lifetime_c12*(1-f13)+lifetime_c13*f13;
const double capprob_c = 1-lifetime_c/mulife;

const double lifetime_o16 = 1796.e-6;
const double capprob_o16 = 1-lifetime_o16/mulife;

const double probrat = capprob_o16/capprob_c;


double gc_vessel_betan_count_high_guess_with_mu_selection_and_beta_energy_and_direction_efficiencies();
double looseH50();
double gc_vessel_betan_count_low_guess_with_mu_selection_and_beta_energy_and_direction_efficiencies();
double looseI50();
double looseG50();
double gc_vessel_beta_count_high_guess_with_mu_selection_and_beta_energy_efficiencies();
double gc_base_area();
double gc_wall_area();
double gc_useful_vessel_fraction();
double oxygen_mass_halfimmersed();
double halfimmersed_acrylic_o_over_c_number_ratio();
double halfimmersed_acrylic_o_over_all_scint_c_capture_ratio();
double halfimmersed_acrlyic_preefficiency_captures_low_guess();
double gc_vessel_beta_count_low_guess_with_mu_selection();
double gc_vessel_beta_count_low_guess_with_mu_selection_and_beta_energy_efficiencies();
double oxygen_mass_immersed_in_acrylic();
double immersed_acrylic_o_over_all_scint_c_number_ratio();
double immersed_acrylic_o_over_all_scint_c_capture_ratio();
double immersed_acrlyic_preefficiency_captures_low_guess();
double n_c12capgc();
double liters_in_gc();
double kg_in_gc();
double kg_ppo_mass_gc();
double looseI3();
double looseI16();
double looseI17();
double kg_in_target();
double looseC4();
double looseD7();
double looseF4();
double looseH4();
double looseC3();
double looseH3();
double liters_in_target();
double looseH2();
double looseH5();
double looseI10();
double looseI12();
double looseI14();
double looseI19();
double looseJ14();
double looseJ19();
double looseB26();
double gc_vessel_beta_count_low_guess_with_mu_selection_and_beta_energy_and_direction_efficiencies();
double looseC26();
double gc_vessel_beta_count_high_guess_with_mu_selection_and_beta_energy_and_direction_efficiencies();
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
const double scint_density = 0.804;
const double THF_fraction_by_weight = 0.005;
const double scint_carbon_mass_fraction = 12/(12 + 1.92);
const double g_per_l_gd = 0.99;
const double g_per_l_oxygen_in_target = g_per_l_gd * 6.*oxygenmass/gdmass;
const double kg_ppo_mass_total = 120.;
const double muon_selection_efficiency_at_gc_vessel = 0.4;

const double immersed_acrlyic_low_guess_energy_efficiency = 0.85,
             immersed_acrlyic_high_guess_energy_efficiency = 0.95;
const double halfimmersed_acrlyic_low_guess_energy_efficiency = 0.75,
             halfimmersed_acrlyic_high_guess_energy_efficiency = 0.85;
const double halfimmersed_acrlyic_low_guess_dir_efficiency = 0.43,
             halfimmersed_acrlyic_high_guess_dir_efficiency = 0.47;
const double portion_of_useful_gc_wall = 0.9;

double gc_vessel_beta_count_high_guess_with_mu_selection_and_beta_energy_efficiencies(){ return gc_vessel_beta_count_low_guess_with_mu_selection_and_beta_energy_efficiencies()*2.*halfimmersed_acrlyic_high_guess_energy_efficiency; }

double gc_base_area(){ return TMath::Pi() * gc_inner_r*gc_inner_r; }

double gc_wall_area(){ return TMath::Pi()*2.*gc_inner_r*gc_inner_h*2.; }

double gc_useful_vessel_fraction(){ return (portion_of_useful_gc_wall*gc_wall_area()+gc_base_area())/(2.*gc_base_area()+gc_wall_area()); }

double oxygen_mass_halfimmersed()
{
  return gc_useful_vessel_fraction()*mass_halfimmersed_acrylic*
     (oxygenmass*2./(oxygenmass*2.+carbonmass*5.+8.*hydrogenmass));
}

double halfimmersed_acrylic_o_over_c_number_ratio()
{
  return oxygen_mass_halfimmersed()/ // XXX missing carbonmass/oxygenmass, right?
    ((kg_in_target()+kg_in_gc())*scint_carbon_mass_fraction);
}

double halfimmersed_acrylic_o_over_all_scint_c_capture_ratio(){ return halfimmersed_acrylic_o_over_c_number_ratio()*probrat; }

double halfimmersed_acrlyic_preefficiency_captures_low_guess(){ return halfimmersed_acrylic_o_over_all_scint_c_capture_ratio()*n_c12cap; }

double gc_vessel_beta_count_low_guess_with_mu_selection(){ return halfimmersed_acrlyic_preefficiency_captures_low_guess()*muon_selection_efficiency_at_gc_vessel; }

double gc_vessel_beta_count_low_guess_with_mu_selection_and_beta_energy_efficiencies(){ return gc_vessel_beta_count_low_guess_with_mu_selection()*halfimmersed_acrlyic_low_guess_energy_efficiency; }

double oxygen_mass_immersed_in_acrylic()
{
  return mass_immersed_acrylic*
    (oxygenmass*2./(oxygenmass*2+carbonmass*5.+8.*hydrogenmass));
} 

double immersed_acrylic_o_over_c_number_ratio()
{
  return oxygen_mass_immersed_in_acrylic() // XXX missing carbonmass/oxygenmass, right?
    /((kg_in_target()+kg_in_gc())*scint_carbon_mass_fraction);
}

double immersed_acrylic_o_over_all_scint_c_capture_ratio(){
  return immersed_acrylic_o_over_c_number_ratio()*probrat; }

double immersed_acrlyic_preefficiency_captures_low_guess(){
  return immersed_acrylic_o_over_all_scint_c_capture_ratio()*n_c12cap; }

double n_c12capgc(){ return n_c12cap-n_c12captarget; }

double liters_in_gc() {
  return (gc_inner_r*gc_inner_r*gc_inner_h*2*TMath::Pi()
   + 2./3. * gc_inner_r*lidslope*gc_inner_r*gc_inner_r*TMath::Pi())
   /1000000.-liters_in_target();
}

double kg_in_gc() { return scint_density*liters_in_gc(); }

double kg_ppo_mass_gc()  { return kg_ppo_mass_total-kg_oxygen_ppo_nt; }

double looseI3()  { return kg_ppo_mass_gc()*looseC3()*1000.; }

double looseI16() { return looseI3()/(kg_in_gc() * scint_carbon_mass_fraction)/1000. * carbonmass/oxygenmass; }

double looseI17() { return looseI16()*probrat; }

double kg_in_target() { return liters_in_target()*scint_density; }

double looseC4() { return oxygenmass/(oxygenmass + 8*hydrogenmass + 4.*carbonmass); }

double looseD7() { return scint_density*THF_fraction_by_weight*1000.; }

double looseF4() { return looseD7()*looseC4(); }

double looseH4() { return looseF4()*liters_in_target(); }

double looseC3() { return oxygenmass/(oxygenmass + nitrogenmass + 11.*hydrogenmass+15*carbonmass); }

double looseH3() { return kg_oxygen_ppo_nt*looseC3()*1000.; }

double liters_in_target(){ return (nt_inner_r*nt_inner_r*nt_inner_h*2*TMath::Pi()
         + 2./3. * nt_inner_r*lidslope*nt_inner_r*nt_inner_r*TMath::Pi())/1000000; }

double looseH2(){ return g_per_l_oxygen_in_target*liters_in_target(); }

double looseH5(){ return looseH2() + looseH3() + looseH4(); }

double looseI10(){ return looseH5()/(kg_in_target() * scint_carbon_mass_fraction) * carbonmass/oxygenmass / 1000.; }

double looseI12(){ return looseI10()*probrat; }

double looseI14(){ return looseI12() * n_c12captarget; }

double looseI19(){ return looseI17() * n_c12capgc(); }

double looseJ14(){ return looseI14() * 2.; }

double looseJ19(){ return looseI19()*2.; }

double looseB26(){ return immersed_acrlyic_preefficiency_captures_low_guess()*immersed_acrlyic_low_guess_energy_efficiency; }

double gc_vessel_beta_count_low_guess_with_mu_selection_and_beta_energy_and_direction_efficiencies(){ return gc_vessel_beta_count_low_guess_with_mu_selection_and_beta_energy_efficiencies()*halfimmersed_acrlyic_low_guess_dir_efficiency;}

double looseC26(){ return immersed_acrlyic_preefficiency_captures_low_guess()*2*immersed_acrlyic_high_guess_energy_efficiency;}

double gc_vessel_beta_count_high_guess_with_mu_selection_and_beta_energy_and_direction_efficiencies(){ return gc_vessel_beta_count_high_guess_with_mu_selection_and_beta_energy_efficiencies()*halfimmersed_acrlyic_high_guess_dir_efficiency;}

double looseG49(){ return looseI14() + looseI19() + looseB26() + gc_vessel_beta_count_low_guess_with_mu_selection_and_beta_energy_and_direction_efficiencies(); }

double looseH49(){ return looseJ14() + looseJ19() + looseC26() + gc_vessel_beta_count_high_guess_with_mu_selection_and_beta_energy_and_direction_efficiencies(); }

double looseI49(){ return (looseG49() + looseH49())/2.; }

double looseI50(){ return (looseG50()+looseH50())/2; }

double gc_vessel_betan_count_high_guess_with_mu_selection_and_beta_energy_and_direction_efficiencies(){ return gc_vessel_beta_count_high_guess_with_mu_selection_and_beta_energy_and_direction_efficiencies()*halfimmersed_acrlyic_high_guess_dir_efficiency; }

double looseH50(){ return looseJ14() + looseJ19() + looseC26() + gc_vessel_betan_count_high_guess_with_mu_selection_and_beta_energy_and_direction_efficiencies(); }

double gc_vessel_betan_count_low_guess_with_mu_selection_and_beta_energy_and_direction_efficiencies(){ return gc_vessel_beta_count_low_guess_with_mu_selection_and_beta_energy_and_direction_efficiencies()*halfimmersed_acrlyic_low_guess_dir_efficiency; }

double looseG50() { return looseI14() + looseI19() + looseB26() + gc_vessel_betan_count_low_guess_with_mu_selection_and_beta_energy_and_direction_efficiencies(); }

void dcfluids_finalfit()
{
  const double n_o16cap_beta = looseI49();
  const double n_o16cap_betan = looseI50();

  printf("TECHNOTE 5.3: Gaussian central value of number of effective "
    "beta O-16 captures per day: %.1f\n", n_o16cap_beta);
  printf("TECHNOTE 5.3: Gaussian central value of number of effective "
    "beta-n O-16 captures per day: %.1f\n", n_o16cap_betan);

  printf("const double n_o16cap_beta  = %f; ", n_o16cap_beta);
  printf("const double n_o16cap_betan = %f;\n", n_o16cap_betan);

  
}
