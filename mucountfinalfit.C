/* Meant to be run #included in other files */

double mucountfinalfit_cut(const char * const cut,
                           const bool returnerr = false)
{
  TFile *_file0 = TFile::Open(rootfile3up, "read");
  TTree * t = (TTree *)_file0->Get("t");

  const int rawcount = t->GetEntries(cut);
  printf("Raw count: %d\n", rawcount);

  const double val = rawcount*mum_frac*(1-mum_contamination);
  const double err = sqrt(pow(rawcount*mum_frac_err,2)
                         +pow(rawcount*mum_frac*mum_contamination_err,2));
  printf("Translated to mu- & corrected for contamination: %f +- %f\n",
    val, err);

  return returnerr? err: val;
}

void mucountfinalfit()
{
  mucountfinalfit_cut(
    "ndecay == 0 && mx**2+my**2 < 1050**2 && mz > -1175 && "
    "abs(fez + 62*ivdedx/2 - 8847.2) < 1000 && rchi2 < 2");

  puts("LESSPOS:");

  mucountfinalfit_cut(
    "ndecay == 0 && mx**2+my**2 < 900**2 && mz > -900 && "
    "abs(fez + 62*ivdedx/2 - 8847.2) < 1000 && rchi2 < 2");

  puts("LESSSLANT:");

  mucountfinalfit_cut(
    "ndecay == 0 && mx**2+my**2 < 1050**2 && mz > -1175 && "
    "abs(fez + 62*ivdedx/2 - 8847.2) < 600 && rchi2 < 2");

  puts("LESSCHI2:");

  mucountfinalfit_cut(
    "ndecay == 0 && mx**2+my**2 < 1050**2 && mz > -1175 && "
    "abs(fez + 62*ivdedx/2 - 8847.2) < 1000 && rchi2 < 1.25");
}
