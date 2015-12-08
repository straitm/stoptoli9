#include <stdio.h>
#include "musicchargeratio_finalfit.out.h"

/* Meant to be run #included in other files */

ve mucountfinalfit_cut(const char * const cut, const bool far)
{
  printf("Counting muons...\n");
  TFile *_file0 = TFile::Open(far?rootfile3up:rootfile3up_near, "read");
  TTree * t = (TTree *)_file0->Get("t");

  const int rawcount = t->GetEntries(cut);
  printf("Raw count: %d\n", rawcount);

  if(!far)
    fprintf(stderr,"WARNING! The mu- fraction and muon contamination\n"
                   "are those for the FD, even though you are looking\n"
                   "at ND data.  The numbers are surely different.\n");

  ve answer;
  answer.val = rawcount*mum_frac*(1-mum_contamination);
  answer.err = sqrt(pow(rawcount*mum_frac_err,2)
                   +pow(rawcount*mum_frac*mum_contamination_err,2));
  printf("Translated to mu- & corrected for contamination: %f +- %f\n",
    answer.val, answer.err);

  delete t;
  delete _file0;

  return answer;
}
