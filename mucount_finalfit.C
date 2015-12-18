#include <stdio.h>
#include "musicchargeratio_finalfit_far.out.h"
#include "musicchargeratio_finalfit_near.out.h"

/* Meant to be run #included in other files */

ve mucountfinalfit_cut(const char * const cut, const bool far)
{
  printf("Counting muons...\n");
  TFile *_file0 = TFile::Open(far?rootfile3up:rootfile3up_near, "read");
  TTree * t = (TTree *)_file0->Get("t");

  const int rawcount = t->GetEntries(cut);
  printf("Raw count: %d\n", rawcount);

  if(!far)
    fprintf(stderr,"WARNING! The ND muon contamination numbers "
                   "are dummy values\n");

  const double mumf = far? mum_frac: mum_frac_near;
  const double mumf_e = far? mum_frac_err: mum_frac_near_err;
  const double con = far? mum_contamination:
                          mum_contamination_near;
  const double con_e = far? mum_contamination_err:
                            mum_contamination_err_near;

  ve answer;
  answer.val = rawcount*mumf*(1-con);
  answer.err = sqrt(pow(rawcount*con_e,2) +pow(rawcount*mumf*con_e,2));
  printf("Translated to mu- & corrected for contamination: %f +- %f\n",
    answer.val, answer.err);

  delete t;
  delete _file0;

  return answer;
}
