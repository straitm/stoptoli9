#include "consts.h"
#include "TFile.h"
#include "TTree.h"
#include "stdio.h"

/*
 * Prints the message once with the requested precision and in RED, then
 * again with all digits in the default color, starting with the first
 * number.
 *
 * The message must only have %f substitutions.
 */
static void printtwice(const char * const msg, const int digits, ...)
{
  char * bmsg = (char *)malloc(strlen(msg)+100); // Ha!
  char * pmsg = (char *)malloc(strlen(msg)+100); // Ha!
  
  // Just for fun...
  char * pmp = pmsg;
  char * bmp = bmsg;
  bool gotone = false;
  for(unsigned int i = 0; i <= strlen(msg); i++){
    switch(msg[i]){
      case '\0':
        *pmp++ = '\0';
        *bmp++ = '\0';
        break;
      case '%':
        gotone = true;
        *pmp++ = '%';
        *bmp++ = '%';
        *pmp++ = '.';
        *pmp++ = digits+'0';
        break;
      default:
        *pmp++ = msg[i];
        if(gotone) *bmp++ = msg[i];
        break;
    }
  }
  
  va_list ap;
  va_start(ap, digits);
  printf(RED);
  vprintf(pmsg, ap);
  printf(CLR);

  va_start(ap, digits);
  vprintf(bmsg, ap);
}

double mucountfinalfit_cut(const char * const cut,
                           const bool returnerr = false)
{
  TFile *_file0 = TFile::Open(rootfile3up, "read");
  TTree * t = (TTree *)_file0->Get("t");

  const int rawcount = t->GetEntries(cut);
  printf("Raw count: %d\n", rawcount);

  const double mum_frac = 0.4410,
               mum_frac_err = 0.0032,
               contamination = 0.0028,
               contamination_err = 0.0019;

  const double val = rawcount*mum_frac*(1-contamination);
  const double err = sqrt(pow(rawcount*mum_frac_err,2)
                         +pow(rawcount*mum_frac*contamination_err,2));
  printtwice("Translated to mu- & corrected for contamination: %f +- %f\n", 0,
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
