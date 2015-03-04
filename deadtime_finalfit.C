#include "consts.h"
#include "deadtime.C"

void deadtime_finalfit()
{
  TFile * f = new TFile(rootfile3up, "read");
  TTree * t = (TTree *)f->Get("t");

  t->Draw("eff(fq, mx, my, mz, 0) >> h(1000, 0, 1)", "ndecay == 0 && !(abs(mz) < 1233 + 0.03*(1154 - sqrt(mx*mx+my*my)) && mx**2+my**2 < 1154**2) && dt < 30 && e < 15 && miche < 12 && !earlymich && e > 4");
  h=h;

  printf("%sGC efficiency: %.2f%s\n", RED, h->GetMean()*100, CLR);

  t->Draw("eff(fq, mx, my, mz, 0) >> h", "ndecay == 0 && (abs(mz) < 1233 + 0.03*(1154 - sqrt(mx*mx+my*my)) && mx**2+my**2 < 1154**2) && dt < 30 && e < 15 && miche < 12 && !earlymich && e > 4");

  printf("%sTarget efficiency: %.2f%s\n", RED, h->GetMean()*100, CLR);
}
