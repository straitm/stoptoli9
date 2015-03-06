#include "consts.h"
#include "deadtime.C"

void deadtime_finalfit()
{
  TFile * f = new TFile(rootfile3up /* no, should not be muinfo */, "read");
  TTree * t = (TTree *)f->Get("t");

  t->Draw("eff(fq, mx, my, mz, 0) >> h(1000, 0, 1)", "ndecay == 0 && !(abs(mz) < 1233 + 0.03*(1154 - sqrt(mx*mx+my*my)) && mx**2+my**2 < 1154**2) && dt < 30 && e < 15 && miche < 12 && !earlymich && e > 4");
  h=h;

  printf("%sGC efficiency: %.2f%s\n", RED, h->GetMean()*100, CLR);

  t->Draw("eff(fq, mx, my, mz, 0) >> h", "ndecay == 0 && (abs(mz) < 1233 + 0.03*(1154 - sqrt(mx*mx+my*my)) && mx**2+my**2 < 1154**2) && dt < 30 && e < 15 && miche < 12 && !earlymich && e > 4");

  printf("%sTarget efficiency: %.2f%s\n", RED, h->GetMean()*100, CLR);

  t->Draw("eff(fq, mx, my, mz, 0) >> h", "ndecay == 0 && !(abs(mz) < 1068 && mx**2+my**2 < 1068**2) && (abs(mz) < 1233 + 0.03*(1154 - sqrt(mx*mx+my*my)) && mx**2+my**2 < 1154**2) && dt < 30 && e < 15 && miche < 12 && !earlymich && e > 4");

  printf("%sHe-6 T4 (outermost) efficiency: %.2f%s\n", RED, h->GetMean()*100, CLR);

  t->Draw("eff(fq, mx, my, mz, 0) >> h", "ndecay == 0 && !(abs(mz) < 933 && mx**2+my**2 < 933**2) && (abs(mz) < 1068 && mx**2+my**2 < 1068**2) && dt < 30 && e < 15 && miche < 12 && !earlymich && e > 4");

  printf("%sHe-6 T3 efficiency: %.2f%s\n", RED, h->GetMean()*100, CLR);

  t->Draw("eff(fq, mx, my, mz, 0) >> h", "ndecay == 0 && !(abs(mz) < 740 && mx**2+my**2 < 740**2) && (abs(mz) < 933 && mx**2+my**2 < 933**2) && dt < 30 && e < 15 && miche < 12 && !earlymich && e > 4");

  printf("%sHe-6 T2 efficiency: %.2f%s\n", RED, h->GetMean()*100, CLR);

  t->Draw("eff(fq, mx, my, mz, 0) >> h", "ndecay == 0 && (abs(mz) < 740 && mx**2+my**2 < 740**2) && dt < 30 && e < 15 && miche < 12 && !earlymich && e > 4");

  printf("%sHe-6 T1 (innermost) efficiency: %.2f%s\n", RED, h->GetMean()*100, CLR);

  printf("These are ONLY the dt efficiencies.  Don't forget the dr efficiencies!\n");
}
