#include "consts.h"
#include "deadtime.C"

unsigned int factorial(const unsigned int n)
{
  if(n < 2) return 1;
  unsigned int a = 2;
  for(int m = 3; m <= n; m++) a*=m;
  return a;
}

void doit(TTree * t, const int ntrue, const int nseen, const int early)
{
  printf("\nEfficiencies for seeing %d neutron%s out of %d, %s early ones:\n", nseen,
          nseen == 1?"":"s", ntrue, early?"INCLUDING":"excluding");

  const unsigned int comb = factorial(ntrue)/factorial(nseen)/factorial(ntrue-nseen);

  t->Draw(Form("%d * (1-eff(fq,mx,my,mz,%d))**%d * eff(fq, mx, my, mz, %d)**%d >> h(1000,0,1)", comb, early, ntrue-nseen, early, nseen), "");
  h=h;
  printf("%sOverall efficiency: %.2f%s\n", RED, h->GetMean()*100, CLR);

  const char * drawstring = "%d * (1-eff(fq,mx,my,mz,%d))**%d * eff(fq, mx, my, mz, %d)**%d >> h";

  t->Draw(Form(drawstring, comb, early, ntrue-nseen, early, nseen),
    "(abs(mz) < 1229 + 0.03*(1150 - sqrt(mx*mx+my*my)) && mx**2+my**2 < 1150**2)");
  printf("%sTarget/GC efficiency: %.2f, ", RED, h->GetMean()*100);

  t->Draw(Form(drawstring, comb, early, ntrue-nseen, early, nseen),
    "!(abs(mz) < 1237 + 0.03*(1158 - sqrt(mx*mx+my*my)) && mx**2+my**2 < 1158**2)");
  printf("%.2f%s\n", h->GetMean()*100, CLR);

  const double spillout = 0.41;
  
  t->Draw(Form("%d * (1-eff(fq,mx,my,mz,%d)*%f)**%d * (eff(fq, mx, my, mz, %d)*%f)**%d >> h",
               comb, early, spillout, ntrue-nseen, early, spillout, nseen),
    "!(abs(mz) < 1586 + 0.03*(1508 - sqrt(mx*mx+my*my)) && mx**2+my**2 < 1508**2)");
  printf("%sGC vessel efficiency: %.2f%s\n", RED, h->GetMean()*100, CLR);

  // Target acrylic region defined as T4 plus 200mm of GC.  Force events to be 
  // in the acrylic for the deadtime calculation
  t->Draw(Form("%d * (1-eff(fq,1154,0,0,%d,0))**%d * eff(fq,1154,0,0,%d,0)**%d >> h",
               comb, early, ntrue-nseen, early, nseen),
    "(abs(mz) < 1429 + 0.03*(1350 - sqrt(mx*mx+my*my)) && mx**2+my**2 < 1350**2) && "
    "!(abs(mz) < 1068 && mx**2+my**2 < 1068**2) && (abs(mz) < 1233 + 0.03*(1154 - sqrt(mx*mx+my*my)) && mx**2+my**2 < 1154**2)");
  printf("%sTarget vessel efficiency (no early H-n): %.2f%s\n", RED, h->GetMean()*100, CLR);

  t->Draw(Form(drawstring, comb, early, ntrue-nseen, early, nseen),
    "!(abs(mz) < 1068 && mx**2+my**2 < 1068**2) && (abs(mz) < 1229 + 0.03*(1150 - sqrt(mx*mx+my*my)) && mx**2+my**2 < 1150**2)");
  printf("%sHe-6 T4 (outermost)/T3/T2/T1 efficiencies: %.2f, ", RED, h->GetMean()*100);

  t->Draw(Form(drawstring, comb, early, ntrue-nseen, early, nseen),
    "!(abs(mz) < 933 && mx**2+my**2 < 933**2) && (abs(mz) < 1068 && mx**2+my**2 < 1068**2)");
  printf("%.2f, ", h->GetMean()*100);

  t->Draw(Form(drawstring, comb, early, ntrue-nseen, early, nseen),
    "!(abs(mz) < 740 && mx**2+my**2 < 740**2) && (abs(mz) < 933 && mx**2+my**2 < 933**2)");
  printf("%.2f, ", h->GetMean()*100);

  t->Draw(Form(drawstring, comb, early, ntrue-nseen, early, nseen),
    "(abs(mz) < 740 && mx**2+my**2 < 740**2)");
  printf("%.2f%s\n", h->GetMean()*100, CLR);

  t->Draw(Form(drawstring, comb, early, ntrue-nseen, early, nseen),
    "(mz > -1175 && mx**2+my**2 < 1050**2 && chi2 < 2 && abs(fez + 62*ivdedx/2 - 8847.2) < 1000)");
  printf("%sHigh-purity cuts: %.2f%s\n", RED, h->GetMean()*100, CLR);

/*  t->Draw(Form(drawstring, comb, early, ntrue-nseen, early, nseen), "fq < 1.8e6");
  printf("%sEfficiency when fq < 1.8e6: %.2f%s\n", RED, h->GetMean()*100, CLR); */

  puts("");
}

void deadtime_finalfit(const int ntrue)
{
  gErrorIgnoreLevel = kError;
  TFile * f = new TFile(//rootfile3up, "read");
     "/cp/s4/strait/fullfido-300s-3-25MeV-20150219-b12deadtimepass.root", "read");
  TTree * t = (TTree *)f->Get("t");

/*  TFile * temp = new TFile("/tmp/tmp.root", "recreate");
  TTree * sel = t->CopyTree("ndecay == 0 && dt < 30 && e < 15 && miche < 12 && !earlymich && e > 4");
  sel->Write();
  temp->Close();
   return; */

  printf("\n\n\n");
  for(int e = 0; e < 2; e++) for(int i = 1; i <= ntrue; i++) doit(t, ntrue, i, e);

  printf("These are ONLY the dt efficiencies.  Don't forget the dr efficiencies!\n");
}
