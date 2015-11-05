#include <fstream>
#include <stdio.h>
#include "consts.h"
#include "TTree.h"
#include "TError.h"
#include "TFile.h"
#include "TH1.h"

#include "deadtime.C"

unsigned int factorial(const unsigned int n)
{
  if(n < 2) return 1;
  unsigned int a = 2;
  for(unsigned int m = 3; m <= n; m++) a*=m;
  return a;
}

void doit(TTree * t, const int ntrue, const int nseen, const int early)
{
  const string t0_cut = "(mx**2+my**2 < 740**2 && abs(mz) < 740)";
  const string t1_cut = "(mx**2+my**2 < 933**2 && abs(mz) < 933)";
  const string t2_cut = "(mx**2+my**2 < 1068**2 && abs(mz) < 1068)";
  const char * const target_cut =
    "(mx**2+my**2 < 1154**2 && "
    "abs(mz) < 1229 + 0.03*(1154 - sqrt(mx**2+my**2)))";

  printf("\nEfficiencies for seeing %d neutron%s out of %d, %s early ones:\n", nseen,
          nseen == 1?"":"s", ntrue, early?"INCLUDING":"excluding");

  const unsigned int comb = factorial(ntrue)/factorial(nseen)/factorial(ntrue-nseen);
  const char * eff = early?"eff(fq,mx,my,mz,1)":
                           "eff(fq,mx,my,mz,0)";

  TH1D * h = new TH1D("h", "", 1000,0,1);

  fprintf(stderr,"drawing %lld entries...\n", t->GetEntries());
  
  t->Draw(Form("%d * (1-%s)**%d * %s**%d >> h",
    comb, eff, ntrue-nseen, eff, nseen), "");
  fprintf(stderr, "%sOverall efficiency: %.2f%s\n", RED, h->GetMean()*100, CLR);

  const char * drawstring = "%d * (1-%s)**%d * %s**%d >> h";

  t->Draw(Form(drawstring, comb, eff, ntrue-nseen, eff, nseen), target_cut);
  printf("%sTarget/GC efficiency: %.2f, ", RED, h->GetMean()*100);

  t->Draw(Form(drawstring, comb, eff, ntrue-nseen, eff, nseen),
               ("!" + string(target_cut)).c_str());
  printf("%.2f%s\n", h->GetMean()*100, CLR);

  const double spillout = 0.41;
  
  t->Draw(Form("%d * (1-%s*%f)**%d * (%s*%f)**%d >> h",
               comb, eff, spillout, ntrue-nseen, eff, spillout, nseen),
    "!(abs(mz) < 1586 + 0.03*(1508 - sqrt(mx*mx+my*my)) && mx**2+my**2 < 1508**2)");
  printf("%sGC vessel efficiency: %.2f%s\n", RED, h->GetMean()*100, CLR);

  // Target acrylic region defined as T4 plus 200mm of GC.  Force events to be 
  // in the acrylic for the deadtime calculation
  t->Draw(Form("%d * (1-eff(fq,1154,0,0,%d,0))**%d * eff(fq,1154,0,0,%d,0)**%d >> h",
               comb, early, ntrue-nseen, early, nseen),
    "(abs(mz) < 1429 + 0.03*(1350 - sqrt(mx*mx+my*my)) && mx**2+my**2 < 1350**2) && "
    "!(abs(mz) < 1068 && mx**2+my**2 < 1068**2) && (abs(mz) < 1233 + 0.03*(1154 - sqrt(mx*mx+my*my)) && mx**2+my**2 < 1154**2)");
  printf("%sTarget vessel efficiency (no early H-n): %.2f%s\n", RED, h->GetMean()*100, CLR);

  t->Draw(Form(drawstring, comb, eff, ntrue-nseen, eff, nseen),
    "!(abs(mz) < 1068 && mx**2+my**2 < 1068**2) && (abs(mz) < 1229 + 0.03*(1150 - sqrt(mx*mx+my*my)) && mx**2+my**2 < 1150**2)");
  printf("%sHe-6 T4 (outermost)/T3/T2/T1 efficiencies: %.2f, ", RED, h->GetMean()*100);

  t->Draw(Form(drawstring, comb, eff, ntrue-nseen, eff, nseen),
    "!(abs(mz) < 933 && mx**2+my**2 < 933**2) && (abs(mz) < 1068 && mx**2+my**2 < 1068**2)");
  printf("%.2f, ", h->GetMean()*100);

  t->Draw(Form(drawstring, comb, eff, ntrue-nseen, eff, nseen),
    "!(abs(mz) < 740 && mx**2+my**2 < 740**2) && (abs(mz) < 933 && mx**2+my**2 < 933**2)");
  printf("%.2f, ", h->GetMean()*100);

  t->Draw(Form(drawstring, comb, eff, ntrue-nseen, eff, nseen),
    "(abs(mz) < 740 && mx**2+my**2 < 740**2)");
  printf("%.2f%s\n", h->GetMean()*100, CLR);

  t->Draw(Form(drawstring, comb, eff, ntrue-nseen, eff, nseen),
    "(mz > -1175 && mx**2+my**2 < 1050**2 && rchi2 < 2 && abs(fez + 62*ivdedx/2 - 8847.2) < 1000)");
  printf("%sHigh-purity cuts: %.2f%s\n", RED, h->GetMean()*100, CLR);

  t->Draw(Form(drawstring, comb, eff, ntrue-nseen, eff, nseen),
    "(mz > -1175 && mx**2+my**2 < 1050**2 && rchi2 < 1.25 && abs(fez + 62*ivdedx/2 - 8847.2) < 1000)");
  printf("%sHigh-purity lesschi2 cuts: %.2f%s\n", RED, h->GetMean()*100, CLR);

  t->Draw(Form(drawstring, comb, eff, ntrue-nseen, eff, nseen),
    "(mz > -900 && mx**2+my**2 < 900**2 && rchi2 < 2 && abs(fez + 62*ivdedx/2 - 8847.2) < 1000)");
  printf("%sHigh-purity lesspos cuts: %.2f%s\n", RED, h->GetMean()*100, CLR);

  t->Draw(Form(drawstring, comb, eff, ntrue-nseen, eff, nseen),
    "(mz > -1175 && mx**2+my**2 < 1050**2 && rchi2 < 2 && abs(fez + 62*ivdedx/2 - 8847.2) < 600)");
  printf("%sHigh-purity lessslant cuts: %.2f%s\n", RED, h->GetMean()*100, CLR);

  puts("");
}

void deadtime_finalfit(const int ntrue)
{
  TFile * f = new TFile(rootfile3up, "read");
  TTree * t = (TTree *)f->Get("t");

  string myownsource;
  string tmp;
  ifstream infile("deadtime_finalfit.C");
  if(!infile.is_open()){
    fprintf(stderr, "Could not open deadtime_finalfit.C\n");
    exit(1);
  }
  while(infile >> tmp) myownsource += tmp;

  printf("HERE COMES THE SOURCE\n\n\n");
  fputs(myownsource.c_str(), stdout);
  printf("THERE WENT THE SOURCE\n\n\n");


  {
    bool saidoff = false, saidon = false;
    for(int i = 0; i < t->GetListOfBranches()->GetEntries(); i++){
      const char * name = t->GetListOfBranches()->At(i)->GetName();
      if(strstr(myownsource.c_str(), name) == NULL){
        printf("%s %s", saidoff?"":"\nTurning off branch(es): ", name);
        saidoff = true; saidon = false;
        t->SetBranchStatus(name, 0);
      }
      else{
        printf("%s %s", saidon?"":"\nLeaving on branch(es): ", name);
        saidon = true; saidoff = false;
      }
    }
  }
  puts("");
  infile.close();

  for(int e = 0; e < 2; e++) for(int i = 1; i <= ntrue; i++) doit(t, ntrue, i, e);

  printf("These are ONLY the dt efficiencies.  Don't forget the dr efficiencies!\n");
}
