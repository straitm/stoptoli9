#include <math.h>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include "TTree.h"
#include "TError.h"
#include "TFile.h"
#include "TH1.h"
#include <string>
using std::string;

#include "neff.C"

unsigned int factorial(const unsigned int n)
{
  if(n < 2) return 1;
  unsigned int a = 2;
  for(unsigned int m = 3; m <= n; m++) a*=m;
  return a;
}

void reallydoit(TTree * t, const char * const cut,
                const char * const varname, const char * const desc,
                const int ntrue,
                const bool hasdeadtime // i.e. No deadtime assumed for
                                       // FD after FE upgrade, or ND
                                       // anytime. This assumption still
                                       // needs to be checked, however!
)
{
  // Only count each muon once, and verify that it is a stopping muon by
  // selecting michel decays
  TTree * selt = t->CopyTree(Form("%s && ndecay == 0 && miche > 12", cut));

  float fq, fqiv, mx, my, mz;
  selt->SetBranchStatus("*", 0);
  selt->SetBranchStatus("fq", 1);
  selt->SetBranchStatus("fqiv", 1);
  selt->SetBranchStatus("mx", 1);
  selt->SetBranchStatus("my", 1);
  selt->SetBranchStatus("mz", 1);
  selt->SetBranchAddress("fq", &fq);
  selt->SetBranchAddress("fqiv", &fqiv);
  selt->SetBranchAddress("mx", &mx);
  selt->SetBranchAddress("my", &my);
  selt->SetBranchAddress("mz", &mz);

  for(int nseen = 1; nseen <= ntrue; nseen++){
    const unsigned int comb = factorial(ntrue)/factorial(nseen)/factorial(ntrue-nseen);

    for(int early = 0; early < 2; early++){
      double effsum = 0;
      for(int i = 0; i < selt->GetEntries(); i++){
        selt->GetEntry(i);

        // energy only matters if it causes deadtime, otherwise zero it
        // to get the deadtime-free answer.
        const double dteff = neff_dt(hasdeadtime*fq, fqiv, mx, my, mz, early);

        // When there is only one neutron, we can separate out the dr
        // efficiency and handle it elsewhere, which is probably a
        // good thing. But when there is more than one, and you want
        // an overall average, I think that it has to be done here,
        // before the average is taken. I only need multiple neutron
        // efficiencies for the dr800 definition, since I only use the
        // dr1000 definition for betan, which (out of isotopes I look
        // for) only ever has one neutron.
        if(ntrue == 1){
          effsum += comb * pow(1 - dteff, ntrue - nseen) * pow(dteff, nseen);
        }
        else{
          const double dreff = neff_dr_800(mx, my, mz);
          effsum += comb * pow(1 - dteff*dreff, ntrue - nseen)
                         * pow(    dteff*dreff, nseen);
        }
      }
      printf("%s dt %sefficiency for seeing %d neutrons of %d, %s: %.2f%%\n",
             desc,
             ntrue>1?"and dr800 ":"",
             nseen,
             ntrue,
             early?"INCLUDING early ones":"excluding early ones",
             effsum/selt->GetEntries()*100);

      char xofx[7]; // works up to 99
      sprintf(xofx, "%dof%d", nseen, ntrue);

      printf("const double n%seff_dt%s%s%s%s = %.16f;\n",
             ntrue>1?xofx:"",
             ntrue>1?"_dr_800":"",
             varname,
             early?"_wearly":"",
             hasdeadtime?"":"_nodeadtime",
             effsum/selt->GetEntries());
    }
  }
  delete selt;
}

void doit(TTree * t, const int ntrue, const bool hasdeadtime)
{
  printf("Opening temp file...\n"); fflush(stdout);
  char filename[100];
  strcpy(filename, "/tmp/b12like.XXXXXX");
  close(mkstemp(filename));
  TFile * tmpfile = new TFile(filename, "recreate", "", 0);
  if(!tmpfile || tmpfile->IsZombie()){
    fprintf(stderr, "Could not open temp file %s\n", filename);
    exit(1);
  }

  const string t0_cut = "(mx**2+my**2 < 740**2 && abs(mz) < 740)";
  const string t1_cut = "(mx**2+my**2 < 933**2 && abs(mz) < 933)";
  const string t2_cut = "(mx**2+my**2 < 1068**2 && abs(mz) < 1068)";
  const string target_cut =
    "(mx**2+my**2 < 1154**2 && "
    "abs(mz) < 1229 + 0.03*(1154 - sqrt(mx**2+my**2)))";

  reallydoit(t, "1", "", "Overall", ntrue, hasdeadtime); 

  reallydoit(t, target_cut.c_str(), "_targ", "Target", ntrue, hasdeadtime); 
  reallydoit(t, ("!"+target_cut).c_str(), "_gc", "Gamma Catcher", ntrue, hasdeadtime); 
    

  const double spillout = 0.41;

  const double b12gammatargfrac = 
    double(t->GetEntries(("ndecay == 0 && miche > 12 && fq/8300 < 215 && " + target_cut).c_str()))/
           t->GetEntries( "ndecay == 0 && miche > 12 && fq/8300 < 215");

  printf("B-12gamma selected muons (i.e. under 215) that stop in the target: %.1f\n",
         b12gammatargfrac*100);
  if(ntrue == 1 && hasdeadtime) // Need this to only go in to the
                                // .out.h file once
    printf("const double b12gamma215targfraction = %f;\n", b12gammatargfrac);

  reallydoit(t, "fq/8300 < 215", "_215MeV", "E_mu < 215 MeV (B-12 gamma search)", ntrue, hasdeadtime);

  reallydoit(t, t0_cut.c_str(),                        "_t0", "He-6 T0 region", ntrue, hasdeadtime);
  reallydoit(t, (t1_cut     + "&&!" + t0_cut).c_str(), "_t1", "He-6 T1 region", ntrue, hasdeadtime);
  reallydoit(t, (t2_cut     + "&&!" + t1_cut).c_str(), "_t2", "He-6 T2 region", ntrue, hasdeadtime);
  reallydoit(t, (target_cut + "&&!" + t2_cut).c_str(), "_t3", "He-6 T3 region", ntrue, hasdeadtime);

}

void neff_dt_finalfit(const int ntrue)
{
  TFile * f = new TFile(rootfilemuinfo, "read");
  TTree * t = (TTree *)f->Get("t");

  string myownsource;
  string tmp;
  ifstream infile("neff_dt_finalfit.C");
  if(!infile.is_open()){
    fprintf(stderr, "Could not open deadtime_finalfit.C\n");
    exit(1);
  }
  while(infile >> tmp) myownsource += tmp;
  infile.close();

  for(int i = 0; i < t->GetListOfBranches()->GetEntries(); i++){
    const char * name = t->GetListOfBranches()->At(i)->GetName();
    if(strstr(myownsource.c_str(), name) == NULL)
      t->SetBranchStatus(name, 0);
  }

  doit(t, ntrue, true);
  doit(t, ntrue, false);

  printf("These are ONLY the dt efficiencies.  Don't forget the dr efficiencies!\n");
}
