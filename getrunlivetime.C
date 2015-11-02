#include "time.h"
#include "TTree.h"
#include "TFile.h"
#include "TError.h"
#include <algorithm>
#include <vector>
using std::vector;
#include <fstream>
using std::ifstream;

#include "reactorpowerbin.C"

double getrunlivetime(const int run)
{
  TFile * rf = new TFile(Form("/cp/s4/dchooz/cheetah/prod-08-05_p01_v2/"
    "reduced.Run%07d_Seq010.root", run), "read");
  if(!rf || rf->IsZombie()){
    fprintf(stderr, "Could not open reduced file for run %d\n", run);
    exit(1);
  }

  TTree * const data = (TTree* const)rf->Get("data");
  if(!data){
    fprintf(stderr, "Could not get 'data' tree for run %d\n", run);
    exit(1);
  }

  double trgtime;
  data->SetBranchStatus("*", 0); // This really helps!
  data->SetBranchStatus("trgtime", 1);
  data->SetBranchAddress("trgtime", &trgtime);
  const Long64_t nentries = data->GetEntriesFast();
  data->GetEntry(nentries-1);

  delete data;
  rf->Close();

  // Time from the first event to the last event -- the first event 
  // always has trgtime zero.
  const double obviousseconds = trgtime/1e9;

  // On average, the run was live for one more averge interval between
  // events, half at the beginning and half at the end.
  const double endseconds = obviousseconds/(nentries-1);

  const double livetime = obviousseconds+endseconds;

  printf("%d %.9lf\n", run, livetime);

  return livetime;
}

int totallivetime_finalfit()
{
  ifstream infile("runlist");
  if(!infile.is_open()){
    fprintf(stderr, "totallivetime_finalfit: Couldn't open runlist\n");
    return 1;
  }
  double livetime = 0;
  double rrmlivetimes[3] = {0};
  double shortesttime = 1e30;
  int run;
  int num_runs = 0;
  while(infile >> run){
    num_runs++;
    const double thislivetime = getrunlivetime(run);
    livetime += thislivetime;
    rrmlivetimes[reactorpowerbin(run)] += thislivetime;
    if(thislivetime < shortesttime) shortesttime = thislivetime;
  }
  printf("const int num_runs = %d;\n", num_runs); 
  printf("const double shortest_run_s = %.9f;\n", shortesttime); 
  printf("const double livetime_s = %.9f;\n", livetime); 
  printf("const double livetime = %.9f;\n", livetime/86400.); 
  printf("const double rrmlivetimes[3] = { %.9f, %.9f, %.9f };\n",
    rrmlivetimes[0]/86400.,
    rrmlivetimes[1]/86400.,
    rrmlivetimes[2]/86400.); 
  printf("TECHNOTE run time: %.3f days;\n", livetime/86400); 
  return 0;
}
