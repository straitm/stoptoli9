#include "time.h"
#include "TTree.h"
#include "TFile.h"
#include "TError.h"
#include <fstream>
using std::ifstream;

double getrunlivetime(const int run)
{
  //gErrorIgnoreLevel = kFatal;
  TFile rf(Form("/cp/s4/dchooz/cheetah/prod-08-05_p01_v2/reduced.Run%07d_Seq010.root", run), "read");
  //TFile rf11(Form("/cp/s4/dchooz/cheetah/seq11/reduced.Run%07d_Seq011.root", run), "read");
  //gErrorIgnoreLevel = kWarning;
  TTree * data = (TTree*)rf.Get("data");
  //if(!data || data->IsZombie()) data = (TTree*)rf11.Get("data");

  double trgtime;
  data->SetBranchAddress("trgtime", &trgtime);
  data->GetEntry(data->GetEntries()-1);

  // Time from the first event to the last event -- the first event 
  // always has trgtime zero.
  const double obviousseconds = trgtime/1e9;

  // On average, the run was live for one more averge interval between
  // events, half at the beginning and half at the end.
  const double endseconds = obviousseconds/(data->GetEntries()-1);

  const double runtime = obviousseconds+endseconds;

  printf("%.9lf\n", runtime);

  delete data;
  rf.Close();
  //rf11.Close();

  return runtime;
}

double totallivetime_finalfit()
{
  ifstream infile("runlist");
  double runtime = 0;
  double shortesttime = 1e30;
  int run;
  int num_runs = 0;
  while(infile >> run){
    num_runs++;
    const double thisruntime = getrunlivetime(run);
    runtime += thisruntime;
    if(thisruntime < shortesttime) shortesttime = thisruntime;
  }
  printf("const int num_runs = %d\n", num_runs); 
  printf("const double runtime_s = %.9f\n", runtime); 
  printf("const double shortest_run_s = %.9f\n", shortesttime); 
  printf("TECHNOTE run time: %.3f days\n", runtime/86400); 
  return runtime;
}
