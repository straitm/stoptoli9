#include "TH1.h"
#include "TF1.h"
#include "TRandom3.h"

static double max(const double a, const double b)
{
  return a>b?a:b;
}

void b12spectrum()
{
  const int nb12b = 3;
  double b12branchbe[nb12b] = {5.7150, 8.9304, 13.3693};
  const double totalb12e = b12branchbe[2];
  double b12branchp[nb12b] = {1.5, 1.23, 97.2};

  {
    double sum = 0;
    for(int i = 0; i < nb12b; i++) sum += b12branchp[i];
    for(int i = 0; i < nb12b; i++) b12branchp[i] /= sum;
    for(int i = 1; i < nb12b; i++) b12branchp[i] += b12branchp[i-1];
  }
    
  const int nb13b = 7;
  double b13branchbe[nb13b] = {3.540, 4.577, 5.890, 9.834, 9.7527, 10.3478, 13.4372};
  const double totalb13e = b13branchbe[6];
  double b13branchp[nb13b] = {0.022, 0.16, 0.094, 0.4, 7.6, 0.3, 92.1};

  {
    double sum = 0;
    for(int i = 0; i < nb13b; i++) sum += b13branchp[i];
    for(int i = 0; i < nb13b; i++) b13branchp[i] /= sum;
    for(int i = 1; i < nb13b; i++) b13branchp[i] += b13branchp[i-1];
  }

  TF1 beta("beta", "(x+0.510998928)*sqrt((x+0.510998928)**2 "
    "- 0.510998928**2) * ([0] - x)**2", 0, 1);


  TH1D * b12spec = new TH1D("b12spec", "", 64, 0, 16);
  TH1D * b13spec = new TH1D("b13spec", "", 64, 0, 16);

  TRandom3 r;
  const int N = 100000;
  for(int c = 0; c < N; c++){
    const double brand = r.Rndm();
    int b = 0;
    for(int i = 0; i < nb12b; i++){
      if(brand < b12branchp[i]){
        b = i;
        break;
      }
    }

    beta.SetParameter(0, b12branchbe[b]);
    beta.SetRange(0, b12branchbe[b]);
   
    const double gammae = totalb12e - b12branchbe[b];
    const double truee = beta.GetRandom() + gammae;

    const double reco_err = pow(0.017,2) + pow(0.077,2) * truee
                          + pow(0.045*truee,2);

    b12spec->Fill(max(0, truee + reco_err));
  }

  for(int c = 0; c < N; c++){
    const double brand = r->Rndm();
    int b = 0;
    for(int i = 0; i < nb13b; i++){
      if(brand < b13branchp[i]){
        b = i;
        break;
      }
    }

    beta.SetParameter(0, b13branchbe[b]);
    beta.SetRange(0, b13branchbe[b]);
   
    const double gammae = totalb13e - b13branchbe[b];
    const double truee = beta.GetRandom() + gammae;

    const double reco_err = pow(0.017,2) + pow(0.077,2) * truee
                          + pow(0.045*truee,2);

    b13spec->Fill(max(0, truee + reco_err));
  }

  b12spec->Draw("e");
  b13spec->SetLineColor(kRed);
  b13spec->SetMarkerColor(kRed);
  b13spec->Draw("esame");
}
