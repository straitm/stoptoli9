#include "TH1.h"
#include "TFile.h"
#include "TF1.h"
#include "TRandom3.h"

static double max(const double a, const double b)
{
  return a>b?a:b;
}

double sv_alpha(const double T, const int Z)
{
  if(T < 1.2) return -0.811 + 4.46e-2*Z + 1.08e-4*Z*Z;
  return -8.46e-2 + 2.48e-2*Z +2.37e-4*Z*Z;
}

double sv_beta(const double T, const int Z)
{
  if(T < 1.2) return 0.673 - 1.82e-2*Z + 6.38e-5*Z*Z;
  return  1.15e-2 + 3.58e-4*Z - 6.17e-5*Z*Z;
}

TF1 * beta = new TF1("beta",
/*   E                            p                       (Q-T)**2 */
"(x+0.51099893)*sqrt((x+0.51099893)**2 - 0.51099893**2) * ([0]-x)**2 *"

/* F(Z,E) as per Schenter and Vogel: */
/*    E         /             p      */
"(x+0.51099893)/sqrt((x+0.51099893)**2 - 0.51099893**2) *"
"exp("
/* S+V's alpha function. Note since Z is always 6, we could evaluate
this in place */
"((x < 1.2*0.51099893)*(-0.811 + 4.46e-2*6 + 1.08e-4*6*6) + (x >= 1.2*0.51099893)*(-8.46e-2 + 2.48e-2*6 +2.37e-4*6*6)) +"
/* S+V's beta function */
"((x < 1.2*0.51099893)*(0.673 - 1.82e-2*6 + 6.38e-5*6*6) + (x >= 1.2*0.51099893)*(1.15e-2 + 3.58e-4*6 - 6.17e-5*6*6))*"
/*                               sqrt(E/m_e - 1)    */
"sqrt((x+0.51099893)/0.51099893 - 1))"
,

0, 1);

void fillit(TH1D * hist, TRandom3 * r, const int nbranch,
            double * branchp, double * branche, const int N)
{
  TF1 * bf[nbranch];
  for(int b = 0; b < nbranch; b++){
    bf[b] = (TF1 *)beta->Clone(Form("bf%d", b));
    bf[b]->SetParameter(0, branche[b]);
    bf[b]->SetRange(0, branche[b]);
  }

  for(int c = 0; c < N; c++){
    if(c%1000 == 999) printf("."), fflush(stdout);
    const double brand = r->Rndm();
    int b = 0;
    for(int i = 0; i < nbranch; i++){
      if(brand < branchp[i]){
        b = i;
        break;
      }
    }

    const double gammae = branche[nbranch-1] - branche[b];
    const double truee = 1.06 /* XXX */ *(bf[b]->GetRandom() + gammae);

    const double reco_err = sqrt(
                               + pow(0.077,2)*truee
                               + pow(0.018,2)*truee*truee
                               + pow(0.017,2)
                               );

    hist->Fill(max(0, truee + r->Gaus(0, reco_err)));
  }
  for(int b = 0; b < nbranch; b++) delete bf[b];
  puts("");
}

void b12spectrum(const char * outfile = NULL)
{
  const int nb12b = 3;
  double b12branchbe[nb12b] = {5.7150, 8.9304, 13.3693};
  double b12branchp[nb12b] = {1.5, 1.23, 97.2};

  {
    double sum = 0;
    for(int i = 0; i < nb12b; i++) sum += b12branchp[i];
    for(int i = 0; i < nb12b; i++) b12branchp[i] /= sum;
    for(int i = 1; i < nb12b; i++) b12branchp[i] += b12branchp[i-1];
  }
    
  const int nb13b = 7;
  double b13branchbe[nb13b] = {3.540, 4.577, 5.890, 9.834, 9.7527, 10.3478, 13.4372};
  double b13branchp[nb13b] = {0.022, 0.16, 0.094, 0.4, 7.6, 0.3, 92.1};

  {
    double sum = 0;
    for(int i = 0; i < nb13b; i++) sum += b13branchp[i];
    for(int i = 0; i < nb13b; i++) b13branchp[i] /= sum;
    for(int i = 1; i < nb13b; i++) b13branchp[i] += b13branchp[i-1];
  }

  const int nli8b = 1;
  double li8branchbe[nli8b] = { 16.00413 };
  double li8branchp[nli8b] = { 1 };


  TH1D * b12spec = new TH1D("b12spec", "", 200, 0, 20);
  TH1D * b13spec = new TH1D("b13spec", "", 200, 0, 20);
  TH1D * li8spec = new TH1D("li8spec", "", 200, 0, 20);

  TRandom3 r(0); // 0 -> UUID seed
  fillit(b12spec, &r, nb12b, b12branchp, b12branchbe, 10000000);
  fillit(b13spec, &r, nb13b, b13branchp, b13branchbe, 10000000);
  fillit(li8spec, &r, nli8b, li8branchp, li8branchbe, 10000000);

  b12spec->Draw("e");
  b13spec->SetLineColor(kRed);
  b13spec->SetMarkerColor(kRed);
  b13spec->Draw("esame");
  li8spec->SetLineColor(kBlue);
  li8spec->SetMarkerColor(kBlue);
  li8spec->Draw("esame");

  if(outfile){
    TFile out(outfile, "create");
    b12spec->Write();
  }
}
