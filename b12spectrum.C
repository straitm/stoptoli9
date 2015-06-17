#include "TH1.h"
#include "TFile.h"
#include "TF1.h"
#include "TRandom3.h"

TRandom3 r(0); // 0 -> UUID seed

TH1D * b12spec = new TH1D("b12spec", "", 200, 0, 20);
TH1D * b13spec = new TH1D("b13spec", "", 200, 0, 20);
TH1D * li8spec = new TH1D("li8spec", "", 200, 0, 20);

const int nb12b = 4;
/*                                            forbidden?              allowed       allowed */
double b12branchbe[nb12b] = {5.7150,                        5.7150,         8.9304, 13.3693};
double b12branchge[nb12b] = {0.12*(13.3693-5.7150), 13.3693-5.7150, 13.3693-8.9304, 0};
double tmp_b12branchp[nb12b] = {1.58*0.9, 1.58*0.1, 1.23, -1};
double b12branchp[nb12b] = {tmp_b12branchp[0], tmp_b12branchp[0], tmp_b12branchp[0], 100 - tmp_b12branchp[0] - tmp_b12branchp[1]- tmp_b12branchp[2] };
bool   b12brancha[nb12b] = {0, 0, 1, 1};
double b12branch_aw[nb12b] = {0, 0, 0.0, 0.0093/0.511};

const int nb13b = 8;
double b13branchbe[nb13b] = {8.4909, 3.540, 4.577, 5.890, 9.834, 9.7527, 10.3478, 13.4372};
double b13branchge[nb13b] = {4.438,  13.4372-3.540, 13.4372-4.577, 13.4372-5.890, 13.4372-9.834, 13.4372-9.7527, 13.4372-10.3478, 0};
double b13branchp[nb13b] = {0.0029, 0.022, 0.16, 0.094, 0.4, 7.6, 0.3, 92.1};
bool   b13brancha[nb13b] = {1, 1, 1, 1, 1, 1, 1, 1}; // no idea
double b13branch_aw[nb13b] = {0, 0, 0, 0, 0, 0, 0, 0}; // no idea

const int nli8b = 1;
double li8branchbe[nli8b] = { 16.00413 - 3.050 }; // but it is broad!
double li8branchge[nli8b] = { 0 };
double li8branchp[nli8b] = { 1 };
bool   li8brancha[nli8b] = { 1 };
double li8branch_aw[nli8b] = { 0 };

static double max(const double a, const double b)
{
  return a>b?a:b;
}

const char * const allowedspectrum = 

"(1 + [1]*x)*" // shape correction for second-forbidden transitions...?

/*   E                            p                       (Q-T)**2 */
"(x+0.51099893)*sqrt((x+0.51099893)**2 - 0.51099893**2) * ([0]-x)**2 *" // pE
//"(x)*sqrt((x+0.51099893)**2 - 0.51099893**2) * ([0]-x)**2 *" // pT

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
"sqrt((x+0.51099893)/0.51099893 - 1))";

// I think this is a fairly extreme but reasonable case
const char * const forbiddencorrection = "x";


TF1 * abeta = new TF1("abeta", allowedspectrum, 0, 1);
TF1 * fbeta = new TF1("fbeta", Form("(%s) * (%s)", allowedspectrum, forbiddencorrection), 0, 1);

void fillit(TH1D * hist, const int nbranch, double * branchp,
            double * bbranche, double * gbranche, bool * abranch,
   double * branch_aw,
            const int N, const double escale)
{
  TF1 * bf[nbranch];
  for(int b = 0; b < nbranch; b++){
    bf[b] = (TF1 *)(abranch[b]?abeta:fbeta)->Clone(Form("bf%d", b));
    bf[b]->SetParameter(0, bbranche[b]);
    bf[b]->SetParameter(1, branch_aw[b]);
    bf[b]->SetRange(0, bbranche[b]);
  }

  for(int c = 0; c < N; c++){
    if(c%100000 == 99999) printf("."), fflush(stdout);
    const double brand = r.Rndm();
    int b = 0;
    for(int i = 0; i < nbranch; i++){
      if(brand < branchp[i]){
        b = i;
        break;
      }
    }

    const double truee = escale*(bf[b]->GetRandom() + gbranche[b]);

    const double reco_err = sqrt(
                               + pow(0.077+0.002*2,2)*truee
                               + pow(0.018+0.001*10,2)*truee*truee
                               + pow(0.017+0.011*2,2)
                               );

    hist->Fill(max(0, truee + r.Gaus(0, reco_err)));
  }
  for(int b = 0; b < nbranch; b++) delete bf[b];
  puts("");
}

void fillb12(const double escale)
{
  fillit(b12spec, nb12b, b12branchp, b12branchbe, b12branchge, b12brancha, b12branch_aw, 10000000, escale);
}

void norm()
{
  {
    double sum = 0;
    for(int i = 0; i < nb12b; i++) sum += b12branchp[i];
    for(int i = 0; i < nb12b; i++) b12branchp[i] /= sum;
    for(int i = 1; i < nb12b; i++) b12branchp[i] += b12branchp[i-1];
  }
    
  {
    double sum = 0;
    for(int i = 0; i < nb13b; i++) sum += b13branchp[i];
    for(int i = 0; i < nb13b; i++) b13branchp[i] /= sum;
    for(int i = 1; i < nb13b; i++) b13branchp[i] += b13branchp[i-1];
  }
}

void b12spectrum()
{
  norm();
  fillit(b12spec, nli8b, li8branchp, li8branchbe, li8branchge,li8brancha, li8branch_aw, 1e6,1);
}
