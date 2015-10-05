#include "consts.h"
#include "TRandom3.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TFile.h"
#include "TTree.h"

// PRD 76, 052003 eq (18).  Returns mu+/mu- ratio.
// This uses plain old costheta, no star.
// In PRD 83, 032011, eq(7) is the same, but it uses costheta-star.
// Probably doesn't matter, since low-angle muons aren't very
// important to us.
static double chargeratio(const double energy, const double costheta, 
   const double fpiplus, const double fkplus)
{
  return
  (          fpiplus/(1+(1.1*energy * costheta)/115.)
   + 0.054 * fkplus /(1+(1.1*energy * costheta)/850.) )
  /
  (          (1-fpiplus)/(1+(1.1*energy * costheta)/115.)
   + 0.054 * (1-fkplus )/(1+(1.1*energy * costheta)/850.) );
}

// Returns the fraction of mu- at a given surface energy and angle.
static double fracmuminus(const double energy, const double costheta,
const double fpiplus, const double fkplus)
{
  return 1/(chargeratio(energy, costheta, fpiplus, fkplus) + 1);
}


void musicchargeratiofinalfit()
{
  TFile * music = new TFile(rootfilemusic);
  TTree * t = (TTree *)music->Get("mu_initial");

  float costheta, energy;
  t->SetBranchAddress("cz_ini", &costheta);
  t->SetBranchAddress("Emu_ini", &energy);

  TRandom3 ran;
  // using my result without going through the derivation again
  const double l3csyst = 0.32/44.12;

  // From MINOS ND charge ratio paper, eq 8
  const double minossyst =
    (1/2.266 - 1/(2.266 + sqrt(pow(0.001,2)+pow(0.0145,2))))*2.266;

  // 0.0: Use MINOS ND data
  // 1.0: Use L3+C data
  const double weights[3] = { 0.0,
    (1/pow(l3csyst,2))/(1/pow(minossyst,2)+1/pow(l3csyst,2)),
    1.0 };

  for(int exper = 0; exper < 3; exper++){
    const double weight = weights[exper];

    TH1D results("results", "", 200, 0.43, 0.45);

    for(int trial = 0; trial < 10000; trial++){
      // Since the MINOS papers don't quote errors, assume that we are
      // in sig-fig territory -- equally likely between the posts. Sigh.
      const double fpiplus = 0.545 + ran.Rndm()*0.01; 
      const double fkplus  = ran.Rndm() < weight?
                           0.665 + ran.Rndm()*0.01:
                           0.695 + ran.Rndm()*0.01;

      double fracsum = 0;

      for(int i = 0; i < t->GetEntries(); i++){
        t->GetEntry(i);
        fracsum += fracmuminus(energy, costheta, fpiplus, fkplus);
      }

      const double frac = fracsum/t->GetEntries();

      results.Fill(frac);
    }

    if     (exper == 0) printf("MINOS ND: ");
    else if(exper == 1) printf("Combined: ");
    else                printf("L3+C:     ");
    printf("%.4f +- %.4f (sig-figs) or +- %.4f (their systematic (shadier))\n",
           results.GetMean(),
           results.GetRMS(),
           results.GetMean() * (exper == 0? minossyst:
                                 exper == 2? l3csyst:
                     1/sqrt(1/pow(l3csyst,2) + 1/pow(minossyst,2)))    );
  }
}
