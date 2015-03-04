#include "consts.h"

void li8finalfit(const int nn)
{
  const char * const RED     = "\033[31;1m"; // bold red
  const char * const CLR      = "\033[m"    ; // clear

  TFile * fiel = new TFile(rootfile3up, "read");
  TTree * t = (TTree *) fiel->Get("t");

  const char * const cut =
   nn == 1? // terrible
  "!earlymich && miche<12 && dist<400 && latennear==1 && e>5 && e<14 && timeleft>1e5 && b12like < 0.02":
   nn == -1?
  "!earlymich && miche<12 && dist<400 &&                 e>5 && e<14 && timeleft>1e5 && b12like < 0.02":
  "!earlymich && miche<12 && dist<400 && latennear==0 && e>5 && e<14 && timeleft>1e5 && b12like < 0.02";

  TCanvas * c1 = new TCanvas;
  c1->SetLogy();

  t->Draw("dt/1000 >> hfit(10000, 0.001, 100)", cut);

  TF1 * ee->= new TF1("ee", "[0]*exp(-x*log(2)/0.0202) + "
               "[1]*exp(-x*log(2)/[2]) + "
               "[3]*exp(-x*log(2)/7.13) + "
               "[4]", 0, 100);

  ee->SetParameters(1, 1, 0.8399, 1, 1);
  ee->SetParLimits(3, 0, 1);
  ee->FixParameter(2, 0.8399);
  hfit->Fit("ee", "l");
  ee->ReleaseParameter(2);
  hfit->Fit("ee", "le");

  printf("%sli8 lifetime: %f +%f %f%s\n", RED,
         ee->GetParameter(2), gMinuit->fErp[2], gMinuit->fErn[2], CLR);

  ee->FixParameter(2, 0.8399);

  hfit->Fit("ee", "le");

  t->Draw("dt/1000 >> hdisp(40, 0.1, 20.1)", cut, "e");

  TF1 * eedisp = ee->Clone("eedisp");
  eedisp->SetNpx(400);
  eedisp->SetLineColor(kRed);

  int tomult[4] = { 0, 1, 3, 4};
  const double mult = hdisp->GetBinWidth(1)/hfit->GetBinWidth(1);
  for(int i = 0; i < 4; i++)
    eedisp->SetParameter(tomult[i], eedisp->GetParameter(tomult[i])*mult);
  eedisp->Draw("same");

  TF1 * b12 = new TF1("b12", "[0]*exp(-x*log(2)/0.0202)" , 0, 100);
  TF1 * li8 = new TF1("li8", "[0]*exp(-x*log(2)/0.8399)", 0, 100);
  TF1 * n16 = new TF1("n16", "[0]*exp(-x*log(2)/7.13)"  , 0, 100);
  TF1 * acc = new TF1("acc", "[0]", 0, 100);

  b12->SetNpx(400);

  TF1 * parts[4] = { b12, li8, n16, acc };


  b12->SetParameter(0, eedisp->GetParameter(0));
  li8->SetParameter(0, eedisp->GetParameter(1));
  n16->SetParameter(0, eedisp->GetParameter(3));
  acc->SetParameter(0, eedisp->GetParameter(4));

  hdisp->GetYaxis()->SetRangeUser(acc->GetParameter(0)/50,
                                  acc->GetParameter(0)*10);

  for(int i = 0; i < 4; i++){
    parts[i]->SetLineStyle(7);
    parts[i]->SetLineWidth(2);
    parts[i]->Draw("Same");
  } 

  const double Nfound = li8->Integral(0, 20)/hdisp->GetBinWidth(1);
  const double Nerrup = Nfound * gMinuit->fErp[1]/ee->GetParameter(1);
  const double Nerrlo = Nfound * gMinuit->fErn[1]/ee->GetParameter(1);

  printf("%sN found: %f +%f %f%s\n", RED, Nfound, Nerrup, Nerrlo, CLR);

  const double eff = 1
    * 0.981 // subsequent muons
    * 0.977 // previous muons
    * 0.9156 // delta r
    * 0.9709 // 100s from end of run
    * 0.71 // energy
    * (nn == 1?overalllateneff:1) // neutron, terrible code
    * 0.906 // b12likelihood
  ;

  const double captures = (nn == 0?n_c12cap:nn==-1?n_c12cap+n_c13cap:n_c13cap)*livetime;

  const double toprob = 1./captures/eff;

  printf("Efficiency: %.2f%%\n", eff*100);
  printf("%sProb per C-%s: %g +%g %g%s\n", 
      RED, nn==1?"13":n==-1?"nat":"12", toprob*Nfound, toprob*Nerrup, toprob*Nerrlo, CLR);
/*
  TCanvas * c3 = new TCanvas;

  const char * const cutmich =
    "!earlymich && dist < 400 && e > 4 && e < 14 && latennear == 2"
    "&& timeleft > 100e3  && dt/1000 > 0.25 && dt/1000 < 4";

  t->Draw("miche >> emichhist(319, 0.25, 80)", cutmich, "ehist");
*/

}
