#include "consts.h"

void b12gammafinalfit()
{
  TFile * f = new TFile(rootfile0up, "read");
  TTree * t = (TTree *)f->Get("t");

  const double n_c12cap_shortrack = 177.427947;

  const double eff = 1
    * 0.962 // subsequent muons with 1ms veto
    * 0.977 // previous muons
    * 0.9999418 // 200ms from end of run
    * 0.82 // B-12 beta decay energy cut, eval'd for short tracks
    * 0.9177 // dist cut for B-12 beta decay
    * (exp(-3.008*log(2)/2.028)-exp(-5.008*log(2)/2.028)) // B-12 capture time
    * (exp(-2.000*log(2)/20.20)-exp(-200.0*log(2)/20.20)) // B-12 beta decay time
    * 0.95 // michd ???????????????????????????????????
  ;

  t->Draw("miche >> ehist(150, 0.7, 15.7)",
         "!earlymich && "
         "latennear==0 && "
         "ndecay == 0 && "
         "e> 4 && e < 15 && "
         "dt >2 && dt < 200 && "
         "timeleft>200 && "
         "dist < 400 && "
         "gclen < 1500 && "
         "micht >= 3008 && micht < 5008 &&"
         "michd < 1000"
      , "e");

  TH1 * ehist = gROOT->FindObject("ehist");
  ehist->Draw("e");

  TF1 *gg = new TF1("gg", "[0] +"
      " [3]/[2]* exp(-(((x- [4]*[1])/[2])**2)/2) +"
      " [5]/[2]* exp(-(((x- [6]*[1])/[2])**2)/2) +"
      " [7]/[2]* exp(-(((x- [8]*[1])/[2])**2)/2) +"
      " [9]/[2]* exp(-(((x-[10]*[1])/[2])**2)/2) +"
      "[11]/[2]* exp(-(((x-[12]*[1])/[2])**2)/2) +"
      "[13]/[2]* exp(-(((x-[14]*[1])/[2])**2)/2)  "
      , 0,15);

  gg->SetParName(0, "accidentals");
  gg->SetParName(1, "energyscale");
  gg->SetParName(2, "resolution");
  gg->SetParName(3, "n1");
  gg->SetParName(4, "e1");
  gg->SetParName(5, "n2");
  gg->SetParName(6, "e2");
  gg->SetParName(7, "n3");
  gg->SetParName(8, "e3");
  gg->SetParName(9, "n4");
  gg->SetParName(10, "e4");
  gg->SetParName(11, "n5");
  gg->SetParName(12, "e5");
  gg->SetParName(13, "n6");
  gg->SetParName(14, "e6");
  gg->SetNpx(400);

  gg->FixParameter(1, 1);
  gg->FixParameter(2, 0.4);

  gg->FixParameter(4, 0.95314);
  gg->FixParameter(6, 1.67365);
  gg->FixParameter(8, 2.6208);
  gg->FixParameter(10,2.723);
  gg->FixParameter(12,3.759);
  gg->FixParameter(14,9.040);

  gg->SetParLimits(0, 0, 100);
  for(int i = 3; i <= 13; i+=2) gg->SetParLimits(i, 0, 100);

  ehist->Fit("gg", "l");
  gg->ReleaseParameter(2);
  ehist->Fit("gg", "l");
  gg->ReleaseParameter(1);
  ehist->Fit("gg", "liem");

  TF1 * mypeak = new TF1("mypeak",
    "[2]/[1]* exp(-(((x- [3]*[0])/[1])**2)/2)", 0, 15);
  mypeak->SetNpx(400);

  for(int i = 0; i < 2; i++) mypeak->SetParameter(i, gg->GetParameter(i+1));
  for(int i = 0; i < 2; i++) mypeak->SetParameter(i+2, gg->GetParameter(i+7));

  mypeak->Draw("same");

  printf("%sEvents, eff corrected: %.1f%s\n",
         RED, mypeak->Integral(0, 15)/eff, CLR);
}
