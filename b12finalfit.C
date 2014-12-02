#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TH1.h"
#include "TROOT.h"
#include <stdio.h>

const char * const RED     = "\033[31;1m"; // bold red
const char * const CLR      = "\033[m"    ; // clear

static void printfr(const char * const msg, ...)
{
  va_list ap;
  va_start(ap, msg);
  printf(RED);
  vprintf(msg, ap);
  printf(CLR);
}

double eff = 0.82 * 0.977 * 0.962 * 0.9709;
double b12frac = 0.186;
double livedays = 489.509;
double ferror1213 = 0.01;
double ferrorecut = 0.01;
double ferrorb12 = 0.7/18.6;
double f12 = 0.9893;
double ferrorc13 = 0.08/1.07;
double p13op12 = 0.93;
double ferrorp13op12 = 0.06/0.93;
 
void all(TTree * t, TF1 * ee)
{
  t->Draw("dt/1000 >> h(10000, 0.001, 100)",
          "timeleft > 100e3 && miche < 12 && e > 4 && e < 14.5 "
          "&& !earlymich", "e");
  TH1D * h = (TH1D*)gROOT->FindObject("h");
  h->Fit("ee", "ql");
  h->Fit("ee", "l");
  TF1 e("e", "[0]*exp(-x*log(2)/0.0202)", 0, 100);
  e.SetParameter(0, ee->GetParameter(0));

  const double ferrorfit = ee->GetParError(0)/ee->GetParameter(0);
  const double rawintegral = e.Integral(0, 10)/h->GetBinWidth(1);

  double ferror = sqrt(pow(ferrorfit, 2) + pow(ferrorecut, 2));
  printf("b12 raw %f +- %f\n", rawintegral, ferrorfit * rawintegral);
  double ferror2 = sqrt(pow(ferrorfit,  2)+ pow(ferrorecut, 2)+
                        pow(ferror1213, 2)+ pow(ferrorb12, 2));
  printfr("per day c12 cap rate %f +- %f\n",
         rawintegral / eff / b12frac / livedays * f12,
         rawintegral / eff / b12frac / livedays * f12 * ferror2);

  double ferror3 = sqrt(pow(ferrorfit,  2)+ pow(ferrorecut, 2)+
                        pow(ferror1213, 2)+ pow(ferrorb12, 2)+
                        pow(ferrorc13, 2) + pow(ferrorp13op12, 2));

  printfr("per day c13 cap rate %f +- %f\n",
         rawintegral / eff / b12frac / livedays * (1-f12) * p13op12 ,
         rawintegral / eff / b12frac / livedays * (1-f12) * p13op12 * ferror3);
}

void targ(TTree * t, TF1 * ee)
{
  printf("\nIn the target only:\n");
  t->Draw("dt/1000 >> h(10000, 0.001, 100)",
          "dx**2 + dy**2 < 1150**2 && "
          "abs(dz) < 1229+0.03*(1150-sqrt(dx**2+dy**2)) &&"
          "timeleft > 100e3 && miche < 12 && "
          "e > 4 && e < 14.5 && !earlymich", "e");
  TH1D * h = (TH1D*)gROOT->FindObject("h");
  h->Fit("ee", "l");
  TF1 e("e", "[0]*exp(-x*log(2)/0.0202)", 0, 100);
  e.SetParameter(0, ee->GetParameter(0));

  const double ferrorfit = ee->GetParError(0)/ee->GetParameter(0);
  const double rawintegral = e.Integral(0, 10)/h->GetBinWidth(1);

  double ferror = sqrt(pow(ferrorfit, 2) + pow(ferrorecut, 2));
  printf("b12 raw %f +- %f\n", rawintegral, ferrorfit * rawintegral);
  double ferror2 = sqrt(pow(ferrorfit, 2)+ pow(ferrorecut, 2)+
                        pow(ferror1213, 2)+ pow(ferrorb12, 2));
  printfr("per day c12 cap rate %f +- %f\n",
         rawintegral / eff / b12frac / livedays * f12,
         rawintegral / eff / b12frac / livedays * f12 * ferror2);

  double ferror3 = sqrt(pow(ferrorfit,  2)+ pow(ferrorecut, 2)+
                        pow(ferror1213, 2)+ pow(ferrorb12, 2)+
                        pow(ferrorc13, 2) + pow(ferrorp13op12, 2));

  printfr("per day c13 cap rate %f +- %f\n",
         rawintegral / eff / b12frac / livedays * (1-f12) * p13op12 ,
         rawintegral / eff / b12frac / livedays * (1-f12) * p13op12 * ferror3);
}

void innertarg(TTree * t, TF1 * ee, const double limr, const double limz)
{
  printf("\nIn the inner %.1fx%.1f target only:\n", limr, limz);
  t->Draw(Form("dt/1000 >> h%.0f(10000, 0.001, 100)", limr),
          Form("dx**2 + dy**2 < %f**2 && "
          "abs(dz) < %f &&"
          "timeleft > 100e3 && miche < 12 && "
          "e > 4 && e < 14.5 && !earlymich", limr, limz), "e");
  TH1D * h = (TH1D*)gROOT->FindObject(Form("h%.0f", limr));
  h->Fit("ee", "l");
  TF1 e("e", "[0]*exp(-x*log(2)/0.0202)", 0, 100);
  e.SetParameter(0, ee->GetParameter(0));

  const double ferrorfit = ee->GetParError(0)/ee->GetParameter(0);
  const double rawintegral = e.Integral(0, 10)/h->GetBinWidth(1);

  double ferror = sqrt(pow(ferrorfit, 2) + pow(ferrorecut, 2));
  printf("b12 raw %f +- %f\n", rawintegral, ferrorfit * rawintegral);
  double ferror2 = sqrt(pow(ferrorfit, 2)+ pow(ferrorecut, 2)+
                        pow(ferror1213, 2)+ pow(ferrorb12, 2));
  printfr("per day c12 cap rate %f +- %f (stat) %f (tot)\n",
         rawintegral / eff / b12frac / livedays * f12,
         rawintegral / eff / b12frac / livedays * f12 * ferrorfit,
         rawintegral / eff / b12frac / livedays * f12 * ferror2);

  double ferror3 = sqrt(pow(ferrorfit,  2)+ pow(ferrorecut, 2)+
                        pow(ferror1213, 2)+ pow(ferrorb12, 2)+
                        pow(ferrorc13, 2) + pow(ferrorp13op12, 2));

  printfr("per day c13 cap rate %f +- %f\n",
         rawintegral / eff / b12frac / livedays * (1-f12) * p13op12 ,
         rawintegral / eff / b12frac / livedays * (1-f12) * p13op12 * ferror3);
}

void b12finalfit()
{
  TFile *_file0 = TFile::Open(
    "/cp/s4/strait/fullfido-100s-3-25MeV-20141022.root");
  TTree * t = (TTree *)_file0->Get("t");
  TF1 * ee = new TF1("ee", "[0]*exp(-x*log(2)/0.0202) + "
     "[1]*exp(-x*log(2)/0.8399) + "
     "[2]*exp(-x*log(2)/[3]) + "
     "[4]*exp(-x*log(2)/[5]) + "
     "[6]*exp(-x*log(2)/[7]) + "
     "[8]", 0, 100);
  ee->SetParLimits(1, 0, 10);

  ee->SetParLimits(2, 0, 30);
  ee->SetParameter(3, 7.13);
  ee->SetParLimits(3, 1.5, 10);

  ee->SetParLimits(4, 0, 30);
  ee->SetParameter(5, 15);
  ee->SetParLimits(5, 10, 20);

  ee->SetParLimits(6, 0, 30);
  ee->SetParameter(7, 0.2);
  ee->SetParLimits(7, 0.05, 0.4);
  ee->FixParameter(6, 0);
  ee->FixParameter(7, 0);

  all(t, ee);
  ee->FixParameter(4, 0);
  ee->FixParameter(5, 15);
  targ(t, ee);
  innertarg(t, ee, 1045, 1127);
  ee->FixParameter(2, 0);
  ee->FixParameter(3, 7.13);
  innertarg(t, ee, 913, 985);
  innertarg(t, ee, 724, 781);
}
