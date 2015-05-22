#include <iostream>
#include <fstream>
#include <algorithm>
#include "consts.h"

#include "TCanvas.h"
#include "TColor.h"
#include "TF1.h"
#include "TFile.h"
#include "TGaxis.h"
#include "TGraph.h"
#include "TH1.h"
#include "TLatex.h"
#include "TMarker.h"
#include "TMath.h"
#include "TMinuit.h"
#include "TROOT.h"
#include "TTree.h"
#include "TRandom3.h"

TCanvas * c1 = new TCanvas("c1", "", 0, 0, 800, 600);

static TRandom3 superrand;

struct parres{
  double val, ep, en;
  bool fix;

  // !!! parameters are assumed to never have limits set !!!
  parres(double _val, double _ep, double _en, bool _fix)
  {
    if(!_fix && (_ep == -54321 || _ep==0))
      fprintf(stderr, "Can't handle lack of positive errors!\n");

    if(!_fix && (_en == -54321 || _en==0)){
      fprintf(stderr, "Faking it given lack of negative errors\n");
      _en = _val;
    }

    val = _val; ep = fabs(_ep), en = fabs(_en),
    fix = _fix;
  }

  double rand() const
  {
    if(fix) return val;
    if(superrand.Rndm() < 0.5) return val+superrand.Rndm()*ep*1.01;
    else                       return val-superrand.Rndm()*en*1.01;
  }

  double min() const
  {
    if(fix) return val;
    return val-en;
  }

  double max() const
  {
    if(fix) return val;
    return val+ep;
  }
};

static double getpar(int i)
{
  double answer, dum;
  gMinuit->GetParameter(i, answer, dum);
  return answer;
}

void fixat(int i, float v)
{
  gMinuit->Command(Form("REL %d", i));
  gMinuit->Command(Form("SET PAR %d %f", i, v));
  gMinuit->Command(Form("FIX %d", i));
}

static void envelope(const double mult)
{
  std::vector<parres> bestpar;
  int par = 5;

  for(int i = 0; i < par; i++)
    bestpar.push_back(parres(getpar(i), gMinuit->fErp[i - (i >= 3)],
      gMinuit->fErn[i - (i >= 3)], i == 2));

  double bestlike;

  {
    double pars[par];
    for(int i = 0; i < par; i++) pars[i] = bestpar[i].val;
    gMinuit->Eval(par, NULL, bestlike, pars, 0);
  }

  vector<TF1 *> bounds;

  int tomult[4] = { 0, 1, 3, 4};

  for(int t = 0; t < par*2; t++){
    const bool lo = t%2;
    const int par1 = t/2;
    if(bestpar[par1].fix) continue;
    if(lo) printf("%d min/max\n", par1);
    for(int X = 0; X < par; X++){
      gMinuit->Command(Form("REL %d\n", X+1));
      gMinuit->Command(Form("SET PAR %d %f",X+1,bestpar[X].val));
    }
    fixat(3, bestpar[2].val);
    fixat(par1+1, lo?bestpar[par1].min():bestpar[par1].max());
    
    for(int y = 0; y < 2; y++) gMinuit->Command("MIGRAD"); 
    const double dchi2 = gMinuit->fAmin-bestlike; 
    if(dchi2 > 0.5 + 0.01 || dchi2 < -1e-5){ 
      printf("MIGRAD gave a delta of %g :-(\n", dchi2); 
    }else{ 
      printf("OK!\n"); 
      double ppar[par];
      for(int V = 0; V < par; V++) ppar[V] = getpar(V);  
      TF1 * f = new TF1(Form("g%d", t), "[0]*exp(-x*log(2)/0.0202) + "
                             "[1]*exp(-x*log(2)/[2]) + "
                             "[3]*exp(-x*log(2)/7.13) + "
                             "[4]", 0, 100);
      for(int i = 0; i < f->GetNpar(); i++)
        f->SetParameter(i, ppar[i]);
      for(int i = 0; i < 4; i++)
        f->SetParameter(tomult[i], f->GetParameter(tomult[i])*mult);
      bounds.push_back(f);
      bounds[bounds.size()-1]->SetLineColor(kRed);
      bounds[bounds.size()-1]->SetLineWidth(1);
      bounds[bounds.size()-1]->Draw("same");
      c1->Modified(); c1->Update();
    }
  }


  for(int t = 0; t < par*(par-1)/2; t++){
    int npar1 = 1, npar2 = 2; 
    for(int moose = 0; moose < t; moose++){
      npar2++;
      if(npar2 > par){
        npar1++;
        npar2 = npar1+1;
      }    
    }    
    if(bestpar[npar1-1].fix || bestpar[npar2-1].fix) continue;

    for(int X = 0; X < par; X++){
      gMinuit->Command(Form("REL %d\n", X+1));
      gMinuit->Command(Form("SET PAR %d %f",X+1,bestpar[X].val));
    }
    fixat(3, bestpar[2].val);

    printf("%d %d contour\n", npar1, npar2);
    gMinuit->Command("MIGRAD");
    gMinuit->Command("Set print 0");
    gMinuit->fGraphicsMode = false;
    gMinuit->Command(Form("MNCONT %d %d 20", npar1, npar2));
    gMinuit->fGraphicsMode = true;
    gMinuit->Command(Form("MNCONT %d %d 20", npar1, npar2));
    TGraph * gr = gMinuit->GetPlot()?(TGraph*)
                  ((TGraph*)gMinuit->GetPlot())->Clone():NULL;
    if(!gr){
      printf("Could not get contour!\n"); continue;
    }    
    gMinuit->Command("Set print -1");
    for(int squirrel = 0; squirrel < gr->GetN(); squirrel++){
      printf("  contour point %d: %f %f\n", squirrel,
             gMinuit->fXpt[squirrel], gMinuit->fYpt[squirrel]);
      fixat(npar1, gr->GetX()[squirrel]);
      fixat(npar2, gr->GetY()[squirrel]);

      for(int y = 0; y < 2; y++) gMinuit->Command("MIGRAD"); 
      const double dchi2 = gMinuit->fAmin - bestlike; 
      if(dchi2 > 0.5 + 0.01 || dchi2 < -1e-5){ 
        printf("MIGRAD gave a delta of %g :-(\n", dchi2); 
      }else{ 
        printf("OK!\n"); 
        double ppar[par];
        for(int V = 0; V < par; V++) ppar[V] = getpar(V);  
        TF1 * f = new TF1(Form("g%d", t), "[0]*exp(-x*log(2)/0.0202) + "
                               "[1]*exp(-x*log(2)/[2]) + "
                               "[3]*exp(-x*log(2)/7.13) + "
                               "[4]", 0, 100);
        for(int i = 0; i < f->GetNpar(); i++)
          f->SetParameter(i, ppar[i]);
        for(int i = 0; i < 4; i++)
          f->SetParameter(tomult[i], f->GetParameter(tomult[i])*mult);

        bounds.push_back(f);
        bounds[bounds.size()-1]->SetLineColor(kBlack);
        bounds[bounds.size()-1]->SetLineWidth(1);
        bounds[bounds.size()-1]->Draw("same");
        c1->Modified(); c1->Update();
      }    
    }    
    delete gr;
  }

  const int nrand = 1000;
  for(int curve = 0; curve < nrand; curve++){
    TF1 * f = new TF1(Form("f%d", curve), "[0]*exp(-x*log(2)/0.0202) + "
                           "[1]*exp(-x*log(2)/[2]) + "
                           "[3]*exp(-x*log(2)/7.13) + "
                           "[4]", 0, 100);
    double like, pars[par];
    do{
      for(int i = 0; i < f->GetNpar(); i++) pars[i] = bestpar[i].rand();
      gMinuit->Eval(par, NULL, like, pars, 0);
    }while(fabs(like - bestlike - 0.5) > 0.01);

    if(curve%10 == 0) printf("%d/%d\n", curve, nrand);

    for(int i = 0; i < f->GetNpar(); i++) f->SetParameter(i, pars[i]);

    for(int i = 0; i < 4; i++)
      f->SetParameter(tomult[i], f->GetParameter(tomult[i])*mult);
    bounds.push_back(f);
    bounds[bounds.size()-1]->SetLineColor(kGreen+2);
    bounds[bounds.size()-1]->SetLineWidth(1);
    bounds[bounds.size()-1]->Draw("same");
    c1->Modified(); c1->Update();
  }

  const double low = 0.1, high = 20;

  TGraph ghigh, glow;
  const int xpoints = 100;
  const double llow = low?low:0.001;
  for(int ix = -1; ix <= xpoints; ix++){
    const double x = ix == -1?0:
       exp(log(llow) + double(ix)*(log(high)-log(llow))/xpoints);
    double lowest = 1e50, highest = -1e50;
    int besttl = -1, bestth = -1;
    for(unsigned int t = 0; t < bounds.size(); t++){
      const double y = bounds[t]->Eval(x);
      if(y < lowest) besttl=t, lowest = y;
      if(y > highest)bestth=t, highest= y;
    }

    ghigh.SetPoint(ghigh.GetN(), x, highest);
    glow .SetPoint( glow.GetN(), x,  lowest);
  }

  TGraph * gall = new TGraph();
  for(int p = glow.GetN()-1; p >= 0; p--)
    gall->SetPoint(gall->GetN(),  glow.GetX()[p],  glow.GetY()[p]);
  for(int p = 0; p < ghigh.GetN(); p++)
    gall->SetPoint(gall->GetN(), ghigh.GetX()[p], ghigh.GetY()[p]);

  gall->SetFillColor(TColor::GetColor("#ccccff"));
  gall->SetFillStyle(1001);
  gall->Draw("f");
  gall->SavePrimitive(cout);
}

TTree ibd;

bool isibd(const int run, const int prompttrig)
{
  return ibd.GetEntries(Form("run==%d && trig==%d", run, prompttrig));
}

void li8finalfit(const int nn)
{
  TFile * fiel = new TFile(rootfile3up, "read");
  TTree * t = (TTree *) fiel->Get("t");

  //ibd.ReadFile("/cp/s4/strait/li9ntuples/Hprompts20141119", "run:trig");
  //ibd.ReadFile("/cp/s4/strait/li9ntuples/Gdprompts20140925");

  const char * const cut =
   nn == 1? // terrible
  "!earlymich && miche<12 && dist<400 && latennear==1 && e>5 && e<14 && timeleft>1e5 && b12like < 0.02": // && !isibd(run, trig)":
   nn == -1?
  "!earlymich && miche<12 && dist<400 &&                 e>5 && e<14 && timeleft>1e5 && b12like < 0.02": // && !isibd(run, trig)":
  "!earlymich && miche<12 && dist<400 && latennear==0 && e>5 && e<14 && timeleft>1e5 && b12like < 0.02"; // && !isibd(run, trig)";

  t->Draw("dt/1000 >> hfit(10000, 0.001, 100)", cut);
  TH1D * hfit = (TH1D *)gROOT->FindObject("hfit");

  TF1 * ee = new TF1("ee", "[0]*exp(-x*log(2)/0.0202) + "
               "[1]*exp(-x*log(2)/[2]) + "
               "[3]*exp(-x*log(2)/7.13) + "
               "[4]", 0, 100);

  ee->SetParameters(1, 1, 0.8399, 1, 1);
  ee->SetParLimits(3, 0, 1);
  ee->FixParameter(2, 0.8399);
  hfit->Fit("ee", "l");
  ee->ReleaseParameter(2);
  hfit->Fit("ee", "l"); // "le"

  printf("%sli8 lifetime: %f +%f %f%s\n", RED,
         ee->GetParameter(2), gMinuit->fErp[2], gMinuit->fErn[2], CLR);

  ee->FixParameter(2, 0.8399);

  hfit->Fit("ee", "le");
  gMinuit->Command("show min");

  t->Draw("dt/1000 >> hdisp(40, 0.1, 20.1)", cut, "e");
  TH1D * hdisp = (TH1D *)gROOT->FindObject("hdisp");

  TF1 * eedisp = (TF1 *)ee->Clone("eedisp");
  eedisp->SetNpx(400);
  eedisp->SetLineColor(kBlack);
  eedisp->SetLineStyle(kDashed);

  int tomult[4] = { 0, 1, 3, 4};
  const double mult = hdisp->GetBinWidth(1)/hfit->GetBinWidth(1);
  for(int i = 0; i < 4; i++)
    eedisp->SetParameter(tomult[i], eedisp->GetParameter(tomult[i])*mult);

  envelope(mult);
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

  printf("%sN found, before efficiency: %f +%f %f%s\n",
         RED, Nfound, Nerrup, Nerrlo, CLR);

  const double eff = 1
    * 0.981 // subsequent muons
    * 0.977 // previous muons
    * wholedet_dist400eff // delta r
    * 0.9709 // 100s from end of run
    * 0.71 // energy
    * (nn == 1?neff_dr_800_avg*neff_dt_avg:1) // neutron, terrible code
    * 0.906 // b12likelihood
  ;

  const double captures = (nn == 0?n_c12cap:nn==-1?n_c12cap+n_c13cap:n_c13cap)*livetime;

  const double toprob = 1./captures/eff;

  printf("neutron Efficiency: %.2f%%\n",  (nn == 1?neff_dr_800_avg*neff_dt_avg:1)*100);
  printf("Total Efficiency: %.2f%%\n", eff*100);
  printf("%sProb per C-%s: %g +%g %g%s\n", 
      RED, nn==1?"13":nn==-1?"nat":"12", toprob*Nfound, toprob*Nerrup, toprob*Nerrlo, CLR);
/*
  TCanvas * c3 = new TCanvas;

  const char * const cutmich =
    "!earlymich && dist < 400 && e > 4 && e < 14 && latennear == 2"
    "&& timeleft > 100e3  && dt/1000 > 0.25 && dt/1000 < 4";

  t->Draw("miche >> emichhist(319, 0.25, 80)", cutmich, "ehist");
*/

}
