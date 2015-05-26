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

struct ve{
  double val, eup, elo;
  ve(const double val_, const double eup_, const double elo_)
  {
    val = val_, eup = eup_, elo = elo_;
  }
};

ve li8finalfit(const int nn, const bool excludeibd = false,
               const bool quiet = false)
{
  static int stacklevel = 0;
  stacklevel++;
  TFile * fiel = new TFile(rootfile3up, "read");
  TTree * t = (TTree *) fiel->Get("t");

  ibd.ReadFile("/cp/s4/strait/li9ntuples/Hprompts20141119", "run:trig");
  ibd.ReadFile("/cp/s4/strait/li9ntuples/Gdprompts20140925");

  char cut[1000];
  snprintf(cut, 999, "!earlymich && miche<12 && dist<400 %s "
           "&& e>5 && e<14 && timeleft>1e5 && b12like < 0.02 %s",
    nn == 1? "&& latennear==1": nn == -1?  "": "&& latennear==0",
    excludeibd?"&& !isibd(run, trig)":"");

  TH1D * hfit = new TH1D(Form("hfit%d", nn), "", 10000, 0.001, 100);

  t->Draw(Form("dt/1000 >> hfit%d", nn), cut);

  TF1 * ee = new TF1("ee", "[0]*exp(-x*log(2)/0.0202) + "
               "[1]*exp(-x*log(2)/[2]) + "
               "[3]*exp(-x*log(2)/7.13) + "
               "[4]", 0, 100);

  ee->SetParameters(1, 1, 0.8399, 1, 1);
  ee->SetParLimits(3, 0, 1);
  ee->FixParameter(2, 0.8399);
  hfit->Fit("ee", "lq");
  ee->ReleaseParameter(2);
  hfit->Fit("ee", "leq");
  if(!quiet) gMinuit->Command("show min");

  if(!quiet) printf("%sli8 lifetime: %f +%f %f%s\n", RED,
         ee->GetParameter(2), gMinuit->fErp[2], gMinuit->fErn[2], CLR);

  ee->FixParameter(2, 0.8399);

  hfit->Fit("ee", "leq");
  if(!quiet) gMinuit->Command("show min");
  
  if(!quiet){
    TH1D * hdisp = new TH1D(Form("hdisp%d", nn), "", 40, 0.1, 20.1);
    t->Draw(Form("dt/1000 >> hdisp%d", nn), cut, "e");

    TF1 * eedisp = (TF1 *)ee->Clone("eedisp");
    eedisp->SetNpx(400);
    eedisp->SetLineColor(kBlack);
    eedisp->SetLineStyle(kDashed);

    int tomult[4] = { 0, 1, 3, 4};
    const double mult = hdisp->GetBinWidth(1)/hfit->GetBinWidth(1);
    for(int i = 0; i < 4; i++)
      eedisp->SetParameter(tomult[i], eedisp->GetParameter(tomult[i])*mult);

    //envelope(mult);
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
  }

  const double eff = 1
    * 0.981 // subsequent muons
    * 0.977 // previous muons
    * wholedet_dist400eff // delta r
    * 0.9709 // 100s from end of run
    * 0.71 // energy
    * (nn == 1?neff_dr_800_avg*neff_dt_avg:1) // neutron, terrible code
    * 0.906 // b12likelihood
  ;

  const ve Nfound_uncorr(
    ee->GetParameter(1)/0.8399/hfit->GetBinWidth(1),
    ee->GetParameter(1)/0.8399/hfit->GetBinWidth(1) *
      gMinuit->fErp[1]/ee->GetParameter(1),
    ee->GetParameter(1)/0.8399/hfit->GetBinWidth(1) *
      gMinuit->fErn[1]/ee->GetParameter(1));

  if(!quiet) printf("%sN found, before efficiency: %.3f +%.3f %.3f%s\n",
         RED, Nfound_uncorr.val, Nfound_uncorr.eup, Nfound_uncorr.elo, CLR);
  const double captures = (nn == 0?n_c12cap:nn==-1?n_c12cap+n_c13cap:n_c13cap)*livetime;

  ve Nfound_corr1(0, 0, 0);
  if(nn == 0){
    for(int i = 0; i < stacklevel; i++) printf(" ");
    printf("Finding correction for 1-neutron events with inefficiency...\n");
    ve correction = li8finalfit(1, excludeibd, true);
    correction.val *= (1 - neff_dr_800_avg*neff_dt_avg);
    correction.eup *= (1 - neff_dr_800_avg*neff_dt_avg);
    correction.elo *= (1 - neff_dr_800_avg*neff_dt_avg);
    for(int i = 0; i < stacklevel; i++) printf(" ");
    printf("Subtracting %.3g events for 1-neutron events\n", correction.val);
    Nfound_corr1.val = Nfound_uncorr.val - correction.val;
    Nfound_corr1.elo = sqrt(pow(Nfound_uncorr.elo, 2) - pow(correction.elo, 2));
    Nfound_corr1.eup = sqrt(pow(Nfound_uncorr.eup, 2) - pow(correction.eup, 2));
  }
  else{
    Nfound_corr1 = Nfound_uncorr;
  }

  ve Nfound_corr2(0, 0, 0);
  if(!excludeibd){
    for(int i = 0; i < stacklevel; i++) printf(" ");
    printf("Finding correction for IBD events...\n");
    ve correction = li8finalfit(nn, true, true);
    
    correction.val = Nfound_corr1.val - correction.val;

    // zero-neutron assumed to be all Li-9, one-neutron N-17
    const double bn_prob = nn == 0? 0.492: 0.951;
    correction.val /= bn_prob;

    for(int i = 0; i < stacklevel; i++) printf(" ");
    printf("Subtracting %.3g events for IBD events\n", correction.val);
    Nfound_corr2.val = Nfound_uncorr.val - correction.val;
    Nfound_corr2.elo = Nfound_corr1.elo; // neglect error on correction
    Nfound_corr2.eup = Nfound_corr1.eup;
  }
  else{
    Nfound_corr2 = Nfound_corr1;
  }

  const double toprob = 1./captures/eff;

  if(!quiet){
    printf("neutron Efficiency: %.2f%%\n",  (nn == 1?neff_dr_800_avg*neff_dt_avg:1)*100);
    printf("Total Efficiency: %.2f%%\n", eff*100);
    printf("%sProb per C-%s: %.3g +%.3g %.3g%s\n", 
        RED, nn==1?"13":nn==-1?"nat":"12", toprob*Nfound_corr2.val,
        toprob*Nfound_corr2.eup, toprob*Nfound_corr2.elo, CLR);
    printf("%sSystematic error: %.3g%s\n", RED,
      (nn==1?n_c13cap_err/n_c13cap:n_c12cap_err/n_c12cap)*toprob*Nfound_corr2.val, CLR);
  }

  stacklevel--;
  return Nfound_corr2;
}
