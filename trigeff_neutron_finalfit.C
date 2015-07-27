#include "TError.h"
#include "TMarker.h"
#include "TCanvas.h"
#include "TArrow.h"
#include "TTree.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TRandom3.h"
#include "TH1.h"
#include "TF1.h"
#include "TMinuit.h"
#include <vector>
using std::vector;
#include <string>
using std::string;
using std::pair;

TTree * t = NULL;

TRandom3 r(0);

TF1 * makef()
{
  // [0] t_cpt H
  // [1] t_cpt Gd
  // [2] t_th
  // [3] norm H
  // [4] norm Gd
  // [5] bg
  
  const string alphah = "([0]*[2]/([0] - [2]))";

  const string alphagd = "([1]*[2]/([1] - [2]))";

  const string f =
    "abs([3])/([0] - [2]) * exp(-x/[0]) * "
    "(1-exp(-x/" + alphah + ")) + "
    "abs([4])/([1] - [2]) * exp(-x/[1]) * "
    "(1-exp(-x/" + alphagd + ")) + abs([5])";

  TF1 * two = new TF1("two", f.c_str(), 0.01, 100);

  two->FixParameter(0, 190); // capture H
  two->FixParameter(1, 26.15); // capture Gd
  two->FixParameter(2, 0.1902); // epithermal
  two->SetParameter(3, 1000);
  two->SetParameter(4, 1000);
  two->SetParameter(5, 1);

  two->SetParLimits(3, 0, 1e5);
  two->SetParLimits(4, 0, 1e5);

  two->SetParName(0, "hcaptime");
  two->SetParName(1, "gdcaptime");
  two->SetParName(2, "thtime");
  two->SetParName(3, "hnorm");
  two->SetParName(4, "gdnorm");
  two->SetParName(5, "bg");

  two->SetLineColor(kRed);
 
  return two;
}

TTree * readfile()
{
  puts("Reading file");

  TFile * f = new TFile("/cp/s4/strait/fullfido-neutron-20150722.d/all-withpost3rdpub.root", "read");
  TTree * tree = (TTree *)f->Get("t");
  tree->SetName("tree");

  //new TTree("t", "t");
  //tt->ReadFile("/cp/s4/strait/fullfido-neutron-20150722.d/all",
  //             "i/I:nn/I:fq/F:t/F:e/F:ov/I:fqiv/F:x/F:y/F:z/F");
  return tree;
}

vector<double> getbestfit(TF1 * f)
{
  vector<double> ans;
  for(int i = 0; i < f->GetNpar(); i++)
    ans.push_back(f->GetParameter(i));
  return ans;
}

const double mult = 3; // number of sigma wide the box is

vector<double> gethigbounds(TF1 * f)
{
  vector<double> ans;
  for(int i = 0; i < f->GetNpar(); i++){
    ans.push_back(f->GetParameter(i) + f->GetParError(i)*mult);
    double lowlim, higlim;
    f->GetParLimits(i, lowlim, higlim);

    const bool limits = !(lowlim == 0 && higlim == 0);

    if(limits && ans[i] > higlim) ans[i] = higlim;
  }
  return ans;
}

vector<double> getlowbounds(TF1 * f)
{
  vector<double> ans;
  for(int i = 0; i < f->GetNpar(); i++){
    ans.push_back(f->GetParameter(i) - f->GetParError(i)*mult);
    double lowlim, higlim;
    f->GetParLimits(i, lowlim, higlim);

    const bool limits = !(lowlim == 0 && higlim == 0);

    if(limits && ans[i] < lowlim) ans[i] = lowlim;
  }
  return ans;
}

vector<double> pickpoint(const vector<double> & bestfit,
                         const vector<double> & lowbounds,
                         const vector<double> & higbounds)
{
  vector<double> ans;
  for(unsigned int i = 0; i < bestfit.size(); i++)
    ans.push_back(lowbounds[i] + r.Rndm()*(higbounds[i]+lowbounds[i]));
  return ans; 
}

// Estimate the error on the fit using a MC method
pair<double,double> mcerror(TF1 * fin, const double low,
                            const double hig)
{
  vector<double> bestfit = getbestfit(fin);
  vector<double> lowbounds = getlowbounds(fin);
  vector<double> higbounds = gethigbounds(fin);

  for(unsigned int i = 0; i < lowbounds.size(); i++)
    printf("%d bounds: %f, %f\n", i, lowbounds[i], higbounds[i]);

  const double bestintegral = fin->Integral(low, hig);

  double bestlike = 0;
  int npar = fin->GetNpar();
  gMinuit->Eval(npar, NULL, bestlike, &(bestfit[0]), 0);

  TH1D * soup = new TH1D("soup", "", 1000, 0, bestintegral*2);

  TF1 * f = (TF1 *)fin->Clone("fcopy");

  const int maxtrial = 1000000;
  const int desired  =     100;
  vector<double> vals;
  double like = 0;
  for(int i = 0; i < maxtrial || (puts("bored now") && 0); i++){
    vals = pickpoint(bestfit, lowbounds, higbounds);
    gMinuit->Eval(npar, NULL, like, &(vals[0]), 0);
    const double rat = exp(-like+bestlike);
    if(r.Rndm() < rat){
      f->SetParameters(&(vals[0]));
      soup->Fill(f->Integral(low, hig));
      if(int(soup->GetEntries())%(desired/80) == 0){
        printf("."); fflush(stdout);
      }
      if(soup->GetEntries() > desired) break;
    }
  }

  new TCanvas;

  soup->Draw("e");

  const int midbin = soup->GetXaxis()->FindBin(bestintegral);
  const double halfmidval = soup->GetBinContent(midbin)/2;

  const double above = soup->Integral(midbin+1, soup->GetNbinsX())
                     + halfmidval;
  int higlimbin = midbin+1;
  for(; higlimbin < soup->GetNbinsX(); higlimbin++)
    if((halfmidval+soup->Integral(midbin+1, higlimbin))/above > 0.6827)
      break;
  higlimbin--;

  const double below = soup->Integral(1, midbin-1)
                     + soup->GetBinContent(midbin)/2;
  int lowlimbin = midbin-1;
  for(; lowlimbin > 0; lowlimbin--)
    if((halfmidval+soup->Integral(lowlimbin, midbin-1))/below > 0.6827)
      break;
  lowlimbin++;

  const double anslow = (midbin-lowlimbin)*soup->GetXaxis()->GetBinWidth(1);
  const double anshig = (higlimbin-midbin)*soup->GetXaxis()->GetBinWidth(1);

  TMarker * m = new TMarker(bestintegral, 0, kFullCircle);
  m->SetMarkerSize(1.2);

  TArrow * a = new TArrow(bestintegral-anslow, 0,
                          bestintegral+anshig, 0, 0.02, "<|>");

  soup->Fit("gaus", "l", "e");
  m->Draw();
  a->Draw();

  return pair<double,double>(anslow/bestintegral, anshig/bestintegral);
}

void trigeff_neutron_finalfit(const double lowe = 70,
                              const double highe = 215,
                              const double lowt = 20)
{
  if(!t) t = readfile();

  TF1 * two = makef();

  TCanvas * c1 = new TCanvas("c1", "c1");

  TH1D * h = new TH1D("h", "", 1600, 0, 800);
  t->Draw("t/1000 >> h",
    Form("!ov && fqiv > 225e3 && fqiv < 400e3 && "
         "fq/8300 > %lf && fq/8300 < %lf && e > 2 && e < 2.6 &&"
         "nn > 1" // crucial for rejecting 
                  // stopping muons (in IV?) as through-going +
                  // Michel decay (in buffer?) as neutron
         , lowe, highe),
    "e");

  h->Fit("two", "l", "", lowt, 800); 

  two->ReleaseParameter(0);
  h->Fit("two", "l", "", lowt, 800); 

  two->ReleaseParameter(1);
  two->SetParLimits(1, 24, 30); // capture Gd
  h->Fit("two", "le", "", lowt, 800); 

  two->ReleaseParameter(2);
  two->SetParLimits(2, 0, 0.73); // Thermalization, fig3b of doc-4028
  h->Fit("two", "le", "", lowt, 800); 

  gMinuit->Command("show min");

  h->GetXaxis()->SetRangeUser(0, 50);

  TF1 * ext = (TF1 *)two->Clone("ext");
  ext->SetLineStyle(kDashed);
  ext->Draw("Same");

  const double eff = h->Integral(7,10)*h->GetBinWidth(1)/two->Integral(3,5);
  const double feffed = 1/sqrt(h->Integral(7,10));
  const pair<double,double> feffem = mcerror(two, 3, 5); // low, high
  const double effe_up=eff*sqrt(pow(feffed,2)+pow(feffem.first, 2)); // flipped!
  const double effe_lo=eff*sqrt(pow(feffed,2)+pow(feffem.second,2));

  TGraphAsymmErrors * gmcerr = new TGraphAsymmErrors;

  gmcerr->SetMarkerColor(kRed);
  gmcerr->SetLineColor(kRed);
  gmcerr->SetLineWidth(3);
  gmcerr->SetMarkerStyle(kFullCircle);
  gmcerr->SetPoint(0, 4.0, two->Integral(3,5)/2);
  gmcerr->SetPointError(0, 0, 0,
                        feffem.first*two->Integral(3,5)/2,
                        feffem.second*two->Integral(3,5)/2);

  gmcerr->SetName("gmcerr");
  c1->cd();
  gmcerr->Draw("p");

  gErrorIgnoreLevel = kWarning;

  c1->SaveAs(Form("trigeff-%.0f.C", lowe/2+highe/2));

  printf("g->SetPoint(i,%3.0f, %.4f); "
         "g->SetPointError(i++,xe,xe, %.4f, %.4f);\n",
         lowe/2+highe/2, eff, effe_lo, effe_up);
}
