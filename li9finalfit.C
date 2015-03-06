#include "consts.h"

const char * const RED     = "\033[31;1m"; // bold red
const char * const CLR      = "\033[m"    ; // clear

const int npar = 10;

const float dist = 300;

// bn decay probability multiplied by the number of captures relative
// to C-12
const double li9ebn = 0.5080,
             he8ebn = 0.1600,
             c16ebn = 0.99  * 88.0/102.5*0.00243*n_o16cap_betan/n_c12cap,
             n17ebn = 0.951 * 88.0/102.5*0.00243*n_o16cap_betan/n_c12cap,
             b13ebn = 0.0029 * n_c13cap/n_c12cap,
             li11ebn = 0.789 * n_c13cap/n_c12cap;
             
const double distcuteffgc = 0.7487,
             distcutefftarg = 0.9202;

// 100s begin-of-run requirement taken into account here
const double denominator = 0.9709*livetime*n_c12cap;

/* DC3rdPub product of muon, light noise, OV, multiplicity,
   neutron (E, t, R), FV and IV efficiencies */
const double Geff_sans_prompt_or_mun =
       distcutefftarg*
       (1-4.49/100.)*
       (1-0.01/100.)*
       (1-0.06/100.)*
       (1-1.06/100.)*
       0.9829*
       (1-0.66/100.)*
       (1-0.04/100.);

const double Heff_sans_prompt_or_mun =
  (
  (n_c12cap - n_c12captarget)      * distcuteffgc
+ n_c12captarget * (1-gd_fraction) * distcutefftarg
  )/  
  (n_c12cap - gd_fraction*n_c12captarget)* 
       (1-1.25*4.49/100.)* // muon - ok, straightforwards scaling
       (1-0.01/100.)* // ? LN - same cut, but not obviously same eff
                      // however, *very* small for Gd, so...
       (1-0.06/100.)* // OV - ok, same
       (1-((8.+9.)/(2.+6.))*1.06/100.)* // Mult - I think this is valid
       0.9512*       // neutron (ANN,E,t,R) -- plot on slide 6 of 
                    // doc-5863 -- I think this is right
       (1-0.806/100.)* // FV -- doc5480
       (1-0.025/100.)* // IV prompt -- doc5813
       (1-0.014/100.); // IV delayed -- doc5813


// But not the actual gd fraction because of geometrical effects
const double expectedgdfrac = gd_fraction*n_c12captarget/n_c12cap;

int whichh = 0;
void drawhist(TTree * tgsel, TTree * thsel,
              const vector< vector<double> > & parsave, 
              const int nbin, const double low, const double high)
{
  whichh++;

  TCanvas * c1 = new TCanvas;
//  c1->Divide(1, 3);
//  c1->cd(1);

  tgsel->Draw(Form("dt/1000 >> hdispg%d(%d,%f,%f)", whichh,nbin,low,high));
  thsel->Draw(Form("dt/1000 >> hdisph%d(%d,%f,%f)", whichh,nbin,low,high));
  thsel->Draw(Form("dt/1000 >> hdisp%d (%d,%f,%f)", whichh,nbin,low,high));
  tgsel->Draw(Form("dt/1000 >> +hdisp%d", whichh), "", "e");

  TH1D * hdispg = gROOT->FindObject(Form("hdispg%d", whichh));
  TH1D * hdisph = gROOT->FindObject(Form("hdisph%d", whichh));
  TH1D * hdisp  = gROOT->FindObject(Form("hdisp%d",  whichh));
  TH1D * hists[3] = { hdispg, hdisph, hdisp };
  //for(int i = 0; i < 3; i++) hists[i]->SetLineWidth(1);

  TF1 * eedisps[parsave.size()];
  for(int i = 0; i < parsave.size(); i++){
    eedisps[i] = new TF1(Form("eedisp%d-%d", i, whichh),
      Form("%f*"
      "([10]*%f*("
         "%f*[2]*(1-[9])/0.257233*exp(-x/0.257233)+"  // H Li-9
         "%f*[3]*(1-[9])/0.171825*exp(-x/0.171825)+"  // H He-8
         "%f*[4]*(1-[5])/6.020366*exp(-x/6.020366)+"  // H N-17
         "%f*[6]*(1-[5])/1.077693*exp(-x/1.077693)+"  // H C-16
         "%f*[7]*(1-[9])/0.025002*exp(-x/0.025002)+"  // H B-13
         "%f*[8]*(1-[9])/0.012624*exp(-x/0.012624)+"  // H Li-11
         "[0]*(1-[1])) + "                            // H bg
       "[11]*%f*("
         "%f*[2]*[9]/0.257233*exp(-x/0.257233)+" // Gd Li-9
         "%f*[3]*[9]/0.171825*exp(-x/0.171825)+" // Gd He-8
         "%f*[4]*[5]/6.020366*exp(-x/6.020366)+" // Gd N-17
         "%f*[6]*[5]/1.077693*exp(-x/1.077693)+" // Gd C-16
         "%f*[7]*[9]/0.025002*exp(-x/0.025002)+" // Gd B-13
         "%f*[8]*[9]/0.012624*exp(-x/0.012624)+" // Gd Li-11
         "[0]*[1]))",                            // Gd bg
           denominator,
           Heff, li9ebn, he8ebn, n17ebn, c16ebn, b13ebn, li11ebn,
           Geff, li9ebn, he8ebn, n17ebn, c16ebn, b13ebn, li11ebn), 0, 100);
    for(int j = 0; j < eedisps[i]->GetNpar(); j++)
      eedisps[i]->SetParameter(j, parsave[i][j]);
    eedisps[i]->SetLineWidth(2);
  }

  for(int i = 0; i < parsave.size(); i++){
    eedisps[i]->SetParameter(npar, hdisp->GetBinWidth(1));
    eedisps[i]->SetParameter(npar+1, hdisp->GetBinWidth(1));
    eedisps[i]->SetLineColor(kRed);
    eedisps[i]->SetNpx(400);
    eedisps[i]->Draw("Same");
  }
/*
  c1->cd(2);
  hdisph->Draw("e");
  for(int i = 0; i < parsave.size(); i++){
    TF1 * eedisph = (TF1*)eedisps[i]->Clone(Form("eedisph%d-%d", i, whichh));
    eedisph->SetParameter(npar, hdisp->GetBinWidth(1));
    eedisph->SetParameter(npar+1, 0);
    eedisph->Draw("same");
  }

  c1->cd(3);
  hdispg->Draw("e");
  for(int i = 0; i < parsave.size(); i++){
    TF1 * eedispg = (TF1*)eedisps[i]->Clone(Form("eedispg%d-%d", i, whichh));
    eedispg->SetParameter(npar, 0);
    eedispg->SetParameter(npar+1, hdisp->GetBinWidth(1));
    eedispg->Draw("same");
  }
*/
}

int whichc = -1;
TCanvas * cans[100]; // oh no
void contour(TMinuit * mn, const int par1, const int par2,
             const double xrange, const double yrange, const int points,
             const char * const comment)
{
  whichc++;
  TCanvas * c = new TCanvas(Form("c%d_%d%d", whichc, par1, par2),
                            Form("c%d_%d%d", whichc, par1, par2),
                            600, 350);
  mn->Command("MIGRAD");
  const double minx = getpar(mn, par1-1);
  const double miny = getpar(mn, par2-1);

  const int oldprintlevel = mn->fISW[4];
  mn->Command("Set print 0");

  mn->fUp = 2.3; // 68% in 2D
  mn->Command(Form("mncont %d %d %d", par1, par2, points));
  TGraph * sigma_2d =
    mn->GetPlot()?(TGraph*)((TGraph*)mn->GetPlot())->Clone():NULL;

/*
  mn->fUp = 2.61; // 90% in 1D
  mn->Command(Form("mncont %d %d %d", par1, par2, points));
  TGraph * ninty_1d =
    mn->GetPlot()?(TGraph*)((TGraph*)mn->GetPlot())->Clone():NULL;
*/

  mn->fUp = 4.61; // 90%
  mn->Command(Form("mncont %d %d %d", par1, par2, points));
  TGraph * ninty_2d =
    mn->GetPlot()?(TGraph*)((TGraph*)mn->GetPlot())->Clone():NULL;


  if(ninty_2d){
    ninty_2d->SetFillColor(kViolet);
    ninty_2d->Draw("alf");
    ninty_2d->GetXaxis()->SetRangeUser(0, xrange);
    ninty_2d->GetYaxis()->SetRangeUser(0, yrange);
    ninty_2d->GetXaxis()->SetTitle(mn->fCpnam[par1-1]);
    ninty_2d->GetYaxis()->SetTitle(mn->fCpnam[par2-1]);
    ((TGaxis*)(ninty_2d->GetXaxis()))->SetMaxDigits(3);
    ((TGaxis*)(ninty_2d->GetYaxis()))->SetMaxDigits(3);
  }
  //if(ninty_1d) ninty_1d->SetLineColor(kRed),   ninty_1d->Draw("l");
  if(sigma_2d) sigma_2d->SetLineColor(kBlack), sigma_2d->Draw("l");
/*  mn->fUp = 11.83; // 99.73% in 2D
  mn->Command(Form("mncont %d %d %d", par1, par2, points));
  TGraph * threesigma_2d =
    mn->GetPlot()?(TGraph*)((TGraph*)mn->GetPlot())->Clone():NULL;

  if(threesigma_2d) threesigma_2d->SetLineStyle(kDashed),
                    threesigma_2d->Draw("l"); */

  TMarker * best = new TMarker(minx, miny, kStar);
  best->Draw();

  if(comment != ""){
    TLatex * t = new TLatex(0.5, 0.8, comment);
    t->SetTextSize(0.07);
    t->SetNDC(1);
    t->Draw();
  }

  cans[whichc] = c;
  for(int i = 0; i <= whichc; i++){
    cans[i]->Modified();
    cans[i]->Update();
  }

  mn->Command(Form("set print %d", oldprintlevel));
}

void setupmn(TMinuit * mn, const double expectedgdfrac)
{
  for(int i = 0; i < npar; i++) mn->Command(Form("REL %d", i+1));
  mn->Command("SET LIM 1 0 1e-3");

  mn->Command(Form("SET PAR 2 %f", expectedgdfrac));
  mn->Command("Set LIM 2 0 0.5");

  mn->Command("SET LIM 3 0 1e-3");
  mn->Command("SET LIM 4 0 3e-3");
  mn->Command("SET LIM 5 0 10"); // allow 100% production plus
                                 // a large factor for slop.  There's
                                 // a pull term, too.
  mn->Command("SET LIM 6 0 0.5");
  mn->Command("SET LIM 7 0 1");  // as with N-17
  mn->Command("SET LIM 8 0 3");
  mn->Command("SET LIM 9 0 1");

  mn->Command(Form("SET PAR 10 %f", expectedgdfrac));
  mn->Command("Set LIM 10 0 0.5");
}

vector<double> gtimes, htimes;
double Geff, Heff;

void fcn(int & npar, double * gin, double & like, double *par, int flag)
{
  like = (
        denominator*Heff*(
           li9ebn*par[2]*(1-par[9])+// H Li-9
           he8ebn*par[3]*(1-par[9])+// H He-8
           n17ebn*par[4]*(1-par[5])+// H N-17
           c16ebn*par[6]*(1-par[5])+// H C-16
           b13ebn*par[7]*(1-par[9])+// H B-13
          li11ebn*par[8]*(1-par[9])+// H Li-11
           99.999*par[0]*(1-par[1])) + // H bg
         denominator*Geff*(
           li9ebn*par[2]*par[9]+ // Gd Li-9
           he8ebn*par[3]*par[9]+ // Gd He-8
           n17ebn*par[4]*par[5]+ // Gd N-17
           c16ebn*par[6]*par[5]+ // Gd C-16
           b13ebn*par[7]*par[9]+ // Gd B-13
          li11ebn*par[8]*par[9]+ // Gd Li-11
           99.999*par[0]*par[1]));     // Gd bg

  vector<double> * v[2] = { &gtimes, &htimes };
  for(int j = 0; j < 2; j++){
    for(int i = 0; i < v[j]->size(); i++){
      const double x = (*v[j])[i];
  
      double f = j == 1?
        Heff*(
           li9ebn*par[2]*(1-par[9])/0.257233*exp(-x/0.257233)+// H Li-9
           he8ebn*par[3]*(1-par[9])/0.171825*exp(-x/0.171825)+// H He-8
           n17ebn*par[4]*(1-par[5])/6.020366*exp(-x/6.020366)+// H N-17
           c16ebn*par[6]*(1-par[5])/1.077693*exp(-x/1.077693)+// H C-16
           b13ebn*par[7]*(1-par[9])/0.025002*exp(-x/0.025002)+// H B-13
          li11ebn*par[8]*(1-par[9])/0.012624*exp(-x/0.012624)+// H Li-11
           par[0]*(1-par[1])) :                               // H bg
         Geff*(
           li9ebn*par[2]*par[9]/0.257233*exp(-x/0.257233)+ // Gd Li-9
           he8ebn*par[3]*par[9]/0.171825*exp(-x/0.171825)+ // Gd He-8
           n17ebn*par[4]*par[5]/6.020366*exp(-x/6.020366)+ // Gd N-17
           c16ebn*par[6]*par[5]/1.077693*exp(-x/1.077693)+ // Gd C-16
           b13ebn*par[7]*par[9]/0.025002*exp(-x/0.025002)+ // Gd B-13
          li11ebn*par[8]*par[9]/0.012624*exp(-x/0.012624)+ // Gd Li-11
           par[0]*par[1]);                                 // Gd bg

      if(f > 0) like += -log(f);
    }
  }

  like *= 2;

  // pull terms
  like += pow((par[4]-0.5)/0.5, 2); // 50%+-70% for N-17
  like += pow((par[6]-0.05)/0.05, 2); // 5%+-10%  for C-16
}

double getpar(TMinuit * mn, int i)
{
  double answer, dum;
  mn->GetParameter(i, answer, dum);
  return answer;
}

void fixat(TMinuit * mn, int i, float v)
{
  mn->Command(Form("SET LIM %d", i));
  mn->Command(Form("SET PAR %d %f", i, v));
  mn->Command(Form("FIX %d", i));
}

void fixatzero(TMinuit * mn, int i)
{
  mn->Command(Form("SET LIM %d", i));
  mn->Command(Form("SET PAR %d 0", i));
  mn->Command(Form("FIX %d", i));
}

TGraph * lastgraphstyle(TCanvas * can, const int which)
{
  TIter next(can->GetListOfPrimitives());
  TObject * obj = NULL;
  TGraph * tomakered = NULL;
  while(obj=next()) if(dynamic_cast<TGraph*>obj) tomakered=(TGraph*)obj;
  if(!tomakered) return;
  if(which == 0){
    tomakered->SetMarkerColor(kRed);
    tomakered->SetMarkerStyle(kFullDotLarge);
  }
  if(which == 1){
    tomakered->SetMarkerColor(kViolet);
    tomakered->SetMarkerStyle(kFullDotLarge);
  }
  return tomakered;
}

void li9finalfit(int neutrons = -1, int contourmask = 0)
{
  // First factor takes into account the efficiency of selecting a
  // neutron after a muon, second is the prompt energy cut, third as
  // documented above
  // 
  // Neutron efficiencies are assuming any within the Michel window
  // are rejected.
  Geff=(neutrons > 0?pow(neff_dr_1000_gd*neff_dt_targ, neutrons):1)
    *0.996*Geff_sans_prompt_or_mun,
  Heff=(neutrons > 0?pow(neff_dr_1000_h*neff_dt_h, neutrons):1)
    *0.993*Heff_sans_prompt_or_mun;

  ///////////////////////////////////////////////////////////////////
 
  char nodistcut[1000];
  snprintf(nodistcut, 999,
    "%smiche < 12 && !earlymich && prompttime > 100e9 && dt < 100e3",
           neutrons >= 0?Form("nlate==%d&&", neutrons):"");
  char cut[1000];
  snprintf(cut, 999, "%s && dist < %f", nodistcut, dist);

  ////////////////////////////////////////////////////////////////////

  TTree tg("t", "t");
  TTree th("t", "t");
  tg.ReadFile("/cp/s4/strait/li9ntuples/li9-20150219.Gd.ntuple");
  th.ReadFile("/cp/s4/strait/li9ntuples/li9-20150219.H.ntuple");

  TTree * tgsel = tg.CopyTree(cut);
  TTree * thsel = th.CopyTree(cut);


  const double plotsigtime = 0.4;

/*
  TCanvas * dxy = new TCanvas("dxy", "dxy", 200, 200);
  th.Draw("my-dy:mx-dx", Form("%s && dt/1000 > 5  ", cut), ".");
  tg.Draw("my-dy:mx-dx", Form("%s && dt/1000 > 5  ", cut), ".same");
  th.Draw("my-dy:mx-dx", Form("%s && dt/1000 < %f", cut, plotsigtime),"%same");
  lastgraphstyle(dxy, 0);
  tg.Draw("my-dy:mx-dx", Form("%s && dt/1000 < %f", cut, plotsigtime),"%same");
  lastgraphstyle(dxy, 0);

  TCanvas * dxz = new TCanvas("dxz", "dxz", 200, 200);
  th.Draw("mz-dz:mx-dx", Form("%s && dt/1000 > 5  ", cut), ".");
  tg.Draw("mz-dz:mx-dx", Form("%s && dt/1000 > 5  ", cut), ".same");
  th.Draw("mz-dz:mx-dx", Form("%s && dt/1000 < %f", cut, plotsigtime),"%same");
  lastgraphstyle(dxz, 0);
  tg.Draw("mz-dz:mx-dx", Form("%s && dt/1000 < %f", cut, plotsigtime),"%same");
  lastgraphstyle(dxz, 0);

  TCanvas * dyz = new TCanvas("dyz", "dyz", 200, 200);
  th.Draw("mz-dz:my-dy", Form("%s && dt/1000 > 5  ", cut), ".");
  tg.Draw("mz-dz:my-dy", Form("%s && dt/1000 > 5  ", cut), ".same");
  th.Draw("mz-dz:my-dy", Form("%s && dt/1000 < %f", cut, plotsigtime),"%same");
  lastgraphstyle(dyz, 0);
  tg.Draw("mz-dz:my-dy", Form("%s && dt/1000 < %f", cut, plotsigtime),"%same");
  lastgraphstyle(dyz, 0);

  TCanvas * muenergy = new TCanvas("muenergy", "muenergy", 200, 200);
  th.Draw("dedxslant >> muqivbg(25, 0, 12000)", Form("%s && dt/1000 > 5  ", cut), "hist");
  tg.Draw("dedxslant >> +muqivbg", Form("%s && dt/1000 > 5  ", cut), "esame");
  th.Draw("dedxslant >> muqivsig(25, 0, 12000)", Form("%s && dt/1000 < 0.4", cut),"same");
  tg.Draw("dedxslant >> +muqivsig", Form("%s && dt/1000 < 0.4", cut),"same");
*/
  float tim;
  tgsel->SetBranchAddress("dt", &tim);
  thsel->SetBranchAddress("dt", &tim);

  for(int i = 0; i < tgsel->GetEntries(); i++){
    tgsel->GetEntry(i);
    gtimes.push_back((tim - 0.002)/1000);
  }
  for(int i = 0; i < thsel->GetEntries(); i++){
    thsel->GetEntry(i);
    htimes.push_back((tim - 0.002)/1000);
  }

  //////////////////

  TMinuit * mn = new TMinuit(npar);
  mn->SetFCN(fcn);
  int err;
  mn->mnparm(1 -1, "bg",     1e-4,  5e-6, 0, 0, err);
  mn->mnparm(2 -1, "gdfracacc",0.38,0.01, 0, 0, err);
  mn->mnparm(3 -1, "Li-9",   1e-5,  1e-6, 0, 0, err); // n: yes
  mn->mnparm(4 -1, "He-8",    0.1,  1e-3, 0, 0, err); // n: yes
  mn->mnparm(5 -1, "N-17",    0.1,   0.5, 0, 0, err); // n: yes
  mn->mnparm(6 -1,"ocapgdfrac",0.01,1e-2, 0, 0, err);
  mn->mnparm(7 -1, "C-16",   0.01,   0.5, 0, 0, err); // n: yes
  mn->mnparm(8 -1, "B-13",   0.01,     1, 0, 0, err); // n: no
  mn->mnparm(9 -1, "Li-11",  0.01,  3e-3, 0, 0, err); // n: no
  mn->mnparm(10-1, "gdfracsig",0.38, 0.01,0, 0, err);

  mn->Command("set print -10");

  vector< vector<double> > parsaves;
  {
    setupmn(mn, expectedgdfrac);
    const unsigned int nfix = 7;
    const int fix[nfix] = { 3, 4, 5, 6, 7, 8, 9 };
    for(int i = 0; i < nfix; i++){
      mn->Command(Form("SET LIM %d", fix[i]));
      mn->Command(Form("SET PAR %d 0", fix[i]));
      mn->Command(Form("FIX %d", fix[i]));
    }
    mn->Command("MIGRAD");
    mn->Command("MINOS 10000 3");
    mn->Command("SHOW min");
    printf("%sNo Li-9 or other bn isotopes (%.2f)%s\n",
           RED, mn->fAmin, CLR);
    vector<double> parsave;
    for(int i = 0; i < npar; i++) parsave.push_back(getpar(mn, i));
    //parsaves.push_back(parsave);
  }
  const double chi2nothing = mn->fAmin;

  {
    setupmn(mn, expectedgdfrac);
    const unsigned int nfix = 6;
    const int fix[nfix] = { 4, 5, 6, 7, 8, 9 };
    for(int i = 0; i < nfix; i++){
      mn->Command(Form("SET LIM %d", fix[i]));
      mn->Command(Form("SET PAR %d 0", fix[i]));
      mn->Command(Form("FIX %d", fix[i]));
    }
    mn->Command("MIGRAD");
    mn->Command("MINOS 10000 3");
    mn->Command("SHOW min");
    printf("%sLi-9 prob without other isotopes (%.2f): %f %f +%f%s\n",
           RED, mn->fAmin, getpar(mn, 2), mn->fErn[2], mn->fErp[2], CLR);
    vector<double> parsave;
    for(int i = 0; i < npar; i++) parsave.push_back(getpar(mn, i));
    //parsaves.push_back(parsave);
  }
  const double chi2justli9 = mn->fAmin;

  {
    setupmn(mn, expectedgdfrac);
    const unsigned int nfix = 5;
    const int fix[nfix] = { 4, 6, 5, 8, 9 };
    for(int i = 0; i < nfix; i++){
      mn->Command(Form("SET LIM %d", fix[i]));
      mn->Command(Form("SET PAR %d 0", fix[i]));
      mn->Command(Form("FIX %d", fix[i]));
    }
    mn->Command("MIGRAD");
    mn->Command("MINOS 10000 3");
    mn->Command("SHOW min");
    printf("%sLi-9 with C-16 only (%.2f): %f %f +%f%s\n", RED,
      mn->fAmin, getpar(mn, 2), mn->fErn[2], mn->fErp[2], CLR);
  }
  const double chi2_li9_c16 = mn->fAmin;

  {
    setupmn(mn, expectedgdfrac);
    const unsigned int nfix = 5;
    const int fix[nfix] = { 4, 6, 7, 8, 9 };
    for(int i = 0; i < nfix; i++){
      mn->Command(Form("SET LIM %d", fix[i]));
      mn->Command(Form("SET PAR %d 0", fix[i]));
      mn->Command(Form("FIX %d", fix[i]));
    }
    mn->Command("MIGRAD");
    mn->Command("MINOS 10000 3");
    mn->Command("SHOW min");
    printf("%sLi-9 with N-17 only (%.2f): %f %f +%f%s\n", RED,
      mn->fAmin, getpar(mn, 2), mn->fErn[2], mn->fErp[2], CLR);
  }
  const double chi2_li9_n17 = mn->fAmin;

  {
    setupmn(mn, expectedgdfrac);
    const unsigned int nfix = 5;
    const int fix[nfix] = { 5, 6, 7, 8, 9 };
    for(int i = 0; i < nfix; i++){
      mn->Command(Form("SET LIM %d", fix[i]));
      mn->Command(Form("SET PAR %d 0", fix[i]));
      mn->Command(Form("FIX %d", fix[i]));
    }
    mn->Command("MIGRAD");
    mn->Command("MINOS 10000 3");
    mn->Command("SHOW min");
    printf("%sLi-9 prob with only He-8 (%.2f): %f %f +%f%s\n", RED,
      mn->fAmin, getpar(mn, 2), mn->fErn[2], mn->fErp[2], CLR);
  }
  const double chi2_li9_he8 = mn->fAmin;

  {
    setupmn(mn, expectedgdfrac);
    const unsigned int nfix = 1;
    const int fix[nfix] = { 4 };
    for(int i = 0; i < nfix; i++){
      mn->Command(Form("SET LIM %d", fix[i]));
      mn->Command(Form("SET PAR %d 0", fix[i]));
      mn->Command(Form("FIX %d", fix[i]));
    }
    if(neutrons > 0){
      const unsigned int nfixn = 2;
      const int fixn[nfixn] = { 8, 9 };
      for(int i = 0; i < nfixn; i++){
        mn->Command(Form("SET LIM %d", fixn[i]));
        mn->Command(Form("SET PAR %d 0", fixn[i]));
        mn->Command(Form("FIX %d", fixn[i]));
      }
    }
    mn->Command("MIGRAD");
    mn->Command("MIGRAD");
    mn->Command("MINOS 10000 3 5");
    mn->Command("SHOW min");
    printf("%sLi-9 prob/C-12 capture with everything except "
           "He-8 (%.2f): %f %f +%f%s\n",
           RED, mn->fAmin, getpar(mn, 2), mn->fErn[2], mn->fErp[2], CLR);
    printf("%sLi-9 prob/C-13 capture with everything except "
           "He-8 (%.2f): %f %f +%f%s\n",
           RED, mn->fAmin, getpar(mn, 2)*n_c12cap/n_c13cap,
           mn->fErn[2]*n_c12cap/n_c13cap, mn->fErp[2]*n_c12cap/n_c13cap,
           CLR);
    printf("%sLi-9 prob/C-nat capture with everything except "
           "He-8 (%.2f): %f %f +%f%s\n",
           RED, mn->fAmin, getpar(mn, 2)*n_c12cap/(n_c12cap+n_c13cap),
           mn->fErn[2]*n_c12cap/(n_c12cap+n_c13cap),
           mn->fErp[2]*n_c12cap/(n_c12cap+n_c13cap),
           CLR);
  }
  const double chi2_allbut_he8 = mn->fAmin;

  {
    setupmn(mn, expectedgdfrac);
    if(neutrons > 0){
      const unsigned int nfixn = 2;
      const int fixn[nfixn] = { 8, 9 };
      for(int i = 0; i < nfixn; i++){
        mn->Command(Form("SET LIM %d", fixn[i]));
        mn->Command(Form("SET PAR %d 0", fixn[i]));
        mn->Command(Form("FIX %d", fixn[i]));
      }
    }
    mn->Command("MIGRAD");
    mn->Command("MINOS 10000 3 5 8");
    mn->Command("SHOW min");
    printf("%sLi-9 prob with all other isotopes (%.2f): %f %f +%f%s\n",
      RED, mn->fAmin, getpar(mn, 2), mn->fErn[2], mn->fErp[2], CLR);
    vector<double> parsave;
    for(int i = 0; i < npar; i++) parsave.push_back(getpar(mn, i));
    parsaves.push_back(parsave);
  }
  const double chi2_all = mn->fAmin;

  const double maxb13 = getpar(mn, 7) + mn->fErp[7] * sqrt(2.3);
  printf("%s90%% upper limit B-13: %f\n%s", RED, maxb13, CLR);

  if(neutrons < 1){
    setupmn(mn, expectedgdfrac);
    fixat(mn, 8, maxb13);
    mn->Command("MIGRAD");
    mn->Command("MINOS 10000 3 5");
    mn->Command("SHOW min");
    printf("%sLi-9 prob w/everything, B-13 fixed at limit (%.2f): %f %f +%f%s\n", RED,
      mn->fAmin, getpar(mn, 2), mn->fErn[2], mn->fErp[2], CLR);
    vector<double> parsave;
    for(int i = 0; i < npar; i++) parsave.push_back(getpar(mn, i));
    parsaves.push_back(parsave);
  }

  printf("%s", RED);
  printf("Li-9 preferred over nothing by %f\n",
         chi2nothing - chi2justli9 < 0?0:sqrt(chi2nothing - chi2justli9));

  printf("Li-9 + C-16 preferred over just Li-9 by %f\n",
         chi2justli9 - chi2_li9_c16<0?0:sqrt(chi2justli9 - chi2_li9_c16));

  printf("Li-9 + N-17 preferred over just Li-9 by %f\n",
         chi2justli9 - chi2_li9_n17<0?0:sqrt(chi2justli9 - chi2_li9_n17));

  printf("All preferred over Li-9 + N-17 by %f\n",
         chi2_li9_n17 - chi2_all<0?0:sqrt(chi2_li9_n17 - chi2_all));
  printf("All preferred over nothing by %f\n",
         chi2nothing - chi2_all<0?0:sqrt(chi2nothing - chi2_all));
  printf("%s", CLR);

  const int npoint = 100;

  // Li-9 vs. N-17 with nothing else
  if(contourmask & 0x01){
    setupmn(mn, expectedgdfrac);
    fixatzero(mn, 4);
    fixatzero(mn, 7);
    fixatzero(mn, 8);
    fixatzero(mn, 9);
    contour(mn, 3, 5,  0.0079, 2, npoint, "Nothing else");
  }

  // Li-9 vs. C-16 with no He-8
  if(contourmask & 0x02){
    setupmn(mn, expectedgdfrac);
    fixatzero(mn, 4);
    contour(mn, 3, 7, 0.0079, 2, npoint, "No ^{8}He");
  }

  // Li-9 vs. B-13 with no Li-11 or He-8
  if(contourmask & 0x04){
    setupmn(mn, expectedgdfrac);
    fixatzero(mn, 4);
    fixatzero(mn, 9);
    contour(mn, 3, 8,  0.0079, 2, npoint, "No ^{8}He or ^{11}Li");
  }

  // Li-9 vs. He-8
  if(contourmask & 0x08){
    setupmn(mn, expectedgdfrac);
    contour(mn, 3, 4, 0.0079, 0.00149, npoint, "");
  }

  // Li-9 vs. He-8 with nothing else
  if(contourmask & 0x10){
    setupmn(mn, expectedgdfrac);
    fixatzero(mn, 5);
    fixatzero(mn, 6);
    fixatzero(mn, 7);
    fixatzero(mn, 8);
    fixatzero(mn, 9);
    contour(mn, 3, 4, 0.0079, 0.00149, npoint, "Nothing else");
  }

  // B-13 vs. Li-11
  if(contourmask & 0x20){
    setupmn(mn, expectedgdfrac);
    contour(mn, 8, 9, 3, 0.008, npoint, "");
  }

  
  //////////////////////////////////////////////////////////////////////
  drawhist(tgsel, thsel, parsaves, 48, 1, 97);
  drawhist(tgsel, thsel, parsaves, 100, 0, 20);
  drawhist(tgsel, thsel, parsaves, 10, 0.001, 0.501);


  setupmn(mn, expectedgdfrac);
  for(int i = 0; i < 10; i++){
    mn->Command(Form("SET PAR 8 %f", i*0.1));
    mn->Command("FIX 8");
    mn->Command("MIGRAD");
    printf("%f: %f\n", i*0.1, mn->fAmin); 
  }

/*
  TCanvas * xy = new TCanvas("xy", "xy", 200, 200);
  th.Draw("dy/1000:dx/1000", Form("%s && dt/1000 > 10 && dt/1000 < 30", cut), ".");
  tg.Draw("dy/1000:dx/1000", Form("%s && dt/1000 > 10 && dt/1000 < 30", cut), ".same");
  th.Draw("dy/1000:dx/1000", Form("%s && dt/1000 < 0.4", cut),"%same");
  lastgraphstyle(xy, 0);
  tg.Draw("dy/1000:dx/1000", Form("%s && dt/1000 < 0.4", cut),"%same");
  lastgraphstyle(xy, 0);

  TCanvas * rz = new TCanvas("rz", "rz", 200, 200);
  th.Draw("dz:(dx**2+dy**2)", Form("%s && dt/1000 > 10 && dt/1000 < 30", cut), ".");
  tg.Draw("dz:(dx**2+dy**2)", Form("%s && dt/1000 > 10 && dt/1000 < 30", cut), ".same");
  th.Draw("dz:(dx**2+dy**2)", Form("%s && dt/1000 < 0.4", cut),"%same");
  lastgraphstyle(rz, 0);
  tg.Draw("dz:(dx**2+dy**2)", Form("%s && dt/1000 < 0.4", cut),"%same");
  lastgraphstyle(rz, 0);

  TCanvas * dr = new TCanvas("dr", "dr", 200, 200);
  th.Draw("(dist**3)/1e9 >>  distbg(25,0,0.5)",Form("%s && dt/1000 > 5",  nodistcut), "hist");
  th.Draw("(dist**3)/1e9 >> distsig(25,0,0.5)",Form("%s && dt/1000 < 0.4",nodistcut),"esame");
  tg.Draw("(dist**3)/1e9 >>  +distbg",           Form("%s && dt/1000 > 5",  nodistcut), "histsame");
  tg.Draw("(dist**3)/1e9 >> +distsig",           Form("%s && dt/1000 < 0.4",nodistcut),"esame");
  distsig=distsig;
  distbg=distbg;
  distbg->Sumw2();
  distbg->Scale(distsig->Integral(5, 100)/distbg->Integral(5, 100));
  distsig->SetLineColor(kRed);
  distsig->SetMarkerColor(kRed);
  distbg->GetYaxis()->SetRangeUser(0, 25);
*/

/*
  TCanvas * rzh = new TCanvas("rzh", "rzh", 200, 200);
  th.Draw("dz >> dzall(10, -1800, 1800)", cut, "hist");
  th.Draw("dz >> dzn16", Form("%s && dt/1000 > 0.75 && dt/1000 < 5", cut),"samee");
  tg.Draw("dz >> +dzall", cut, "hist");
  tg.Draw("dz >> +dzn16", Form("%s && dt/1000 > 0.75 && dt/1000 < 5", cut),"samee");
  dzall=dzall;
  dzn16=dzn16;
  dzn16->Scale(dzall->Integral()/dzn16->Integral());
*/

/*
  printf("Candidates:\n");
  th.SetScanField(0);
  tg.SetScanField(0);
  th.Scan("run:trig:dt:sqrt(dx*dx+dy*dy):dx:dy:dz", cut);
  tg.Scan("run:trig:dt:sqrt(dx*dx+dy*dy):dx:dy:dz", cut);
*/
}
