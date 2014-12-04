const int npar = 9;

const float dist = 300;

const double livetime = 489.509, nstop = 367.;
const double nstoptarg = 139.0;

// bn decay probability multiplied by the number of captures relative
// to C-12
const double li9ebn = 0.5080,
             he8ebn = 0.1600,
             c16ebn = 0.99  * 88.0/102.5*0.00243*7.0/nstop,
             n17ebn = 0.951 * 88.0/102.5*0.00243*7.0/nstop,
             b13ebn = 0.0029 * 3.7/nstop,
             li11ebn = 0.789 * 3.7/nstop;
             
const double distcuteff = (dist == 400?0.948:dist == 300?0.852:
                           dist == 200?0.565:dist==159?0.376:100000);

// 30s begin-of-run requirement taken into account here
const double denominator = 0.99127*livetime*nstop*distcuteff;

const double gdcapfrac = 0.871;

/* DC3rdPub product of muon, light noise, OV, multiplicity,
   neutron (E, t, R), FV and IV efficiencies */
const double Geff_sans_prompt_or_mun = (1-4.49/100.)*
       (1-0.01/100.)*
       (1-0.06/100.)*
       (1-1.06/100.)*
       0.9829*
       (1-0.66/100.)*
       (1-0.04/100.);

const double Heff_sans_prompt_or_mun =
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


const double expectedgdfrac = gdcapfrac*nstoptarg/nstop;

int whichh = 0;
void drawhist(TTree * tgsel, TTree * thsel,
              const vector< vector<double> > & parsave, const double high,
              const int nbin)
{
  whichh++;

  TCanvas * c1 = new TCanvas;
  c1->Divide(1, 3);
  c1->cd(1);

  tgsel->Draw(Form("dt/1000 >> hdispg%d(%d, 0, %f)", whichh, nbin, high));
  thsel->Draw(Form("dt/1000 >> hdisph%d(%d, 0, %f)", whichh, nbin, high));
  thsel->Draw(Form("dt/1000 >> hdisp%d (%d, 0, %f)", whichh, nbin, high));
  tgsel->Draw(Form("dt/1000 >> +hdisp%d", whichh), "", "e");

  TH1D * hdispg = gROOT->FindObject(Form("hdispg%d", whichh));
  TH1D * hdisph = gROOT->FindObject(Form("hdisph%d", whichh));
  TH1D * hdisp  = gROOT->FindObject(Form("hdisp%d",  whichh));

  TF1 * eedisps[parsave.size()];
  for(int i = 0; i < parsave.size(); i++){
    TF1 * eedisp = new TF1(Form("eedisp%d-%d", i, whichh),
      Form("%f*"
      "([9]*%f*("
         "%f*[2]*(1-[1])/0.257233*exp(-x/0.257233)+"  // H Li-9
         "%f*[3]*(1-[1])/0.171825*exp(-x/0.171825)+"  // H He-8
         "%f*[4]*(1-[5])/6.020366*exp(-x/6.020366)+"  // H N-17
         "%f*[6]*(1-[5])/1.077693*exp(-x/1.077693)+"  // H C-16
         "%f*[7]*(1-[1])/0.025002*exp(-x/0.025002)+"  // H B-13
         "%f*[8]*(1-[1])/0.012624*exp(-x/0.012624)+"  // H Li-11
         "[0]*(1-[1])) + "                            // H bg
       "[10]*%f*("
         "%f*[2]*[1]/0.257233*exp(-x/0.257233)+" // Gd Li-9
         "%f*[3]*[1]/0.171825*exp(-x/0.171825)+" // Gd He-8
         "%f*[4]*[5]/6.020366*exp(-x/6.020366)+" // Gd N-17
         "%f*[6]*[5]/1.077693*exp(-x/1.077693)+" // Gd C-16
         "%f*[7]*[1]/0.025002*exp(-x/0.025002)+" // Gd B-13
         "%f*[8]*[1]/0.012624*exp(-x/0.012624)+" // Gd Li-11
         "[0]*[1]))",                                              // Gd bg
           denominator,
           Heff, li9ebn, he8ebn, n17ebn, c16ebn, b13ebn, li11ebn,
           Geff, li9ebn, he8ebn, n17ebn, c16ebn, b13ebn, li11ebn), 0, 30);
    for(int j = 0; j < eedisp->GetNpar(); j++)
      eedisp->SetParameter(j, parsave[i][j]);
    eedisps[i] = eedisp;
  }

  for(int i = 0; i < parsave.size(); i++){
    TF1 * eedisp = eedisps[i];
    eedisp->SetParameter(npar, hdisp->GetBinWidth(1));
    eedisp->SetParameter(npar+1, hdisp->GetBinWidth(1));
    eedisp->SetLineColor(kRed);
    eedisp->SetNpx(400);
    eedisp->Draw("Same");
    if(eedisp->Eval(0) > hdisp->GetMaximum() + sqrt(hdisp->GetMaximum()))
      hdisp->GetYaxis()->SetRangeUser(0, eedisp->Eval(0)*1.05);
  }

  c1->cd(2);
  hdisph->Draw("e");
  for(int i = 0; i < parsave.size(); i++){
    TF1 * eedisph = (TF1*)eedisps[i]->Clone(Form("eedisph%d-%d", i, whichh));
    eedisph->SetParameter(npar, hdisp->GetBinWidth(1));
    eedisph->SetParameter(npar+1, 0);
    eedisph->Draw("same");
    if(eedisph->Eval(0) > hdisph->GetMaximum() + sqrt(hdisph->GetMaximum()))
      hdisph->GetYaxis()->SetRangeUser(0, eedisph->Eval(0)*1.05);
  }

  c1->cd(3);
  hdispg->Draw("e");
  for(int i = 0; i < parsave.size(); i++){
    TF1 * eedispg = (TF1*)eedisps[i]->Clone(Form("eedispg%d-%d", i, whichh));
    eedispg->SetParameter(npar, 0);
    eedispg->SetParameter(npar+1, hdisp->GetBinWidth(1));
    eedispg->Draw("same");
    if(eedispg->Eval(0) > hdispg->GetMaximum() + sqrt(hdispg->GetMaximum()))
      hdispg->GetYaxis()->SetRangeUser(0, eedispg->Eval(0)*1.05);
  }
}

void contour(TMinuit * mn, const int par1, const int par2,
             const double xrange, const double yrange, const int points,
             const char * const comment)
{
 return;
  TCanvas * c = new TCanvas(Form("c%d%d", par1, par2),
                            Form("c%d%d", par1, par2), 600, 350);
  mn->Command("MIGRAD");
  const double minx = getpar(mn, par1-1);
  const double miny = getpar(mn, par2-1);

  mn->Command("Set print 0");
  mn->Command("Set strategy 2");

  mn->fUp = 1.0/2; // 68% in 1D
  mn->Command(Form("mncont %d %d %d", par1, par2, points));
  TGraph * sigma_1d =
    mn->GetPlot()?(TGraph*)((TGraph*)mn->GetPlot())->Clone():NULL;

  mn->fUp = 2.3/2; // 90% in 1D
  mn->Command(Form("mncont %d %d %d", par1, par2, points));
  TGraph * ninty_1d =
    mn->GetPlot()?(TGraph*)((TGraph*)mn->GetPlot())->Clone():NULL;

  mn->fUp = 4.61/2; // 90%
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
    ((TGaxis*)(ninty_2d->GetXaxis()))->SetMaxDigits(2);
    ((TGaxis*)(ninty_2d->GetYaxis()))->SetMaxDigits(2);
  }
  if(ninty_1d) ninty_1d->SetLineColor(kRed),   ninty_1d->Draw("l");
  if(sigma_1d) sigma_1d->SetLineColor(kBlack), sigma_1d->Draw("l");

  TMarker * best = new TMarker(minx, miny, kStar);
  best->Draw();

  if(comment != ""){
    TLatex * t = new TLatex(0.5, 0.8, comment);
    t->SetTextSize(0.07);
    t->SetNDC(1);
    t->Draw();
  }

  c->Modified();
  c->Update();
  c->Modified();
  c->Update();
}

void setupmn(TMinuit * mn, const double expectedgdfrac)
{
  for(int i = 0; i < npar; i++) mn->Command(Form("REL %d", i+1));
  mn->Command("SET LIM 1 0 1e-3");

  mn->Command(Form("SET PAR 2 %f", expectedgdfrac));
  mn->Command("FIX 2");

  mn->Command("SET LIM 3 0 1e-3");
  mn->Command("SET LIM 4 0 1e-3");
  mn->Command("SET LIM 5 0 20"); // allow 100% production plus
                                // a factor of two for the O capture
                                // normalization plus a factor of ten
                                // for other complications.
  mn->Command("SET LIM 6 0 0.5");
  mn->Command("SET LIM 7 0 10");
  mn->Command("SET LIM 8 0 3");
  mn->Command("SET LIM 9 0 1");
}

vector<double> gtimes, htimes;
double Geff, Heff;

void fcn(int & npar, double * gin, double & like, double *par, int flag)
{
  like = (
        denominator*Heff*(
           li9ebn*par[2]*(1-par[1])+// H Li-9
           he8ebn*par[3]*(1-par[1])+// H He-8
           n17ebn*par[4]*(1-par[5])+// H N-17
           c16ebn*par[6]*(1-par[5])+// H C-16
           b13ebn*par[7]*(1-par[1])+// H B-13
          li11ebn*par[8]*(1-par[1])+// H Li-11
           29.999*par[0]*(1-par[1])) + // H bg
         denominator*Geff*(
           li9ebn*par[2]*par[1]+ // Gd Li-9
           he8ebn*par[3]*par[1]+ // Gd He-8
           n17ebn*par[4]*par[5]+ // Gd N-17
           c16ebn*par[6]*par[5]+ // Gd C-16
           b13ebn*par[7]*par[1]+ // Gd B-13
          li11ebn*par[8]*par[1]+ // Gd Li-11
           29.999*par[0]*par[1]));     // Gd bg

  //printf("%d %f\n", gtimes.size()+htimes.size(), like);

  vector<double> * v[2] = { &gtimes, &htimes };
  for(int j = 0; j < 2; j++){
    for(int i = 0; i < v[j]->size(); i++){
      const double x = (*v[j])[i];
  
      double f = j == 1?
        Heff*(
           li9ebn*par[2]*(1-par[1])/0.257233*exp(-x/0.257233)+// H Li-9
           he8ebn*par[3]*(1-par[1])/0.171825*exp(-x/0.171825)+// H He-8
           n17ebn*par[4]*(1-par[5])/6.020366*exp(-x/6.020366)+// H N-17
           c16ebn*par[6]*(1-par[5])/1.077693*exp(-x/1.077693)+// H C-16
           b13ebn*par[7]*(1-par[1])/0.025002*exp(-x/0.025002)+// H B-13
          li11ebn*par[8]*(1-par[1])/0.012624*exp(-x/0.012624)+// H Li-11
           par[0]*(1-par[1])) :                               // H bg
         Geff*(
           li9ebn*par[2]*par[1]/0.257233*exp(-x/0.257233)+ // Gd Li-9
           he8ebn*par[3]*par[1]/0.171825*exp(-x/0.171825)+ // Gd He-8
           n17ebn*par[4]*par[5]/6.020366*exp(-x/6.020366)+ // Gd N-17
           c16ebn*par[6]*par[5]/1.077693*exp(-x/1.077693)+ // Gd C-16
           b13ebn*par[7]*par[1]/0.025002*exp(-x/0.025002)+ // Gd B-13
          li11ebn*par[8]*par[1]/0.012624*exp(-x/0.012624)+ // Gd Li-11
           par[0]*par[1]);                                 // Gd bg

      like += -log(f);
    }
  }

  like *= 2;
}

double getpar(TMinuit * mn, int i)
{
  double answer, dum;
  mn->GetParameter(i, answer, dum);
  return answer;
}

void li9finalfit(bool neutron = false)
{
  const char * const RED     = "\033[31;1m"; // bold red
  const char * const CLR      = "\033[m"    ; // clear

  // First factor takes into account the efficiency of selecting a
  // neutron after a muon, second is the prompt energy cut, third as
  // documented above
  Geff=(neutron?0.97*0.64:1)*0.996*Geff_sans_prompt_or_mun,
  Heff=(neutron?0.97*0.93:1)*0.993*Heff_sans_prompt_or_mun;

  ///////////////////////////////////////////////////////////////////
 
  char cut[1000];
  snprintf(cut, 999, "%sdist < %f  && miche < 12 && !earlymich && dt < 30000",
           dist, neutron?"n==1&&":" ");

  ////////////////////////////////////////////////////////////////////

  TTree tg("t", "t");
  TTree th("t", "t");
  tg.ReadFile("li9-20141203.Gd.ntuple");
  th.ReadFile("li9-20141203.H.ntuple");

  TTree * tgsel = tg.CopyTree(cut);
  TTree * thsel = th.CopyTree(cut);

  float tim;
  tgsel->SetBranchAddress("dt", &tim);
  thsel->SetBranchAddress("dt", &tim);

  for(int i = 0; i < tgsel->GetEntries(); i++){
    tgsel->GetEntry(i);
    gtimes.push_back(tim/1000);
  }
  for(int i = 0; i < thsel->GetEntries(); i++){
    thsel->GetEntry(i);
    htimes.push_back(tim/1000);
  }

  //////////////////

  TMinuit * mn = new TMinuit(npar);
  mn->SetFCN(fcn);
  int err;
  mn->mnparm(1 -1, "bg",     1e-4,  5e-6, 0, 0, err);
  mn->mnparm(2 -1, "gdfrac", 0.38,  0.01, 0, 0, err);
  mn->mnparm(3 -1, "Li-9",   1e-5,  1e-6, 0, 0, err);
  mn->mnparm(4 -1, "He-8",    0.1,  1e-3, 0, 0, err);
  mn->mnparm(5 -1, "N-17",    0.1,   0.5, 0, 0, err);
  mn->mnparm(6 -1,"ocapgdfrac",0.01,1e-2, 0, 0, err);
  mn->mnparm(7 -1, "C-16",   0.01,   0.5, 0, 0, err);
  mn->mnparm(8 -1, "B-13",   0.01,     1, 0, 0, err);
  mn->mnparm(9 -1, "Li-11",  0.01,  3e-3, 0, 0, err);

  vector< vector<double> > parsaves;
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
    printf("%sLi-9 prob without other isotopes: %f %f +%f%s\n", RED,
      getpar(mn, 2), mn->fErn[1], mn->fErp[1], CLR);
    vector<double> parsave;
    for(int i = 0; i < npar; i++) parsave.push_back(getpar(mn, i));
    parsaves.push_back(parsave);
  }

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
    printf("%sLi-9 prob with only He-8: %f %f +%f%s\n", RED,
      getpar(mn, 2), mn->fErn[1], mn->fErp[1], CLR);
  }

  {
    setupmn(mn, expectedgdfrac);
    const unsigned int nfix = 1;
    const int fix[nfix] = { 4 };
    for(int i = 0; i < nfix; i++){
      mn->Command(Form("SET LIM %d", fix[i]));
      mn->Command(Form("SET PAR %d 0", fix[i]));
      mn->Command(Form("FIX %d", fix[i]));
    }
    mn->Command("MIGRAD");
    mn->Command("MIGRAD");
    mn->Command("MINOS 10000 3 5");
    mn->Command("SHOW min");
    printf("%sLi-9 prob with nuisance isotopes, but not He-8: %f %f +%f%s\n", RED,
      getpar(mn, 2), mn->fErn[1], mn->fErp[1], CLR);
  }

  {
    setupmn(mn, expectedgdfrac);
    mn->Command("MIGRAD");
    mn->Command("MINOS 10000 3 5");
    mn->Command("SHOW min");
    printf("%sLi-9 prob with all other isotopes: %f %f +%f%s\n", RED,
      getpar(mn, 2), mn->fErn[1], mn->fErp[1], CLR);
    vector<double> parsave;
    for(int i = 0; i < npar; i++) parsave.push_back(getpar(mn, i));
    parsaves.push_back(parsave);
  }

  const int npoint = 40;

  // Li-9 vs. N-17 with no He-8
  setupmn(mn, expectedgdfrac);
  mn->Command("SET PAR 4 0");
  mn->Command("FIX 4");
  contour(mn, 3, 5,  0.0079, 30, npoint, "No He-8");

  // Li-9 vs. C-16 with no He-8
  setupmn(mn, expectedgdfrac);
  mn->Command("SET PAR 4 0");
  mn->Command("FIX 4");
  contour(mn, 3, 7,  0.0079, 10, npoint, "No He-8");

  // Li-9 vs. B-13 with no Li-11 or He-8
  setupmn(mn, expectedgdfrac);
  mn->Command("SET PAR 4 0");
  mn->Command("FIX 4");
  mn->Command("SET PAR 9 0");
  mn->Command("FIX 9");
  contour(mn, 3, 8,  0.0079, 2, npoint, "No He-8 or Li-11");

  // Li-9 vs. He-8
  setupmn(mn, expectedgdfrac);
  contour(mn, 3, 4, 0.0079, 0.00149, npoint, "");

  // B-13 vs. Li-11
  setupmn(mn, expectedgdfrac);
  contour(mn, 8, 9, 3, 0.004, npoint, "");

  
  //////////////////////////////////////////////////////////////////////
 
  drawhist(tgsel, thsel, parsaves, 30, 30);
  drawhist(tgsel, thsel, parsaves, 5, 25);
  drawhist(tgsel, thsel, parsaves, 0.2, 20);

return;
  printf("Candidates:\n");
  th.SetScanField(0);
  tg.SetScanField(0);
  th.Scan("run:trig:dt:sqrt(dx*dx+dy*dy):dx:dy:dz", cut);
  tg.Scan("run:trig:dt:sqrt(dx*dx+dy*dy):dx:dy:dz", cut);
}
