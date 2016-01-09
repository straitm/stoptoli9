#include "consts.h"
#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TH1.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TMinuit.h"

TH1D * mt = NULL;

static bool far = true;

double ff(double * ta, double * par)
{
  const double mupeff = par[4], mumeff = par[4] - 0.0039;

  double t = *ta;

  const double mulife_ns = mulife*1e6;

  return par[0]*mt->GetBinWidth(1)*(
      mupeff*(1-par[1])/mulife_ns * exp(-t/mulife_ns)
    + mumeff*par[1]*(par[3]/mulife_ns-0.0003 /* correction for non-C-12 */)
                          /(par[7]+par[8]+par[9]*par[10])
                          *(par[7]/par[3] * exp(-t/par[3])
                           +par[8]/par[6] * exp(-t/par[6])
                           +par[9]/par[5] * exp(-t/par[5])
                           +par[10]/lifetime_gd*1e-6*exp(-t/lifetime_gd*1e-6)))
       + fabs(par[11]); // accidental background
}

void fcn(int & npar, double * gin, double & like, double *par, int flag)
{
  like = 0;
  
  for(int i = (far?40:71); i <= (far?78:192); i++){
    const double data = mt->GetBinContent(i);
    const double t0 = mt->GetBinLowEdge(i)-par[2];
    const double t1 = mt->GetBinLowEdge(i+1)-par[2];
    const int nint = 100;
    double mc = 0;
    for(int j = 0; j <= nint; j++){
      double t = t0 + (t1 - t0)*j/double(nint);
      mc += ff(&t, par)/(nint + 1);
    }

    like += mc - data + data*log(data/mc);
  }
  like *= 2;

  // pulls
  //like += pow((par[1] - 0.4409)/0.0032, 2);
  like += pow((par[2] - (-12))/6, 2);
  like += pow((par[3] - lifetime_c12*1e6)/lifetime_c12_err*1e-6, 2);
  like += pow((par[4] - 0.9888)/0.003, 2);
  like += pow((par[5] - lifetime_o16*1e6)/lifetime_o16_err*1e-6, 2);
  like += pow((par[6] - lifetime_c13*1e6)/lifetime_c13_err*1e-6, 2);
  like += pow((par[9] - 0.002)/0.001, 2);
  like += pow((par[10] - 0.00022)/0.00007, 2);
}

void chargeratio_finalfit(const bool far_ = true)
{
  far=far_;
  if(far){
    mt = new TH1D("mt", "", 90, 0, 5760);
    mt->SetBinContent(3,65);
    mt->SetBinContent(4,2984);
    mt->SetBinContent(5,20636);
    mt->SetBinContent(6,29609);
    mt->SetBinContent(7,29440);
    mt->SetBinContent(8,20978);
    mt->SetBinContent(9,15129);
    mt->SetBinContent(10,16941);
    mt->SetBinContent(11,21610);
    mt->SetBinContent(12,27522);
    mt->SetBinContent(13,30980);
    mt->SetBinContent(14,29539);
    mt->SetBinContent(15,27178);
    mt->SetBinContent(16,25689);
    mt->SetBinContent(17,26164);
    mt->SetBinContent(18,26085);
    mt->SetBinContent(19,26367);
    mt->SetBinContent(20,26247);
    mt->SetBinContent(21,24760);
    mt->SetBinContent(22,24078);
    mt->SetBinContent(23,23158);
    mt->SetBinContent(24,22481);
    mt->SetBinContent(25,22369);
    mt->SetBinContent(26,22082);
    mt->SetBinContent(27,21547);
    mt->SetBinContent(28,20657);
    mt->SetBinContent(29,19765);
    mt->SetBinContent(30,19097);
    mt->SetBinContent(31,18546);
    mt->SetBinContent(32,18170);
    mt->SetBinContent(33,17546);
    mt->SetBinContent(34,17181);
    mt->SetBinContent(35,16508);
    mt->SetBinContent(36,15845);
    mt->SetBinContent(37,15348);
    mt->SetBinContent(38,14826);
    mt->SetBinContent(39,14609);
    mt->SetBinContent(40,14059);
    mt->SetBinContent(41,13608);
    mt->SetBinContent(42,13166);
    mt->SetBinContent(43,12840);
    mt->SetBinContent(44,12535);
    mt->SetBinContent(45,12129);
    mt->SetBinContent(46,11698);
    mt->SetBinContent(47,11372);
    mt->SetBinContent(48,10793);
    mt->SetBinContent(49,10746);
    mt->SetBinContent(50,10429);
    mt->SetBinContent(51,9870);
    mt->SetBinContent(52,9797);
    mt->SetBinContent(53,9616);
    mt->SetBinContent(54,9136);
    mt->SetBinContent(55,8745);
    mt->SetBinContent(56,8738);
    mt->SetBinContent(57,8379);
    mt->SetBinContent(58,8147);
    mt->SetBinContent(59,7745);
    mt->SetBinContent(60,7554);
    mt->SetBinContent(61,7471);
    mt->SetBinContent(62,7160);
    mt->SetBinContent(63,7014);
    mt->SetBinContent(64,6768);
    mt->SetBinContent(65,6483);
    mt->SetBinContent(66,6380);
    mt->SetBinContent(67,6364);
    mt->SetBinContent(68,5964);
    mt->SetBinContent(69,5956);
    mt->SetBinContent(70,5764);
    mt->SetBinContent(71,5543);
    mt->SetBinContent(72,5296);
    mt->SetBinContent(73,5145);
    mt->SetBinContent(74,5127);
    mt->SetBinContent(75,4971);
    mt->SetBinContent(76,4672);
    mt->SetBinContent(77,4605);
    mt->SetBinContent(78,4610);
    mt->SetBinContent(79,4369);
    mt->SetBinContent(80,4354);
    mt->SetBinContent(81,4473);
    mt->SetBinContent(82,4566);
    mt->SetBinContent(83,4160);
    mt->SetBinContent(84,2702);
    mt->SetBinContent(85,1494);
    mt->SetBinContent(86,696);
    mt->SetBinContent(87,5);
  }
  else{
   mt = new TH1D("mt", "", 192, 0, 12288);
   const char * const rootfilenearmich = "/cp/s4/strait/ndfido/nearmich.root";
   TFile * tf = new TFile(rootfilenearmich, "read");
   if(!tf || tf->IsZombie()){
     fprintf(stderr, "Could not open %s\n", rootfilenearmich);
     exit(1);
   }
   TTree * tt = dynamic_cast<TTree *>(tf->Get("t"));
   if(!tt){
     fprintf(stderr, "Could not get t TTree from %s\n", rootfilenearmich);
     exit(1);
   }

   tt->Draw("micht >> mkt(192, 0, 12288)", "ndecay == 0 && miche > 30 && miche < 60");
           //"&& rchi2 < 4 && mx**2+my**2 < 1050**2 && mz > -1175");
   TH1F * mkt = (TH1F *)gROOT->FindObject("mkt");
   for(int i = 1; i <= mkt->GetNbinsX(); i++)
     mt->SetBinContent(i, mkt->GetBinContent(i));
   printf("%.0f michels\n", mt->Integral());
  }


  mt->Draw("e");

  const unsigned int npar = 12;
  TMinuit * mn = new TMinuit(npar);
  mn->SetFCN(fcn);

  int err;
  mn->mnparm(0, "nmuon", 1.6e6, 1e4, 0, 0, err);
  mn->mnparm(1, "mumfrac", 0.44, 0.01, 0, 1, err);
  mn->mnparm(2, "michtoff", 12, 0.01, 0, 0, err);
  mn->mnparm(3, "c12time", lifetime_c12*1e6, lifetime_c12_err*1e6, 0, 0, err);
  mn->mnparm(4, "michpluseff", 0.9888, 0.01, 0, 1, err);
  mn->mnparm(5, "o16time", lifetime_o16*1e6, lifetime_o16_err*1e6, 0, 0, err);
  mn->mnparm(6, "c13time", lifetime_c13*1e6, lifetime_c13_err*1e6, 0, 0, err);
  mn->mnparm(7, "c12afrac", 0.98767, 0.01, 0, 1, err);
  mn->mnparm(8, "c13afrac", 0.01088, 0.01, 0, 1, err);
  mn->mnparm(9, "o16afrac", 0.002, 0.01, 0, 0.1, err);
  mn->mnparm(10, "gdafrac", 0.0002, 0.01, 0, 1, err);
  mn->mnparm(11, "acc",      10, 0.01, 0, 1000, err);

  mn->Command("SET STRATEGY 2");

  mn->Command("FIX 2");
  mn->Command("FIX 3");
  mn->Command("FIX 4");
  mn->Command("FIX 5");
  mn->Command("FIX 6");
  mn->Command("FIX 7");
  mn->Command("FIX 11");

  mn->Command("FIX 8");
  mn->Command("FIX 9");
  
  mn->Command("MIGRAD");

  mn->Command("REL 2");
  mn->Command("REL 3");
  mn->Command("REL 4");
  mn->Command("REL 5");
  mn->Command("REL 6");
  mn->Command("REL 7");
  mn->Command("REL 11");

  while(mn->Command("MIGRAD") == 4);

  mn->Command("MINOS");
  mn->Command("SHOW MINOS");

  TF1 * result = new TF1("result", ff, 0, 100000, npar);
  for(unsigned int i = 0; i < npar; i++){
    double par, e;
    mn->GetParameter(i, par, e);
    result->SetParameter(i, par);
  }
  result->Draw("same");

  double par, e;
  mn->GetParameter(0, par, e);

  printf("number of mu-: %f %f +%f\n", par, mn->fErn[0], mn->fErp[0]);
  if(far){
    const double nexpected = 1628874.;
    printf("percent of expected: %.2f %.2f +%.2f\n",
            par*100/nexpected, mn->fErn[0]*100/nexpected,
                               mn->fErp[0]*100/nexpected);
  }

}
