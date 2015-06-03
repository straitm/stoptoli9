#include "consts.h"

TH1D * mt = new TH1D("mt", "", 80, 0, 5120);

double ff(double * ta, double * par)
{
  const double mupeff = 0.9888, mumeff = 0.9849;

  double t = *ta;

  return par[0]*(
      mupeff*(1-par[1])/2196.9811 * exp(-t/2196.9811)
    + mumeff*  par[1] *0.9228*(0.98767/2028. * exp(-t/2028.)
                              +0.01088/2037. * exp(-t/2037.)
                              +0.00118/1796. * exp(-t/1796.)
                              +0.00027/80.58 * exp(-t/80.58)));
}

void fcn(int & npar, double * gin, double & like, double *par, int flag)
{
  like = 0;
  
  for(int i = 40; i <= 78; i++){
    const double data = mt->GetBinContent(i);
    const double t0 = mt->GetBinLowEdge(i);
    const double t1 = mt->GetBinLowEdge(i+1);
    const int nint = 10;
    double mc = 0;
    for(int j = 0; j <= nint; j++){
      const double t = t0 + (t1 - t0)*j/double(nint);
      mc += ff(&t, par)/(nint + 1);
    }

    like += mc - data + data*log(data/mc);
  }
  like *= 2;

  //like += 
}

void chargeratiofinalfit()
{
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


  mt->Draw("e");

  TMinuit * mn = new TMinuit(2);
  mn->SetPrintLevel(-1);
  mn->SetFCN(fcn);

  int err;
  mn->mnparm(0, "nmuon", 1e6, 1e4, 0, 0, err);
  mn->mnparm(1, "mumfrac", 0.44, 0.01, 0, 1, err);

  mn->Command("MIGRAD");
  mn->Command("MIGRAD");
  mn->Command("MINOS");
  mn->Command("SHOW MINOS");

  TF1 * result = new TF1("result", ff, 0, 10000, 2);
  for(int i = 0; i < 2; i++){
    double par, e;
    mn->GetParameter(i, par, e);
    result->SetParameter(i, par);
  }
  result->Draw("same");

  double par, e;
  mn->GetParameter(0, par, e);

  printf("number of mu-: %f +- %f\n", par/mt->GetBinWidth(1), e/mt->GetBinWidth(1));

}
