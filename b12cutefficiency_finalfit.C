#include "consts.h"

void b12cutefficiency_finalfit()
{
  TFile *_file0 = TFile::Open(rootfile0up, "Read");
  TH1D * offtime = new TH1D("offtime", "", 32, 0, 16);
  TH1D * b12e    = new TH1D("b12e",    "", 32, 0, 16);
  t->Draw("e >>    b12e", "       dt > 1 && dt < 201 && dist < 200 && miche < 5 && e > 0 && timeleft > 10e3", "hist");
  t->Draw("e >> offtime", "0.05*(abs(dt-8000) < 2000 && dist < 200 && miche < 5 && e > 0 && timeleft > 10e3)", "histsame");
  offtime->Sumw2();
  b12e->Sumw2();

  TH1F * b12bgsubed = (TH1F *)b12e->Clone("b12bgsubed");
  b12bgsubed->Sumw2();
  b12bgsubed->Add(offtime, -1);
  b12bgsubed->Draw("samee");

  TH1D *mc = new TH1D("mc","",64,0,16);
  mc->SetBinContent(0,8.98654e-05);
  mc->SetBinContent(1,0.4772014);
  mc->SetBinContent(2,1.273401);
  mc->SetBinContent(3,2.165104);
  mc->SetBinContent(4,3.173662);
  mc->SetBinContent(5,4.336578);
  mc->SetBinContent(6,5.599399);
  mc->SetBinContent(7,6.955147);
  mc->SetBinContent(8,8.28891);
  mc->SetBinContent(9,9.774658);
  mc->SetBinContent(10,11.32163);
  mc->SetBinContent(11,12.69009);
  mc->SetBinContent(12,14.25593);
  mc->SetBinContent(13,15.72911);
  mc->SetBinContent(14,17.12121);
  mc->SetBinContent(15,18.57431);
  mc->SetBinContent(16,19.91585);
  mc->SetBinContent(17,21.14464);
  mc->SetBinContent(18,22.40733);
  mc->SetBinContent(19,23.50755);
  mc->SetBinContent(20,24.53792);
  mc->SetBinContent(21,25.38844);
  mc->SetBinContent(22,26.1637);
  mc->SetBinContent(23,26.96112);
  mc->SetBinContent(24,27.5278);
  mc->SetBinContent(25,28.01203);
  mc->SetBinContent(26,28.32374);
  mc->SetBinContent(27,28.56479);
  mc->SetBinContent(28,28.61508);
  mc->SetBinContent(29,28.6225);
  mc->SetBinContent(30,28.44862);
  mc->SetBinContent(31,28.04431);
  mc->SetBinContent(32,27.59352);
  mc->SetBinContent(33,26.94085);
  mc->SetBinContent(34,26.362);
  mc->SetBinContent(35,25.53488);
  mc->SetBinContent(36,24.79255);
  mc->SetBinContent(37,23.60384);
  mc->SetBinContent(38,22.56014);
  mc->SetBinContent(39,21.45009);
  mc->SetBinContent(40,20.06593);
  mc->SetBinContent(41,18.82833);
  mc->SetBinContent(42,17.31324);
  mc->SetBinContent(43,15.81148);
  mc->SetBinContent(44,14.3762);
  mc->SetBinContent(45,12.86794);
  mc->SetBinContent(46,11.30934);
  mc->SetBinContent(47,9.775316);
  mc->SetBinContent(48,8.350178);
  mc->SetBinContent(49,6.909903);
  mc->SetBinContent(50,5.589053);
  mc->SetBinContent(51,4.310231);
  mc->SetBinContent(52,3.216128);
  mc->SetBinContent(53,2.248103);
  mc->SetBinContent(54,1.449251);
  mc->SetBinContent(55,0.842222);
  mc->SetBinContent(56,0.4434239);
  mc->SetBinContent(57,0.1873616);
  mc->SetBinContent(58,0.06912553);
  mc->SetBinContent(59,0.02313549);
  mc->SetBinContent(60,0.004590106);
  mc->SetBinContent(61,0.0004501984);
  mc->SetBinContent(62,4.357642e-07);
  mc->SetBinContent(63,8.98654e-05);
  mc->SetBinContent(64,4.357642e-07);
  mc->Rebin(2);

  mc->Scale(b12bgsubed->Integral(2, 32)/mc->Integral(2,32));

  mc->Draw("samehist");
  mc->SetLineColor(kRed);
  offtime->SetLineColor(kGreen+2);

  b12bgsubed->SetBinContent(1, mc->GetBinContent(1));
  b12bgsubed->SetBinError(1, mc->GetBinContent(1)/4.);

  const double cut = 4;

  double bb[3] = {0, cut, 16};
  TH1D * cuth = b12bgsubed->Rebin(2, "cut", bb);
  const double eff = cuth->GetBinContent(2)/cuth->Integral(1,2);
  const double error = sqrt((cuth->GetBinContent(2)*cuth->GetBinError(2)/(cuth->Integral(1,2)**2))**2
                          + (cuth->GetBinContent(1)*cuth->GetBinError(1)/(cuth->Integral(1,2)**2))**2);
  printf("Efficiency for %.1f MeV cut for B-12+B-13: %.2f +- %.2f %%\n",
         cut, eff*100, error*100);
};
