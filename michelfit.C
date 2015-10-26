void michelfit()
{
  TH1D * mt = new TH1D("mt", "", 80, 0, 5120);
  t->Draw("micht >> mt",
    "ndecay == 0 && mz > -1175 && mx**2+my**2 < 1050**2 && miche > 20 && miche < 80 && abs(fez + 62*ivdedx/2 - 8847.2) < 1000 && rchi2 < 2", "e")

TF1 * mf2 = new TF1("mf2", "[0]*(3*(x/52.72)^2 - 2*(x/52.72)^3)*((1-(x/52.72))/(3*x/52.72/7))^(2/137.036/3.14159 *(log(105.658*x/52.72/0.510999)-1))", 0, 52.72);
  mf2->SetParameter(0, 1);

  const double minusfrac = 1; res= 0.21, escale= 1.01, norm = 10000000, carbonshift = 0; carbonsmear = 1.5;
  foo->Reset();
  for(int i = 0; i < norm; i++) {
    bool plus = gRandom->Rndm()>minusfrac;
    double truee = (mf2->GetRandom() - 0.511 + plus*1.022 + (!plus)*carbonshift + (!plus)*gRandom->Gaus(carbonsmear))*escale;
    allminus->Fill(truee + sqrt(truee)*gRandom->Gaus(0, res));
  }

  const double mupeff = 0.9888, mumeff = 0.9849;

  TF1 * ee = new TF1("ee", Form("[0]*%f*((1-[1])/2196.9803 * exp(-x/2196.9803) + [1]*%f*0.9228*(0.98767/2028. * exp(-x/2028.)+0.01088/2037. * exp(-x/2037.)+0.00118/1796. * exp(-x/1796.)+0.00027/80.58 * exp(-x/80.58)))", mupeff, mumeff), 0, 1e5);
  ee->SetParLimits(1, 0, 1);


  mt->Fit("ee", "liem", "", 2000, 5000);

}
