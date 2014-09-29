{
  TTree tg("t", "t");
  TTree th("t", "t");
  tg.ReadFile("li9-20140925.Gd.ntuple");
  th.ReadFile("li9-20140925.H.ntuple");

  tg.Draw("dt/1000 >> hfitg(10000, 0.001, 10)", 
          "dist < 300 && miche < 12");
  tg.Draw("dt/1000 >> hdispg(50, 0.001, 10.001)", 
          "dist < 300 && miche < 12");

  TF1 ee("ee", "[0]*exp(-x*log(2)/[1]) + [2]*exp(-x*log(2)/[3]) + [4]", 0, 10);
  ee.FixParameter(1, 0.1783);
  ee.FixParameter(3, 0.1191);
  ee.FixParameter(2, 0);

  hfitg->Fit("ee", "le", "");

  TF1 e("e", "[0]*exp(-x*log(2)/[1])", 0, 10);
  for(int i = 0; i < 2; i++)
    e->SetParameter(i, ee->GetParameter(i));

  printf("Nraw = %.2f\n", e->Integral(0, 10)/hfitg.GetBinWidth(1));
  printf("Nerr+ = %.4f\n", e->Integral(0, 10)/hfitg.GetBinWidth(1) / ee.GetParameter(0) * gMinuit.fErp[0]);
  printf("Nerr- = %.4f\n", e->Integral(0, 10)/hfitg.GetBinWidth(1) / ee.GetParameter(0) * gMinuit.fErn[0]);

  for(int i = 0; i <= 4; i+=2)
    ee.SetParameter(i, ee->GetParameter(i)*hdispg.GetBinWidth(1)/
                                           hfitg.GetBinWidth(1));

  hdispg->Draw("e");

  ee.Draw("Same");


  th.Draw("dt/1000 >> hfitg(10000, 0.001, 10)", 
          "dist < 300 && miche < 12");
  th.Draw("dt/1000 >> hdispg(50, 0.001, 10.001)", 
          "dist < 300 && miche < 12");

  TF1 ee("ee", "[0]*exp(-x*log(2)/[1]) + [2]*exp(-x*log(2)/[3]) + [4]", 0, 10);
  ee.FixParameter(1, 0.1783);
  ee.FixParameter(3, 0.1191);
  ee.FixParameter(2, 0);

  hfitg->Fit("ee", "le", "");

  TF1 e("e", "[0]*exp(-x*log(2)/[1])", 0, 10);
  for(int i = 0; i < 2; i++)
    e->SetParameter(i, ee->GetParameter(i));

  printf("Nraw = %.2f\n", e->Integral(0, 10)/hfitg.GetBinWidth(1));
  printf("Nerr+ = %.4f\n", e->Integral(0, 10)/hfitg.GetBinWidth(1) / ee.GetParameter(0) * gMinuit.fErp[0]);
  printf("Nerr- = %.4f\n", e->Integral(0, 10)/hfitg.GetBinWidth(1) / ee.GetParameter(0) * gMinuit.fErn[0]);

  for(int i = 0; i <= 4; i+=2)
    ee.SetParameter(i, ee->GetParameter(i)*hdispg.GetBinWidth(1)/
                                           hfitg.GetBinWidth(1));

  hdispg->Draw("e");

  ee.Draw("Same");


  th.Draw("dt/1000 >> hfitg(10000, 0.001, 10)", 
          "dist < 300 && miche < 12");
  th.Draw("dt/1000 >> hdispg(50, 0.001, 10.001)", 
          "dist < 300 && miche < 12");
  tg.Draw("dt/1000 >> +hfitg", 
          "dist < 300 && miche < 12");
  tg.Draw("dt/1000 >> +hdispg", 
          "dist < 300 && miche < 12");

  TF1 ee("ee", "[0]*exp(-x*log(2)/[1]) + [2]*exp(-x*log(2)/[3]) + [4]", 0, 10);
  ee.FixParameter(1, 0.1783);
  ee.FixParameter(3, 0.1191);
  ee.FixParameter(2, 0);

  hfitg->Fit("ee", "le", "");

  TF1 e("e", "[0]*exp(-x*log(2)/[1])", 0, 10);
  for(int i = 0; i < 2; i++)
    e->SetParameter(i, ee->GetParameter(i));

  printf("Nraw = %.2f\n", e->Integral(0, 10)/hfitg.GetBinWidth(1));
  printf("Nerr+ = %.4f\n", e->Integral(0, 10)/hfitg.GetBinWidth(1) / ee.GetParameter(0) * gMinuit.fErp[0]);
  printf("Nerr- = %.4f\n", e->Integral(0, 10)/hfitg.GetBinWidth(1) / ee.GetParameter(0) * gMinuit.fErn[0]);

  for(int i = 0; i <= 4; i+=2)
    ee.SetParameter(i, ee->GetParameter(i)*hdispg.GetBinWidth(1)/
                                           hfitg.GetBinWidth(1));

  hdispg->Draw("e");

  ee.Draw("Same");

}
