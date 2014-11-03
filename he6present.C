{
TCanvas c;
c.Divide(1,2, 0.01, 0);

  int r = 2;
  ehistbg_1->Rebin(r);
  ehistbg_2->Rebin(r);
  ehistbg_3->Rebin(r);
  ehistsig_1->Rebin(r);
  ehistsig_2->Rebin(r);
  ehistsig_3->Rebin(r);
  bg_plus_li8_n16_1->Rebin(r);
  bg_plus_li8_n16_2->Rebin(r);
  bg_plus_li8_n16_3->Rebin(r);
  bg_plus_li8_n16_he6_1->Rebin(r);
  bg_plus_li8_n16_he6_2->Rebin(r);
  bg_plus_li8_n16_he6_3->Rebin(r);

  ehistbg_1->Add(ehistbg_2);
  ehistsig_1->Add(ehistsig_2);
  bg_plus_li8_n16_1->Add(bg_plus_li8_n16_2);
  bg_plus_li8_n16_he6_1->Add(bg_plus_li8_n16_he6_2);

/*
  ehistsig_1->Add(ehistsig_3);
  bg_plus_li8_n16_1->Add(bg_plus_li8_n16_3);
  bg_plus_li8_n16_he6_1->Add(bg_plus_li8_n16_he6_3); 
*/

  c.cd(1)->SetLogy();
  c.cd(1)->SetPad(0, 0.4, 0.98, 0.99);
  ((TCanvas *)c.cd(1))->SetLeftMargin(0.14);

  bg_plus_li8_n16_he6_1->GetXaxis()->SetRange(3, bg_plus_li8_n16_he6_1->GetNbinsX());
  bg_plus_li8_n16_he6_1->GetYaxis()->SetRangeUser(0.05, 700);
  bg_plus_li8_n16_he6_1->Draw("hist");
  bg_plus_li8_n16_1->Draw("histsame");
  ehistbg_1->Draw("histsame");
  ehistsig_1->Draw("esame");

  const double size = 0.11;
  bg_plus_li8_n16_he6_1->GetYaxis()->SetTitleSize(size * 2./3);
  bg_plus_li8_n16_he6_1->GetYaxis()->SetLabelSize(size * 2./3);
  bg_plus_li8_n16_he6_1->GetYaxis()->SetTitle("Events/0.25 MeV");
  bg_plus_li8_n16_he6_1->GetYaxis()->CenterTitle();
  

  c.cd(2)->SetPad(0, 0, 0.98, 0.4);
  ((TCanvas *)c.cd(2))->SetBottomMargin(0.24);
  ((TCanvas *)c.cd(2))->SetLeftMargin(0.14);

  TH1D * fit = (TH1D*)bg_plus_li8_n16_he6_1->Clone("fit");
  TH1D * sig = (TH1D*)ehistsig_1->Clone("sig");

  fit->Add(bg_plus_li8_n16_1, -1);
  sig->Add(bg_plus_li8_n16_1, -1);
  sig->GetXaxis()->SetRange(3, ehistsig_1->GetNbinsX());
  sig->GetYaxis()->SetRangeUser(-24, 24);
  sig->Draw("e");
  fit->Draw("histsame][");

  sig->GetXaxis()->SetTitle("Energy (MeV)");
  sig->GetYaxis()->SetTitle("signal - bg");
  sig->GetYaxis()->SetTitleOffset(0.6667);

  sig->GetXaxis()->CenterTitle();
  sig->GetYaxis()->CenterTitle();

  sig->GetXaxis()->SetTitleSize(size);
  sig->GetYaxis()->SetTitleSize(size);
  sig->GetXaxis()->SetLabelSize(size);
  sig->GetYaxis()->SetLabelSize(size);

  sig->GetYaxis()->SetNdivisions(505);

  TF1 f("f","0", 0, 15);
  f.SetNpx(100);
  f.Draw("same");
}
