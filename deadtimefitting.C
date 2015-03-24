{
t.Draw("firstlatenneart/1e3:fq/8300 >> h(70, 0, 700, 800, 0, 800)", "latennear == 1 && ndecay == 0", "colz");
TGraph effh, effg, tim;
TF1 perfg("perfg", "[0] + ([1]*exp(-x/189.) + [2]*exp(-x/25.2)/(1+exp(-(x-3.22873)/3.76334)))", 0, 800);
TF1 perfh("perfh", "[0] + ([1]*exp(-x/189.) + [2]*exp(-x/25.2)/(1+exp(-(x-3.22873)/3.76334)))", 0, 800);

TF1 e("e", "[0] + ([1]*exp(-x/189.) + [2]*exp(-x/25.2)/(1+exp(-(x-3.22873)/3.76334)))/(1 + exp(-(x-[4])/[3]))", 0, 800);
TF1 eg("eg", "[0] + ([1]*exp(-x/189.) + [2]*exp(-x/25.2)/(1+exp(-(x-3.22873)/3.76334)))/(1 + exp(-(x-[4])/[3]))", 0, 800);
TF1 eh("eh", "[0] + ([1]*exp(-x/189.) + [2]*exp(-x/25.2)/(1+exp(-(x-3.22873)/3.76334)))/(1 + exp(-(x-[4])/[3]))", 0, 800);

// Alternative version using erf -- gives same result up to 0.1%
// TF1 e("e", "[0] + ([1]*exp(-x/189.) + [2]*exp(-x/25.2)/(1+exp(-(x-3.22873)/3.76334)))*(0.5+0.5*TMath::Erf((x-[4])/[3]))", 0, 800);
// TF1 eg("eg", "[0] + ([1]*exp(-x/189.) + [2]*exp(-x/25.2)/(1+exp(-(x-3.22873)/3.76334)))*(0.5+0.5*TMath::Erf((x-[4])/[3]))", 0, 800);
// TF1 eh("eh", "[0] + ([1]*exp(-x/189.) + [2]*exp(-x/25.2)/(1+exp(-(x-3.22873)/3.76334)))*(0.5+0.5*TMath::Erf((x-[4])/[3]))", 0, 800);

e.SetParLimits(0, 0, 5);
e.SetParLimits(1, 0, 100);
e.SetParLimits(2, 0, 100);
e.SetParLimits(3, 0.3, 5);
e.SetParLimits(4, 6, 70);


for(int b = 66; b >= 12; b--){ h->ProjectionY("p", b, b)->Draw("e"); p->Fit("e", "l", "e", 6, 800); for(int i = 0; i < 5; i++) { eh.SetParameter(i, e->GetParameter(i)); eg.SetParameter(i, e->GetParameter(i)); perfg.SetParameter(i, e.GetParameter(i)); perfh.SetParameter(i, e.GetParameter(i));} for(int i = 0; i < 3; i++){eg.SetParameter(i, 0); eh.SetParameter(i, 0); perfh.SetParameter(i, 0); perfg.SetParameter(i, 0);} eg.SetParameter(2, 1); eh.SetParameter(1, 1); perfg.SetParameter(2, 1); perfh.SetParameter(1, 1); const double tmpeh = eh.Integral(5.5, 800)/perfh.Integral(5.5, 10000), tmpeg = eg.Integral(5.5, 800)/perfg.Integral(5.5, 10000), energy =  h->GetXaxis()->GetBinCenter(b), tmptim = e->GetParameter(4);  printf("%f %.2f %.2f\n", energy, 100*tmpeh, 100*tmpeg); effh.SetPoint(effh.GetN(), energy, tmpeh); effg.SetPoint(effg.GetN(), energy, tmpeg); tim.SetPoint(tim.GetN(), energy, tmptim); c1->Update(); c1->Modified(); } effg->Draw("ap*"); effh->Draw("p%");

for(int i = 0; i < effh->GetN(); i++) effh.GetY()[i] *= exp(-5.5/189);

for(int i = 0; i < effg->GetN(); i++) effg.GetY()[i] *= (1- 0.0726);

TF1 efffitg("efffitg", "[0]*exp(-x/[1]) + gaus(2)", 100, 700);
efffitg->SetParameters(1, 100, 0.04, 190, 30);
TF1 efffith("efffith", "[0]*exp(-x/[1]) + gaus(2)", 100, 700);
efffith->SetParameters(1, 100, 0.04, 190, 30);

efffith->SetParLimits(2, 0, 0.1);
efffitg->SetParLimits(2, 0, 0.1);
efffith->SetParLimits(3, 150, 220);
efffitg->SetParLimits(3, 150, 220);
efffith->SetParLimits(4, 15, 50);
efffitg->SetParLimits(4, 15, 50);

effg->Fit("efffitg", "", "", 100, 550);
effh->Fit("efffith", "", "", 100, 550);

efffitg->SavePrimitive(cout);
efffith->SavePrimitive(cout);

}
