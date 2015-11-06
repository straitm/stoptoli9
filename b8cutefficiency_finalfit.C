#include "TTree.h"
#include "TH1.h"
#include "TF1.h"
#include "TROOT.h"
#include "consts.h"

void b8cutefficiency_finalfit(const double cutlow = 4)
{
  TH1F *h = new TH1F("h","ctEvisID  {ctX[0]**2+ctX[1]**2 < (1150+250)**2 && abs(ctX[2]) < 1229+250}",180,0,18);
   h->SetBinContent(1,104793);
   h->SetBinContent(2,10176);
   h->SetBinContent(3,3250);
   h->SetBinContent(4,1111);
   h->SetBinContent(5,730);
   h->SetBinContent(6,623);
   h->SetBinContent(7,432);
   h->SetBinContent(8,252);
   h->SetBinContent(9,154);
   h->SetBinContent(10,134);
   h->SetBinContent(11,120);
   h->SetBinContent(12,77);
   h->SetBinContent(13,114);
   h->SetBinContent(14,100);
   h->SetBinContent(15,124);
   h->SetBinContent(16,121);
   h->SetBinContent(17,156);
   h->SetBinContent(18,171);
   h->SetBinContent(19,172);
   h->SetBinContent(20,203);
   h->SetBinContent(21,215);
   h->SetBinContent(22,244);
   h->SetBinContent(23,261);
   h->SetBinContent(24,241);
   h->SetBinContent(25,279);
   h->SetBinContent(26,307);
   h->SetBinContent(27,363);
   h->SetBinContent(28,362);
   h->SetBinContent(29,393);
   h->SetBinContent(30,416);
   h->SetBinContent(31,434);
   h->SetBinContent(32,467);
   h->SetBinContent(33,489);
   h->SetBinContent(34,520);
   h->SetBinContent(35,534);
   h->SetBinContent(36,521);
   h->SetBinContent(37,536);
   h->SetBinContent(38,587);
   h->SetBinContent(39,585);
   h->SetBinContent(40,674);
   h->SetBinContent(41,720);
   h->SetBinContent(42,643);
   h->SetBinContent(43,664);
   h->SetBinContent(44,772);
   h->SetBinContent(45,717);
   h->SetBinContent(46,768);
   h->SetBinContent(47,784);
   h->SetBinContent(48,825);
   h->SetBinContent(49,870);
   h->SetBinContent(50,876);
   h->SetBinContent(51,914);
   h->SetBinContent(52,884);
   h->SetBinContent(53,872);
   h->SetBinContent(54,910);
   h->SetBinContent(55,828);
   h->SetBinContent(56,861);
   h->SetBinContent(57,954);
   h->SetBinContent(58,970);
   h->SetBinContent(59,978);
   h->SetBinContent(60,961);
   h->SetBinContent(61,1009);
   h->SetBinContent(62,1092);
   h->SetBinContent(63,1019);
   h->SetBinContent(64,1055);
   h->SetBinContent(65,1083);
   h->SetBinContent(66,1115);
   h->SetBinContent(67,1175);
   h->SetBinContent(68,1148);
   h->SetBinContent(69,1111);
   h->SetBinContent(70,1105);
   h->SetBinContent(71,1119);
   h->SetBinContent(72,1102);
   h->SetBinContent(73,1085);
   h->SetBinContent(74,1124);
   h->SetBinContent(75,1151);
   h->SetBinContent(76,1151);
   h->SetBinContent(77,1231);
   h->SetBinContent(78,1179);
   h->SetBinContent(79,1169);
   h->SetBinContent(80,1166);
   h->SetBinContent(81,1150);
   h->SetBinContent(82,1202);
   h->SetBinContent(83,1165);
   h->SetBinContent(84,1185);
   h->SetBinContent(85,1101);
   h->SetBinContent(86,1103);
   h->SetBinContent(87,1220);
   h->SetBinContent(88,1086);
   h->SetBinContent(89,1114);
   h->SetBinContent(90,1179);
   h->SetBinContent(91,1151);
   h->SetBinContent(92,1123);
   h->SetBinContent(93,1145);
   h->SetBinContent(94,999);
   h->SetBinContent(95,1024);
   h->SetBinContent(96,1011);
   h->SetBinContent(97,1096);
   h->SetBinContent(98,1033);
   h->SetBinContent(99,976);
   h->SetBinContent(100,981);
   h->SetBinContent(101,1038);
   h->SetBinContent(102,1011);
   h->SetBinContent(103,942);
   h->SetBinContent(104,974);
   h->SetBinContent(105,1049);
   h->SetBinContent(106,931);
   h->SetBinContent(107,943);
   h->SetBinContent(108,986);
   h->SetBinContent(109,853);
   h->SetBinContent(110,917);
   h->SetBinContent(111,863);
   h->SetBinContent(112,837);
   h->SetBinContent(113,862);
   h->SetBinContent(114,795);
   h->SetBinContent(115,855);
   h->SetBinContent(116,830);
   h->SetBinContent(117,752);
   h->SetBinContent(118,751);
   h->SetBinContent(119,698);
   h->SetBinContent(120,719);
   h->SetBinContent(121,647);
   h->SetBinContent(122,643);
   h->SetBinContent(123,636);
   h->SetBinContent(124,570);
   h->SetBinContent(125,541);
   h->SetBinContent(126,570);
   h->SetBinContent(127,544);
   h->SetBinContent(128,500);
   h->SetBinContent(129,460);
   h->SetBinContent(130,506);
   h->SetBinContent(131,480);
   h->SetBinContent(132,392);
   h->SetBinContent(133,397);
   h->SetBinContent(134,426);
   h->SetBinContent(135,340);
   h->SetBinContent(136,332);
   h->SetBinContent(137,289);
   h->SetBinContent(138,273);
   h->SetBinContent(139,269);
   h->SetBinContent(140,224);
   h->SetBinContent(141,230);
   h->SetBinContent(142,201);
   h->SetBinContent(143,167);
   h->SetBinContent(144,186);
   h->SetBinContent(145,149);
   h->SetBinContent(146,152);
   h->SetBinContent(147,121);
   h->SetBinContent(148,111);
   h->SetBinContent(149,96);
   h->SetBinContent(150,91);
   h->SetBinContent(151,78);
   h->SetBinContent(152,55);
   h->SetBinContent(153,60);
   h->SetBinContent(154,29);
   h->SetBinContent(155,26);
   h->SetBinContent(156,23);
   h->SetBinContent(157,22);
   h->SetBinContent(158,23);
   h->SetBinContent(159,12);
   h->SetBinContent(160,8);
   h->SetBinContent(161,3);
   h->SetBinContent(165,1);
   h->SetBinContent(166,1);
   h->SetBinContent(167,1);

   TF1 *abeta = new TF1("abeta","(x > 1.022)*[2]*((1+[1]*(x-1.022))*((x-1.022)+0.51099893)*sqrt(((x-1.022)+0.51099893)^2-0.51099893^2)*([0]-(x-1.022))^2*((x-1.022)+0.51099893)/sqrt(((x-1.022)+0.51099893)^2-0.51099893^2)*exp((((x-1.022)<1.2*0.51099893)*(-0.811+4.46e-2*6+1.08e-4*6*6)+((x-1.022)>=1.2*0.51099893)*(-8.46e-2+2.48e-2*6+2.37e-4*6*6))+(((x-1.022)<1.2*0.51099893)*(0.673-1.82e-2*6+6.38e-5*6*6)+((x-1.022)>=1.2*0.51099893)*(1.15e-2+3.58e-4*6-6.17e-5*6*6))*sqrt(((x-1.022)+0.51099893)/0.51099893-1)))",0,20);

   abeta->SetParameter(1, 0); // forbiddenness shape correction
   abeta->SetParameter(2, 1); // normalization
   abeta->SetParameter(0, 14); // Effective q_beta

   const double lowmc = 3.0; // Must be a bin edge!

   h->Fit("abeta", "li", "", lowmc, 14);

   const double eff = (cutlow > lowmc ? 
      h->Integral(h->FindBin(cutlow), h->GetNbinsX()):
      abeta->Integral(cutlow, lowmc) 
    + h->Integral(h->FindBin(lowmc), h->GetNbinsX()))
     /
     (abeta->Integral(1.022, lowmc)
      + h->Integral(h->FindBin(lowmc), h->GetNbinsX()));


   printf("If your cut is not a multiple of %.2f, results are wrong\n",
          h->GetBinWidth(1));
   printf("TECHNOTE: B-8 MC efficiency for %fMeV is %.2f%%\n", cutlow, eff*100);
  printf("const double b8energyeff%.0fMeV = %.16f;\n", cutlow, eff);
}
