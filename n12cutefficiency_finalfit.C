#include "TTree.h"
#include "TROOT.h"
#include "consts.h"

void n12cutefficiency_finalfit(const double cutlow = 4)
{
  TH1F *h = new TH1F("h","ctEvisID  {ctX[0]**2+ctX[1]**2 < (1150+250)**2 && abs(ctX[2]) < 1229+250}",200,0,20);
   h->SetBinContent(1,11471);
   h->SetBinContent(2,1199);
   h->SetBinContent(3,487);
   h->SetBinContent(4,146);
   h->SetBinContent(5,73);
   h->SetBinContent(6,78);
   h->SetBinContent(7,43);
   h->SetBinContent(8,33);
   h->SetBinContent(9,19);
   h->SetBinContent(10,18);
   h->SetBinContent(11,19);
   h->SetBinContent(12,10);
   h->SetBinContent(13,7);
   h->SetBinContent(14,15);
   h->SetBinContent(15,17);
   h->SetBinContent(16,16);
   h->SetBinContent(17,22);
   h->SetBinContent(18,25);
   h->SetBinContent(19,25);
   h->SetBinContent(20,24);
   h->SetBinContent(21,22);
   h->SetBinContent(22,26);
   h->SetBinContent(23,29);
   h->SetBinContent(24,31);
   h->SetBinContent(25,37);
   h->SetBinContent(26,50);
   h->SetBinContent(27,46);
   h->SetBinContent(28,41);
   h->SetBinContent(29,41);
   h->SetBinContent(30,33);
   h->SetBinContent(31,43);
   h->SetBinContent(32,50);
   h->SetBinContent(33,45);
   h->SetBinContent(34,48);
   h->SetBinContent(35,53);
   h->SetBinContent(36,68);
   h->SetBinContent(37,54);
   h->SetBinContent(38,61);
   h->SetBinContent(39,53);
   h->SetBinContent(40,72);
   h->SetBinContent(41,55);
   h->SetBinContent(42,70);
   h->SetBinContent(43,72);
   h->SetBinContent(44,74);
   h->SetBinContent(45,79);
   h->SetBinContent(46,98);
   h->SetBinContent(47,106);
   h->SetBinContent(48,87);
   h->SetBinContent(49,105);
   h->SetBinContent(50,89);
   h->SetBinContent(51,99);
   h->SetBinContent(52,96);
   h->SetBinContent(53,80);
   h->SetBinContent(54,82);
   h->SetBinContent(55,93);
   h->SetBinContent(56,93);
   h->SetBinContent(57,88);
   h->SetBinContent(58,89);
   h->SetBinContent(59,95);
   h->SetBinContent(60,87);
   h->SetBinContent(61,104);
   h->SetBinContent(62,113);
   h->SetBinContent(63,109);
   h->SetBinContent(64,89);
   h->SetBinContent(65,114);
   h->SetBinContent(66,99);
   h->SetBinContent(67,106);
   h->SetBinContent(68,90);
   h->SetBinContent(69,104);
   h->SetBinContent(70,101);
   h->SetBinContent(71,97);
   h->SetBinContent(72,99);
   h->SetBinContent(73,121);
   h->SetBinContent(74,117);
   h->SetBinContent(75,109);
   h->SetBinContent(76,116);
   h->SetBinContent(77,94);
   h->SetBinContent(78,117);
   h->SetBinContent(79,121);
   h->SetBinContent(80,121);
   h->SetBinContent(81,120);
   h->SetBinContent(82,116);
   h->SetBinContent(83,112);
   h->SetBinContent(84,117);
   h->SetBinContent(85,124);
   h->SetBinContent(86,123);
   h->SetBinContent(87,126);
   h->SetBinContent(88,125);
   h->SetBinContent(89,114);
   h->SetBinContent(90,123);
   h->SetBinContent(91,100);
   h->SetBinContent(92,113);
   h->SetBinContent(93,95);
   h->SetBinContent(94,93);
   h->SetBinContent(95,110);
   h->SetBinContent(96,110);
   h->SetBinContent(97,115);
   h->SetBinContent(98,117);
   h->SetBinContent(99,108);
   h->SetBinContent(100,121);
   h->SetBinContent(101,106);
   h->SetBinContent(102,83);
   h->SetBinContent(103,76);
   h->SetBinContent(104,98);
   h->SetBinContent(105,114);
   h->SetBinContent(106,90);
   h->SetBinContent(107,104);
   h->SetBinContent(108,107);
   h->SetBinContent(109,115);
   h->SetBinContent(110,99);
   h->SetBinContent(111,83);
   h->SetBinContent(112,91);
   h->SetBinContent(113,76);
   h->SetBinContent(114,89);
   h->SetBinContent(115,96);
   h->SetBinContent(116,106);
   h->SetBinContent(117,94);
   h->SetBinContent(118,73);
   h->SetBinContent(119,99);
   h->SetBinContent(120,93);
   h->SetBinContent(121,87);
   h->SetBinContent(122,79);
   h->SetBinContent(123,105);
   h->SetBinContent(124,81);
   h->SetBinContent(125,109);
   h->SetBinContent(126,86);
   h->SetBinContent(127,53);
   h->SetBinContent(128,78);
   h->SetBinContent(129,86);
   h->SetBinContent(130,79);
   h->SetBinContent(131,57);
   h->SetBinContent(132,69);
   h->SetBinContent(133,80);
   h->SetBinContent(134,67);
   h->SetBinContent(135,61);
   h->SetBinContent(136,67);
   h->SetBinContent(137,67);
   h->SetBinContent(138,61);
   h->SetBinContent(139,70);
   h->SetBinContent(140,71);
   h->SetBinContent(141,59);
   h->SetBinContent(142,48);
   h->SetBinContent(143,54);
   h->SetBinContent(144,45);
   h->SetBinContent(145,39);
   h->SetBinContent(146,41);
   h->SetBinContent(147,46);
   h->SetBinContent(148,44);
   h->SetBinContent(149,46);
   h->SetBinContent(150,46);
   h->SetBinContent(151,50);
   h->SetBinContent(152,41);
   h->SetBinContent(153,30);
   h->SetBinContent(154,33);
   h->SetBinContent(155,30);
   h->SetBinContent(156,24);
   h->SetBinContent(157,30);
   h->SetBinContent(158,21);
   h->SetBinContent(159,34);
   h->SetBinContent(160,27);
   h->SetBinContent(161,30);
   h->SetBinContent(162,25);
   h->SetBinContent(163,20);
   h->SetBinContent(164,28);
   h->SetBinContent(165,22);
   h->SetBinContent(166,14);
   h->SetBinContent(167,17);
   h->SetBinContent(168,10);
   h->SetBinContent(169,18);
   h->SetBinContent(170,15);
   h->SetBinContent(171,11);
   h->SetBinContent(172,12);
   h->SetBinContent(173,10);
   h->SetBinContent(174,5);
   h->SetBinContent(175,5);
   h->SetBinContent(176,3);
   h->SetBinContent(177,3);
   h->SetBinContent(178,3);
   h->SetBinContent(179,3);
   h->SetBinContent(180,4);
   h->SetBinContent(181,1);
   h->SetBinContent(182,3);

   TF1 *abeta = new TF1("abeta","(x > 1.022)*[2]*((1+[1]*(x-1.022))*((x-1.022)+0.51099893)*sqrt(((x-1.022)+0.51099893)^2-0.51099893^2)*([0]-(x-1.022))^2*((x-1.022)+0.51099893)/sqrt(((x-1.022)+0.51099893)^2-0.51099893^2)*exp((((x-1.022)<1.2*0.51099893)*(-0.811+4.46e-2*6+1.08e-4*6*6)+((x-1.022)>=1.2*0.51099893)*(-8.46e-2+2.48e-2*6+2.37e-4*6*6))+(((x-1.022)<1.2*0.51099893)*(0.673-1.82e-2*6+6.38e-5*6*6)+((x-1.022)>=1.2*0.51099893)*(1.15e-2+3.58e-4*6-6.17e-5*6*6))*sqrt(((x-1.022)+0.51099893)/0.51099893-1)))",0,20);

   abeta->SetParameter(1, 0); // forbiddenness shape correction
   abeta->SetParameter(2, 1); // normalization
   abeta->SetParameter(0, 15); // Effective q_beta
   abeta->SetParLimits(0, 13, 19); // Effective q_beta

   const double lowmc = 3.0; // Must be a bin edge!

   h->Fit("abeta", "li", "", 2, 13);

   const double qb = abeta->GetParameter(0);

   const double eff = (cutlow > lowmc ? 
      h->Integral(h->FindBin(cutlow), h->GetNbinsX()):
      abeta->Integral(cutlow, lowmc) 
    + h->Integral(h->FindBin(lowmc), h->GetNbinsX()))
     /
     (abeta->Integral(1.022, lowmc)
      + h->Integral(h->FindBin(lowmc), h->GetNbinsX()));


   printf("If your cut is not a multiple of %.2f, results are wrong\n",
          h->GetBinWidth(1));
   printf("N-12 MC efficiency for %fMeV is %.2f%%\n", cutlow, eff*100);
}
