#include "TTree.h"
#include "TH1.h"
#include "TF1.h"
#include "TROOT.h"

void li8cutefficiency_finalfit(const double cutlow = 5)
{
   TH1F *h = new TH1F("h","ctEvisID  {ctX[0]**2+ctX[1]**2 < 1650**2 && abs(ctX[2]) < 1700}",150,0,15);
   h->SetBinContent(1,135793);
   h->SetBinContent(2,25238);
   h->SetBinContent(3,8009);
   h->SetBinContent(4,1680);
   h->SetBinContent(5,424);
   h->SetBinContent(6,303);
   h->SetBinContent(7,263);
   h->SetBinContent(8,308);
   h->SetBinContent(9,299);
   h->SetBinContent(10,339);
   h->SetBinContent(11,391);
   h->SetBinContent(12,398);
   h->SetBinContent(13,426);
   h->SetBinContent(14,490);
   h->SetBinContent(15,471);
   h->SetBinContent(16,547);
   h->SetBinContent(17,523);
   h->SetBinContent(18,547);
   h->SetBinContent(19,574);
   h->SetBinContent(20,601);
   h->SetBinContent(21,673);
   h->SetBinContent(22,755);
   h->SetBinContent(23,742);
   h->SetBinContent(24,794);
   h->SetBinContent(25,856);
   h->SetBinContent(26,886);
   h->SetBinContent(27,905);
   h->SetBinContent(28,983);
   h->SetBinContent(29,1022);
   h->SetBinContent(30,1053);
   h->SetBinContent(31,1045);
   h->SetBinContent(32,1183);
   h->SetBinContent(33,1249);
   h->SetBinContent(34,1255);
   h->SetBinContent(35,1325);
   h->SetBinContent(36,1303);
   h->SetBinContent(37,1296);
   h->SetBinContent(38,1344);
   h->SetBinContent(39,1416);
   h->SetBinContent(40,1378);
   h->SetBinContent(41,1430);
   h->SetBinContent(42,1418);
   h->SetBinContent(43,1492);
   h->SetBinContent(44,1588);
   h->SetBinContent(45,1598);
   h->SetBinContent(46,1653);
   h->SetBinContent(47,1631);
   h->SetBinContent(48,1739);
   h->SetBinContent(49,1774);
   h->SetBinContent(50,1777);
   h->SetBinContent(51,1754);
   h->SetBinContent(52,1797);
   h->SetBinContent(53,1789);
   h->SetBinContent(54,1866);
   h->SetBinContent(55,1826);
   h->SetBinContent(56,1802);
   h->SetBinContent(57,1879);
   h->SetBinContent(58,1819);
   h->SetBinContent(59,1848);
   h->SetBinContent(60,1866);
   h->SetBinContent(61,1942);
   h->SetBinContent(62,1959);
   h->SetBinContent(63,1959);
   h->SetBinContent(64,1883);
   h->SetBinContent(65,1911);
   h->SetBinContent(66,1914);
   h->SetBinContent(67,1935);
   h->SetBinContent(68,1947);
   h->SetBinContent(69,1898);
   h->SetBinContent(70,1966);
   h->SetBinContent(71,1930);
   h->SetBinContent(72,1948);
   h->SetBinContent(73,2113);
   h->SetBinContent(74,2010);
   h->SetBinContent(75,1936);
   h->SetBinContent(76,1925);
   h->SetBinContent(77,1965);
   h->SetBinContent(78,1948);
   h->SetBinContent(79,1830);
   h->SetBinContent(80,1910);
   h->SetBinContent(81,1884);
   h->SetBinContent(82,1865);
   h->SetBinContent(83,1812);
   h->SetBinContent(84,1838);
   h->SetBinContent(85,1784);
   h->SetBinContent(86,1735);
   h->SetBinContent(87,1741);
   h->SetBinContent(88,1676);
   h->SetBinContent(89,1640);
   h->SetBinContent(90,1642);
   h->SetBinContent(91,1609);
   h->SetBinContent(92,1607);
   h->SetBinContent(93,1562);
   h->SetBinContent(94,1479);
   h->SetBinContent(95,1513);
   h->SetBinContent(96,1425);
   h->SetBinContent(97,1328);
   h->SetBinContent(98,1389);
   h->SetBinContent(99,1360);
   h->SetBinContent(100,1332);
   h->SetBinContent(101,1247);
   h->SetBinContent(102,1274);
   h->SetBinContent(103,1253);
   h->SetBinContent(104,1217);
   h->SetBinContent(105,1168);
   h->SetBinContent(106,1045);
   h->SetBinContent(107,1068);
   h->SetBinContent(108,1041);
   h->SetBinContent(109,967);
   h->SetBinContent(110,953);
   h->SetBinContent(111,861);
   h->SetBinContent(112,865);
   h->SetBinContent(113,811);
   h->SetBinContent(114,755);
   h->SetBinContent(115,728);
   h->SetBinContent(116,647);
   h->SetBinContent(117,655);
   h->SetBinContent(118,604);
   h->SetBinContent(119,552);
   h->SetBinContent(120,486);
   h->SetBinContent(121,527);
   h->SetBinContent(122,428);
   h->SetBinContent(123,423);
   h->SetBinContent(124,379);
   h->SetBinContent(125,340);
   h->SetBinContent(126,306);
   h->SetBinContent(127,284);
   h->SetBinContent(128,249);
   h->SetBinContent(129,205);
   h->SetBinContent(130,186);
   h->SetBinContent(131,145);
   h->SetBinContent(132,136);
   h->SetBinContent(133,116);
   h->SetBinContent(134,104);
   h->SetBinContent(135,89);
   h->SetBinContent(136,63);
   h->SetBinContent(137,46);
   h->SetBinContent(138,36);
   h->SetBinContent(139,29);
   h->SetBinContent(140,13);
   h->SetBinContent(141,12);
   h->SetBinContent(142,9);
   h->SetBinContent(143,2);
   h->SetBinContent(144,3);
   h->SetBinContent(145,2);
   h->SetBinContent(146,2);

   TF1 *abeta = new TF1("abeta","[2]*((1+[1]*x)*(x+0.51099893)*sqrt((x+0.51099893)^2-0.51099893^2)*([0]-x)^2*(x+0.51099893)/sqrt((x+0.51099893)^2-0.51099893^2)*exp(((x<1.2*0.51099893)*(-0.811+4.46e-2*6+1.08e-4*6*6)+(x>=1.2*0.51099893)*(-8.46e-2+2.48e-2*6+2.37e-4*6*6))+((x<1.2*0.51099893)*(0.673-1.82e-2*6+6.38e-5*6*6)+(x>=1.2*0.51099893)*(1.15e-2+3.58e-4*6-6.17e-5*6*6))*sqrt((x+0.51099893)/0.51099893-1)))",0,15);

   abeta->FixParameter(1, 0); // forbiddenness shape correction
   abeta->SetParameter(2, 1); // normalization
   abeta->SetParameter(0, 14); // Effective q_beta

   h->Fit("abeta", "li", "", 1, 10);

   const double qb = abeta->GetParameter(0);

   const double eff = abeta->Integral(cutlow, qb)/abeta->Integral(0, qb);
   printf("If your cut is not a multiple of %.2f, results are wrong\n",
          h->GetBinWidth(1));
   printf("TECHNOTE: Li-8 MC efficiency for %fMeV is %.2f%%\n", cutlow, eff*100);
   printf("const double li8energyeff%.0fMeV = %f;\n", cutlow, eff);
}
