#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <sstream>
#include <vector>
using std::vector;
#include <iostream>

#include "TTree.h"
#include "TFile.h"
#include "TError.h"
#include "TH1.h"

#include "search.h"

enum searchtype{ b12, be12, neutron, buffer };

// True if we are processing ND data.
static bool near = true;

static double distcut = 0;
static double maxtime = 1000;
static double minenergy = 4;
static double maxenergy = 14;

// Note change, used to cut events on 1000 and define muons at 5000.
// With both at 5000, the event-by-event cut has no additional
// efficiency impact once we evaluate the subsequent muon veto
// efficiency.
const double fido_qiv_muon_def = 5000;

struct cart{
  double x, y, z;

  cart(){}

  cart(double x_, double y_, double z_){
    x=x_, y=y_, z=z_;
  }
};


struct track{
  float x0, y0, z0, x1, y1, z1;
  double tim;
  int nn;
};

static track maketrack(const float x0, const float y0, const float z0,
                       const float x1, const float y1, const float z1,
                       const float tim, const int nn)
{
  track t;
  t.x0 = x0;
  t.y0 = y0;
  t.z0 = z0;
  t.x1 = x1;
  t.y1 = y1;
  t.z1 = z1;
  t.tim = tim;
  t.nn = nn;
  return t;
}

static float ptol(const track & t, const float x, const float y,
                  const float z)
{
  cart entr, exit, point;
  entr.x = t.x0;
  entr.y = t.y0;
  entr.z = t.z0;

  exit.x = t.x1;
  exit.y = t.y1;
  exit.z = t.z1;

  point.x = x;
  point.y = y;
  point.z = z;

  cart eten;
  {
    cart ete;
    ete.x = entr.x - exit.x;
    ete.y = entr.y - exit.y;
    ete.z = entr.z - exit.z;

    eten.x = ete.x/sqrt(ete.x*ete.x+ete.y*ete.y+ete.z*ete.z);
    eten.y = ete.y/sqrt(ete.x*ete.x+ete.y*ete.y+ete.z*ete.z);
    eten.z = ete.z/sqrt(ete.x*ete.x+ete.y*ete.y+ete.z*ete.z);
  }

  cart etp;
  etp.x = entr.x - point.x;
  etp.y = entr.y - point.y;
  etp.z = entr.z - point.z;

  const double detca = fabs(eten.x*etp.x + eten.y*etp.y + eten.z*etp.z);

  cart closest_app;

  closest_app.x = entr.x - eten.x*detca;
  closest_app.y = entr.y - eten.y*detca;
  closest_app.z = entr.z - eten.z*detca;

  return sqrt(pow(closest_app.x - point.x, 2) +
              pow(closest_app.y - point.y, 2) +
              pow(closest_app.z - point.z, 2));
}

static TH1D *numNeutrons_sig = new TH1D("numNeutrons_sig","signal pdf for neutron multiplicity",100,0,100);
static TH1D *muondist_sig = new TH1D("muondist_sig","signal pdf for muon distance to prompt",100,0,5000);
static TH1D *muondist_bkg = new TH1D("muondist_bkg","background pdf for muon distance to prompt",100,0,5000);
static TH1D *numNeutrons_bkg = new TH1D("numNeutrons_bkg","background pdf for neutron multiplicity",100,0,100);

static bool makelikehists()
{
  numNeutrons_sig->SetBinContent(1,0.472182);
  numNeutrons_sig->SetBinContent(2,0.121737);
  numNeutrons_sig->SetBinContent(3,0.0680238);
  numNeutrons_sig->SetBinContent(4,0.042191);
  numNeutrons_sig->SetBinContent(5,0.0321467);
  numNeutrons_sig->SetBinContent(6,0.0265928);
  numNeutrons_sig->SetBinContent(7,0.0219013);
  numNeutrons_sig->SetBinContent(8,0.0178575);
  numNeutrons_sig->SetBinContent(9,0.0155109);
  numNeutrons_sig->SetBinContent(10,0.0136559);
  numNeutrons_sig->SetBinContent(11,0.0117729);
  numNeutrons_sig->SetBinContent(12,0.0100698);
  numNeutrons_sig->SetBinContent(13,0.00962837);
  numNeutrons_sig->SetBinContent(14,0.00815992);
  numNeutrons_sig->SetBinContent(15,0.00779329);
  numNeutrons_sig->SetBinContent(16,0.00677208);
  numNeutrons_sig->SetBinContent(17,0.0067463);
  numNeutrons_sig->SetBinContent(18,0.00669468);
  numNeutrons_sig->SetBinContent(19,0.00511838);
  numNeutrons_sig->SetBinContent(20,0.00501392);
  numNeutrons_sig->SetBinContent(21,0.0047788);
  numNeutrons_sig->SetBinContent(22,0.00506879);
  numNeutrons_sig->SetBinContent(23,0.00459626);
  numNeutrons_sig->SetBinContent(24,0.00393987);
  numNeutrons_sig->SetBinContent(25,0.00362445);
  numNeutrons_sig->SetBinContent(26,0.00299376);
  numNeutrons_sig->SetBinContent(27,0.00328325);
  numNeutrons_sig->SetBinContent(28,0.00288926);
  numNeutrons_sig->SetBinContent(29,0.00199572);
  numNeutrons_sig->SetBinContent(30,0.0025481);
  numNeutrons_sig->SetBinContent(31,0.00254815);
  numNeutrons_sig->SetBinContent(32,0.00223273);
  numNeutrons_sig->SetBinContent(33,0.00186534);
  numNeutrons_sig->SetBinContent(34,0.00170763);
  numNeutrons_sig->SetBinContent(35,0.00181245);
  numNeutrons_sig->SetBinContent(36,0.0015763);
  numNeutrons_sig->SetBinContent(37,0.00183914);
  numNeutrons_sig->SetBinContent(38,0.00136589);
  numNeutrons_sig->SetBinContent(39,0.00131332);
  numNeutrons_sig->SetBinContent(40,0.00128749);
  numNeutrons_sig->SetBinContent(41,0.00118208);
  numNeutrons_sig->SetBinContent(42,0.000945698);
  numNeutrons_sig->SetBinContent(43,0.000604091);
  numNeutrons_sig->SetBinContent(44,0.000604091);
  numNeutrons_sig->SetBinContent(45,0.000604091);
  numNeutrons_sig->SetBinContent(46,0.000604091);
  numNeutrons_sig->SetBinContent(47,0.000604091);
  numNeutrons_sig->SetBinContent(48,0.000604091);
  numNeutrons_sig->SetBinContent(49,0.000604091);
  numNeutrons_sig->SetBinContent(50,0.000604091);
  numNeutrons_sig->SetBinContent(51,0.000604091);
  numNeutrons_sig->SetBinContent(52,0.000604091);
  numNeutrons_sig->SetBinContent(53,0.000604091);
  numNeutrons_sig->SetBinContent(54,0.000604091);
  numNeutrons_sig->SetBinContent(55,0.000604091);
  numNeutrons_sig->SetBinContent(56,0.000604091);
  numNeutrons_sig->SetBinContent(57,0.000604091);
  numNeutrons_sig->SetBinContent(58,0.000604091);
  numNeutrons_sig->SetBinContent(59,0.000604091);
  numNeutrons_sig->SetBinContent(60,0.000604091);
  numNeutrons_sig->SetBinContent(61,0.000604091);
  numNeutrons_sig->SetBinContent(62,0.000604091);
  numNeutrons_sig->SetBinContent(63,0.000604091);
  numNeutrons_sig->SetBinContent(64,0.000604091);
  numNeutrons_sig->SetBinContent(65,0.000604091);
  numNeutrons_sig->SetBinContent(66,0.000604091);
  numNeutrons_sig->SetBinContent(67,0.000604091);
  numNeutrons_sig->SetBinContent(68,0.000604091);
  numNeutrons_sig->SetBinContent(69,0.000604091);
  numNeutrons_sig->SetBinContent(70,0.000604091);
  numNeutrons_sig->SetBinContent(71,0.000604091);
  numNeutrons_sig->SetBinContent(72,0.000604091);
  numNeutrons_sig->SetBinContent(73,0.000604091);
  numNeutrons_sig->SetBinContent(74,0.000604091);
  numNeutrons_sig->SetBinContent(75,0.000604091);
  numNeutrons_sig->SetBinContent(76,0.000604091);
  numNeutrons_sig->SetBinContent(77,0.000604091);
  numNeutrons_sig->SetBinContent(78,0.000604091);
  numNeutrons_sig->SetBinContent(79,0.000604091);
  numNeutrons_sig->SetBinContent(80,0.000604091);
  numNeutrons_sig->SetBinContent(81,0.000604091);
  numNeutrons_sig->SetBinContent(82,0.000604091);
  numNeutrons_sig->SetBinContent(83,0.000604091);
  numNeutrons_sig->SetBinContent(84,0.000604091);
  numNeutrons_sig->SetBinContent(85,0.000604091);
  numNeutrons_sig->SetBinContent(86,0.000604091);
  numNeutrons_sig->SetBinContent(87,0.000604091);
  numNeutrons_sig->SetBinContent(88,0.000604091);
  numNeutrons_sig->SetBinContent(89,0.000604091);
  numNeutrons_sig->SetBinContent(90,0.000604091);
  numNeutrons_sig->SetBinContent(91,0.000604091);
  numNeutrons_sig->SetBinContent(92,0.000604091);
  numNeutrons_sig->SetBinContent(93,0.000604091);
  numNeutrons_sig->SetBinContent(94,0.000604091);
  numNeutrons_sig->SetBinContent(95,0.000604091);
  numNeutrons_sig->SetBinContent(96,0.000604091);
  numNeutrons_sig->SetBinContent(97,0.000604091);
  numNeutrons_sig->SetBinContent(98,0.000604091);
  numNeutrons_sig->SetBinContent(99,0.000604091);
  numNeutrons_sig->SetBinContent(100,0.000604091);
  numNeutrons_sig->SetBinError(1,0.00807567);
  numNeutrons_sig->SetBinError(2,0.00188551);
  numNeutrons_sig->SetBinError(3,0.00135412);
  numNeutrons_sig->SetBinError(4,0.00106088);
  numNeutrons_sig->SetBinError(5,0.000923336);
  numNeutrons_sig->SetBinError(6,0.000838639);
  numNeutrons_sig->SetBinError(7,0.000760437);
  numNeutrons_sig->SetBinError(8,0.000686424);
  numNeutrons_sig->SetBinError(9,0.00063953);
  numNeutrons_sig->SetBinError(10,0.000599955);
  numNeutrons_sig->SetBinError(11,0.000556959);
  numNeutrons_sig->SetBinError(12,0.000515069);
  numNeutrons_sig->SetBinError(13,0.000503539);
  numNeutrons_sig->SetBinError(14,0.000463532);
  numNeutrons_sig->SetBinError(15,0.000452979);
  numNeutrons_sig->SetBinError(16,0.000422191);
  numNeutrons_sig->SetBinError(17,0.000421372);
  numNeutrons_sig->SetBinError(18,0.000419729);
  numNeutrons_sig->SetBinError(19,0.000367043);
  numNeutrons_sig->SetBinError(20,0.000363259);
  numNeutrons_sig->SetBinError(21,0.000354597);
  numNeutrons_sig->SetBinError(22,0.000365156);
  numNeutrons_sig->SetBinError(23,0.000347711);
  numNeutrons_sig->SetBinError(24,0.000321918);
  numNeutrons_sig->SetBinError(25,0.000308773);
  numNeutrons_sig->SetBinError(26,0.000280641);
  numNeutrons_sig->SetBinError(27,0.000293869);
  numNeutrons_sig->SetBinError(28,0.000275674);
  numNeutrons_sig->SetBinError(29,0.000229143);
  numNeutrons_sig->SetBinError(30,0.000258872);
  numNeutrons_sig->SetBinError(31,0.000258872);
  numNeutrons_sig->SetBinError(32,0.000242331);
  numNeutrons_sig->SetBinError(33,0.000221477);
  numNeutrons_sig->SetBinError(34,0.000211912);
  numNeutrons_sig->SetBinError(35,0.000218335);
  numNeutrons_sig->SetBinError(36,0.000203599);
  numNeutrons_sig->SetBinError(37,0.000219912);
  numNeutrons_sig->SetBinError(38,0.00018954);
  numNeutrons_sig->SetBinError(39,0.000185859);
  numNeutrons_sig->SetBinError(40,0.000183991);
  numNeutrons_sig->SetBinError(41,0.000176322);
  numNeutrons_sig->SetBinError(42,0.000157707);
  numNeutrons_sig->SetBinError(43,0.000126056);
  numNeutrons_sig->SetEntries(12013.9);
  numNeutrons_sig->GetXaxis()->SetTitle("number of neutrons ");
  numNeutrons_sig->GetYaxis()->SetTitle("events");

  muondist_sig->SetBinContent(1,0.0405368);
  muondist_sig->SetBinContent(2,0.0996113);
  muondist_sig->SetBinContent(3,0.113459);
  muondist_sig->SetBinContent(4,0.098744);
  muondist_sig->SetBinContent(5,0.0761948);
  muondist_sig->SetBinContent(6,0.0559451);
  muondist_sig->SetBinContent(7,0.0436028);
  muondist_sig->SetBinContent(8,0.0377785);
  muondist_sig->SetBinContent(9,0.0315165);
  muondist_sig->SetBinContent(10,0.0292702);
  muondist_sig->SetBinContent(11,0.0260321);
  muondist_sig->SetBinContent(12,0.0260774);
  muondist_sig->SetBinContent(13,0.0242399);
  muondist_sig->SetBinContent(14,0.0229207);
  muondist_sig->SetBinContent(15,0.018592);
  muondist_sig->SetBinContent(16,0.0196051);
  muondist_sig->SetBinContent(17,0.0205155);
  muondist_sig->SetBinContent(18,0.0165831);
  muondist_sig->SetBinContent(19,0.0160644);
  muondist_sig->SetBinContent(20,0.0173723);
  muondist_sig->SetBinContent(21,0.0128423);
  muondist_sig->SetBinContent(22,0.0113273);
  muondist_sig->SetBinContent(23,0.0101852);
  muondist_sig->SetBinContent(24,0.0106873);
  muondist_sig->SetBinContent(25,0.011475);
  muondist_sig->SetBinContent(26,0.00981208);
  muondist_sig->SetBinContent(27,0.00985943);
  muondist_sig->SetBinContent(28,0.00464216);
  muondist_sig->SetBinContent(29,0.00856952);
  muondist_sig->SetBinContent(30,0.00564496);
  muondist_sig->SetBinContent(31,0.00777393);
  muondist_sig->SetBinContent(32,0.00572662);
  muondist_sig->SetBinContent(33,0.00332412);
  muondist_sig->SetBinContent(34,0.00488976);
  muondist_sig->SetBinContent(35,0.00278151);
  muondist_sig->SetBinContent(36,0.00507715);
  muondist_sig->SetBinContent(37,0.00513137);
  muondist_sig->SetBinContent(38,0.00370988);
  muondist_sig->SetBinContent(39,0.00176201);
  muondist_sig->SetBinContent(40,0.00262061);
  muondist_sig->SetBinContent(41,0.00277191);
  muondist_sig->SetBinContent(43,0.00265106);
  muondist_sig->SetBinContent(44,0.0016982);
  muondist_sig->SetBinContent(45,0.00204477);
  muondist_sig->SetBinContent(46,8.69766e-05);
  muondist_sig->SetBinContent(47,0.000252734);
  muondist_sig->SetBinContent(48,0.00338045);
  muondist_sig->SetBinContent(49,0.000706991);
  muondist_sig->SetBinContent(50,0.000313267);
  muondist_sig->SetBinContent(51,0.00108307);
  muondist_sig->SetBinContent(54,0.00106081);
  muondist_sig->SetBinContent(55,0.000671631);
  muondist_sig->SetBinContent(56,0.000875449);
  muondist_sig->SetBinContent(57,0.000830845);
  muondist_sig->SetBinContent(58,0.00181829);
  muondist_sig->SetBinContent(59,0.000521821);
  muondist_sig->SetBinContent(60,0.000844004);
  muondist_sig->SetBinContent(62,0.000684597);
  muondist_sig->SetBinContent(65,0.00040473);
  muondist_sig->SetBinContent(66,8.61094e-06);
  muondist_sig->SetBinContent(70,0.000170859);
  muondist_sig->SetBinContent(71,0.000153828);
  muondist_sig->SetBinContent(72,0.000153828);
  muondist_sig->SetBinContent(73,0.000153828);
  muondist_sig->SetBinContent(74,0.000153828);
  muondist_sig->SetBinContent(75,0.000153828);
  muondist_sig->SetBinContent(76,0.000153828);
  muondist_sig->SetBinContent(77,0.000153828);
  muondist_sig->SetBinContent(78,0.000153828);
  muondist_sig->SetBinContent(79,0.000153828);
  muondist_sig->SetBinContent(80,0.000153828);
  muondist_sig->SetBinContent(81,0.000153828);
  muondist_sig->SetBinContent(82,0.000153828);
  muondist_sig->SetBinContent(83,0.000153828);
  muondist_sig->SetBinContent(84,0.000153828);
  muondist_sig->SetBinContent(85,0.000153828);
  muondist_sig->SetBinContent(86,0.000153828);
  muondist_sig->SetBinContent(87,0.000153828);
  muondist_sig->SetBinContent(88,0.000153828);
  muondist_sig->SetBinContent(89,0.000153828);
  muondist_sig->SetBinContent(90,0.000153828);
  muondist_sig->SetBinContent(91,0.000153828);
  muondist_sig->SetBinContent(92,0.000153828);
  muondist_sig->SetBinContent(93,0.000153828);
  muondist_sig->SetBinContent(94,0.000153828);
  muondist_sig->SetBinContent(95,0.000153828);
  muondist_sig->SetBinContent(96,0.000153828);
  muondist_sig->SetBinContent(97,0.000153828);
  muondist_sig->SetBinContent(98,0.000153828);
  muondist_sig->SetBinContent(99,0.000153828);
  muondist_sig->SetBinContent(100,0.000153828);
  muondist_sig->SetBinContent(101,-9.70326e-07);
  muondist_sig->SetBinError(1,0.00105894);
  muondist_sig->SetBinError(2,0.00166514);
  muondist_sig->SetBinError(3,0.00179109);
  muondist_sig->SetBinError(4,0.00169553);
  muondist_sig->SetBinError(5,0.00152662);
  muondist_sig->SetBinError(6,0.00136011);
  muondist_sig->SetBinError(7,0.00125683);
  muondist_sig->SetBinError(8,0.00121802);
  muondist_sig->SetBinError(9,0.00117332);
  muondist_sig->SetBinError(10,0.00117368);
  muondist_sig->SetBinError(11,0.00116351);
  muondist_sig->SetBinError(12,0.00119028);
  muondist_sig->SetBinError(13,0.00119487);
  muondist_sig->SetBinError(14,0.00120666);
  muondist_sig->SetBinError(15,0.00118218);
  muondist_sig->SetBinError(16,0.00121837);
  muondist_sig->SetBinError(17,0.00125032);
  muondist_sig->SetBinError(18,0.00122878);
  muondist_sig->SetBinError(19,0.00124344);
  muondist_sig->SetBinError(20,0.00127476);
  muondist_sig->SetBinError(21,0.00124292);
  muondist_sig->SetBinError(22,0.00124323);
  muondist_sig->SetBinError(23,0.00124326);
  muondist_sig->SetBinError(24,0.0012606);
  muondist_sig->SetBinError(25,0.00127969);
  muondist_sig->SetBinError(26,0.00127091);
  muondist_sig->SetBinError(27,0.0012783);
  muondist_sig->SetBinError(28,0.00122694);
  muondist_sig->SetBinError(29,0.00127321);
  muondist_sig->SetBinError(30,0.00124421);
  muondist_sig->SetBinError(31,0.00126781);
  muondist_sig->SetBinError(32,0.00124654);
  muondist_sig->SetBinError(33,0.00121838);
  muondist_sig->SetBinError(34,0.00123165);
  muondist_sig->SetBinError(35,0.00120223);
  muondist_sig->SetBinError(36,0.00122013);
  muondist_sig->SetBinError(37,0.00121149);
  muondist_sig->SetBinError(38,0.00118553);
  muondist_sig->SetBinError(39,0.00115149);
  muondist_sig->SetBinError(40,0.001148);
  muondist_sig->SetBinError(41,0.00113563);
  muondist_sig->SetBinError(43,0.00110066);
  muondist_sig->SetBinError(44,0.00107066);
  muondist_sig->SetBinError(45,0.00105395);
  muondist_sig->SetBinError(46,0.00100619);
  muondist_sig->SetBinError(47,0.000984358);
  muondist_sig->SetBinError(48,0.00100179);
  muondist_sig->SetBinError(49,0.000937238);
  muondist_sig->SetBinError(50,0.000902209);
  muondist_sig->SetBinError(51,0.000883124);
  muondist_sig->SetBinError(54,0.000778525);
  muondist_sig->SetBinError(55,0.000734502);
  muondist_sig->SetBinError(56,0.000700171);
  muondist_sig->SetBinError(57,0.000660247);
  muondist_sig->SetBinError(58,0.000645235);
  muondist_sig->SetBinError(59,0.000581528);
  muondist_sig->SetBinError(60,0.000554747);
  muondist_sig->SetBinError(62,0.000483615);
  muondist_sig->SetBinError(65,0.0003796);
  muondist_sig->SetBinError(66,0.000332868);
  muondist_sig->SetBinError(70,0.000225135);
  muondist_sig->SetBinError(71,0.000197736);
  muondist_sig->SetBinError(101,2.11742e-07);
  muondist_sig->SetEntries(12178.7);
  muondist_sig->GetXaxis()->SetTitle("distance from mu [mm]");
  muondist_sig->GetYaxis()->SetTitle("events");

  numNeutrons_bkg->SetBinContent(1,0.991708);
  numNeutrons_bkg->SetBinContent(2,0.00667387);
  numNeutrons_bkg->SetBinContent(3,0.000858013);
  numNeutrons_bkg->SetBinContent(4,0.000309747);
  numNeutrons_bkg->SetBinContent(5,0.000142589);
  numNeutrons_bkg->SetBinContent(6,8.15017e-05);
  numNeutrons_bkg->SetBinContent(7,4.88698e-05);
  numNeutrons_bkg->SetBinContent(8,3.38604e-05);
  numNeutrons_bkg->SetBinContent(9,2.44572e-05);
  numNeutrons_bkg->SetBinContent(10,1.89181e-05);
  numNeutrons_bkg->SetBinContent(11,1.425e-05);
  numNeutrons_bkg->SetBinContent(12,1.15697e-05);
  numNeutrons_bkg->SetBinContent(13,8.91181e-06);
  numNeutrons_bkg->SetBinContent(14,7.19199e-06);
  numNeutrons_bkg->SetBinContent(15,6.52193e-06);
  numNeutrons_bkg->SetBinContent(16,4.60109e-06);
  numNeutrons_bkg->SetBinContent(17,4.3554e-06);
  numNeutrons_bkg->SetBinContent(18,3.88635e-06);
  numNeutrons_bkg->SetBinContent(19,3.50665e-06);
  numNeutrons_bkg->SetBinContent(20,3.17162e-06);
  numNeutrons_bkg->SetBinContent(21,2.45689e-06);
  numNeutrons_bkg->SetBinContent(22,2.03252e-06);
  numNeutrons_bkg->SetBinContent(23,1.74216e-06);
  numNeutrons_bkg->SetBinContent(24,1.38479e-06);
  numNeutrons_bkg->SetBinContent(25,1.38479e-06);
  numNeutrons_bkg->SetBinContent(26,1.31779e-06);
  numNeutrons_bkg->SetBinContent(27,1.1391e-06);
  numNeutrons_bkg->SetBinContent(28,1.00509e-06);
  numNeutrons_bkg->SetBinContent(29,9.38085e-07);
  numNeutrons_bkg->SetBinContent(30,7.37067e-07);
  numNeutrons_bkg->SetBinContent(31,7.14732e-07);
  numNeutrons_bkg->SetBinContent(32,7.14732e-07);
  numNeutrons_bkg->SetBinContent(33,4.24372e-07);
  numNeutrons_bkg->SetBinContent(34,4.24372e-07);
  numNeutrons_bkg->SetBinContent(35,5.8072e-07);
  numNeutrons_bkg->SetBinContent(36,3.79701e-07);
  numNeutrons_bkg->SetBinContent(37,3.79701e-07);
  numNeutrons_bkg->SetBinContent(38,4.46707e-07);
  numNeutrons_bkg->SetBinContent(39,4.46707e-07);
  numNeutrons_bkg->SetBinContent(40,2.23354e-07);
  numNeutrons_bkg->SetBinContent(41,3.57366e-07);
  numNeutrons_bkg->SetBinContent(42,2.68024e-07);
  numNeutrons_bkg->SetBinContent(43,2.23354e-07);
  numNeutrons_bkg->SetBinContent(44,2.23354e-07);
  numNeutrons_bkg->SetBinContent(45,2.23354e-07);
  numNeutrons_bkg->SetBinContent(46,2.23354e-07);
  numNeutrons_bkg->SetBinContent(47,2.23354e-07);
  numNeutrons_bkg->SetBinContent(48,2.23354e-07);
  numNeutrons_bkg->SetBinContent(49,2.23354e-07);
  numNeutrons_bkg->SetBinContent(50,2.23354e-07);
  numNeutrons_bkg->SetBinContent(51,2.23354e-07);
  numNeutrons_bkg->SetBinContent(52,2.23354e-07);
  numNeutrons_bkg->SetBinContent(53,2.23354e-07);
  numNeutrons_bkg->SetBinContent(54,2.23354e-07);
  numNeutrons_bkg->SetBinContent(55,2.23354e-07);
  numNeutrons_bkg->SetBinContent(56,2.23354e-07);
  numNeutrons_bkg->SetBinContent(57,2.23354e-07);
  numNeutrons_bkg->SetBinContent(58,2.23354e-07);
  numNeutrons_bkg->SetBinContent(59,2.23354e-07);
  numNeutrons_bkg->SetBinContent(60,2.23354e-07);
  numNeutrons_bkg->SetBinContent(61,2.23354e-07);
  numNeutrons_bkg->SetBinContent(62,2.23354e-07);
  numNeutrons_bkg->SetBinContent(63,2.23354e-07);
  numNeutrons_bkg->SetBinContent(64,2.23354e-07);
  numNeutrons_bkg->SetBinContent(65,2.23354e-07);
  numNeutrons_bkg->SetBinContent(66,2.23354e-07);
  numNeutrons_bkg->SetBinContent(67,2.23354e-07);
  numNeutrons_bkg->SetBinContent(68,2.23354e-07);
  numNeutrons_bkg->SetBinContent(69,2.23354e-07);
  numNeutrons_bkg->SetBinContent(70,2.23354e-07);
  numNeutrons_bkg->SetBinContent(71,2.23354e-07);
  numNeutrons_bkg->SetBinContent(72,2.23354e-07);
  numNeutrons_bkg->SetBinContent(73,2.23354e-07);
  numNeutrons_bkg->SetBinContent(74,2.23354e-07);
  numNeutrons_bkg->SetBinContent(75,2.23354e-07);
  numNeutrons_bkg->SetBinContent(76,2.23354e-07);
  numNeutrons_bkg->SetBinContent(77,2.23354e-07);
  numNeutrons_bkg->SetBinContent(78,2.23354e-07);
  numNeutrons_bkg->SetBinContent(79,2.23354e-07);
  numNeutrons_bkg->SetBinContent(80,2.23354e-07);
  numNeutrons_bkg->SetBinContent(81,2.23354e-07);
  numNeutrons_bkg->SetBinContent(82,2.23354e-07);
  numNeutrons_bkg->SetBinContent(83,2.23354e-07);
  numNeutrons_bkg->SetBinContent(84,2.23354e-07);
  numNeutrons_bkg->SetBinContent(85,2.23354e-07);
  numNeutrons_bkg->SetBinContent(86,2.23354e-07);
  numNeutrons_bkg->SetBinContent(87,2.23354e-07);
  numNeutrons_bkg->SetBinContent(88,2.23354e-07);
  numNeutrons_bkg->SetBinContent(89,2.23354e-07);
  numNeutrons_bkg->SetBinContent(90,2.23354e-07);
  numNeutrons_bkg->SetBinContent(91,2.23354e-07);
  numNeutrons_bkg->SetBinContent(92,2.23354e-07);
  numNeutrons_bkg->SetBinContent(93,2.23354e-07);
  numNeutrons_bkg->SetBinContent(94,2.23354e-07);
  numNeutrons_bkg->SetBinContent(95,2.23354e-07);
  numNeutrons_bkg->SetBinContent(96,2.23354e-07);
  numNeutrons_bkg->SetBinContent(97,2.23354e-07);
  numNeutrons_bkg->SetBinContent(98,2.23354e-07);
  numNeutrons_bkg->SetBinContent(99,2.23354e-07);
  numNeutrons_bkg->SetBinContent(100,2.23354e-07);
  numNeutrons_bkg->SetBinError(1,0.000148829);
  numNeutrons_bkg->SetBinError(2,1.22092e-05);
  numNeutrons_bkg->SetBinError(3,4.37767e-06);
  numNeutrons_bkg->SetBinError(4,2.63027e-06);
  numNeutrons_bkg->SetBinError(5,1.78459e-06);
  numNeutrons_bkg->SetBinError(6,1.34921e-06);
  numNeutrons_bkg->SetBinError(7,1.04476e-06);
  numNeutrons_bkg->SetBinError(8,8.69646e-07);
  numNeutrons_bkg->SetBinError(9,7.39095e-07);
  numNeutrons_bkg->SetBinError(10,6.50032e-07);
  numNeutrons_bkg->SetBinError(11,5.64161e-07);
  numNeutrons_bkg->SetBinError(12,5.08344e-07);
  numNeutrons_bkg->SetBinError(13,4.46149e-07);
  numNeutrons_bkg->SetBinError(14,4.00794e-07);
  numNeutrons_bkg->SetBinError(15,3.81667e-07);
  numNeutrons_bkg->SetBinError(16,3.20573e-07);
  numNeutrons_bkg->SetBinError(17,3.11896e-07);
  numNeutrons_bkg->SetBinError(18,2.94624e-07);
  numNeutrons_bkg->SetBinError(19,2.79861e-07);
  numNeutrons_bkg->SetBinError(20,2.66157e-07);
  numNeutrons_bkg->SetBinError(21,2.34255e-07);
  numNeutrons_bkg->SetBinError(22,2.13066e-07);
  numNeutrons_bkg->SetBinError(23,1.97261e-07);
  numNeutrons_bkg->SetBinError(24,1.75869e-07);
  numNeutrons_bkg->SetBinError(25,1.75869e-07);
  numNeutrons_bkg->SetBinError(26,1.71561e-07);
  numNeutrons_bkg->SetBinError(27,1.59506e-07);
  numNeutrons_bkg->SetBinError(28,1.4983e-07);
  numNeutrons_bkg->SetBinError(29,1.4475e-07);
  numNeutrons_bkg->SetBinError(30,1.28307e-07);
  numNeutrons_bkg->SetBinError(31,1.26348e-07);
  numNeutrons_bkg->SetBinError(32,1.26348e-07);
  numNeutrons_bkg->SetBinError(33,9.73576e-08);
  numNeutrons_bkg->SetBinError(34,9.73576e-08);
  numNeutrons_bkg->SetBinError(35,1.13888e-07);
  numNeutrons_bkg->SetBinError(36,9.20911e-08);
  numNeutrons_bkg->SetBinError(37,9.20911e-08);
  numNeutrons_bkg->SetBinError(38,9.98868e-08);
  numNeutrons_bkg->SetBinError(39,9.98868e-08);
  numNeutrons_bkg->SetBinError(40,7.06306e-08);
  numNeutrons_bkg->SetBinError(41,8.93415e-08);
  numNeutrons_bkg->SetBinError(42,7.7372e-08);
  numNeutrons_bkg->SetBinError(43,7.06306e-08);
  numNeutrons_bkg->SetEntries(4.47716e+07);
  numNeutrons_bkg->GetXaxis()->SetTitle("number of neutrons ");
  numNeutrons_bkg->GetYaxis()->SetTitle("events");

  muondist_bkg->SetBinContent(1,0.000554406);
  muondist_bkg->SetBinContent(2,0.00166561);
  muondist_bkg->SetBinContent(3,0.00277728);
  muondist_bkg->SetBinContent(4,0.0038804);
  muondist_bkg->SetBinContent(5,0.0049752);
  muondist_bkg->SetBinContent(6,0.00610167);
  muondist_bkg->SetBinContent(7,0.00718865);
  muondist_bkg->SetBinContent(8,0.00825752);
  muondist_bkg->SetBinContent(9,0.00934276);
  muondist_bkg->SetBinContent(10,0.0104235);
  muondist_bkg->SetBinContent(11,0.0115412);
  muondist_bkg->SetBinContent(12,0.0126318);
  muondist_bkg->SetBinContent(13,0.0136971);
  muondist_bkg->SetBinContent(14,0.0148229);
  muondist_bkg->SetBinContent(15,0.0158454);
  muondist_bkg->SetBinContent(16,0.0168978);
  muondist_bkg->SetBinContent(17,0.0178583);
  muondist_bkg->SetBinContent(18,0.018782);
  muondist_bkg->SetBinContent(19,0.0196678);
  muondist_bkg->SetBinContent(20,0.0204393);
  muondist_bkg->SetBinContent(21,0.0211741);
  muondist_bkg->SetBinContent(22,0.0219072);
  muondist_bkg->SetBinContent(23,0.0224502);
  muondist_bkg->SetBinContent(24,0.0229785);
  muondist_bkg->SetBinContent(25,0.0234604);
  muondist_bkg->SetBinContent(26,0.0238548);
  muondist_bkg->SetBinContent(27,0.0241646);
  muondist_bkg->SetBinContent(28,0.0243697);
  muondist_bkg->SetBinContent(29,0.0245478);
  muondist_bkg->SetBinContent(30,0.0246471);
  muondist_bkg->SetBinContent(31,0.0246832);
  muondist_bkg->SetBinContent(32,0.0247106);
  muondist_bkg->SetBinContent(33,0.0246256);
  muondist_bkg->SetBinContent(34,0.0244565);
  muondist_bkg->SetBinContent(35,0.0241933);
  muondist_bkg->SetBinContent(36,0.0238688);
  muondist_bkg->SetBinContent(37,0.0234723);
  muondist_bkg->SetBinContent(38,0.0230488);
  muondist_bkg->SetBinContent(39,0.0225687);
  muondist_bkg->SetBinContent(40,0.0220199);
  muondist_bkg->SetBinContent(41,0.0214494);
  muondist_bkg->SetBinContent(42,0.0207921);
  muondist_bkg->SetBinContent(43,0.0201262);
  muondist_bkg->SetBinContent(44,0.0194287);
  muondist_bkg->SetBinContent(45,0.0186376);
  muondist_bkg->SetBinContent(46,0.0178299);
  muondist_bkg->SetBinContent(47,0.016984);
  muondist_bkg->SetBinContent(48,0.0161106);
  muondist_bkg->SetBinContent(49,0.0151701);
  muondist_bkg->SetBinContent(50,0.0142196);
  muondist_bkg->SetBinContent(51,0.0132528);
  muondist_bkg->SetBinContent(52,0.0122223);
  muondist_bkg->SetBinContent(53,0.0112217);
  muondist_bkg->SetBinContent(54,0.0101954);
  muondist_bkg->SetBinContent(55,0.00920432);
  muondist_bkg->SetBinContent(56,0.00823815);
  muondist_bkg->SetBinContent(57,0.0073006);
  muondist_bkg->SetBinContent(58,0.00648598);
  muondist_bkg->SetBinContent(59,0.00572176);
  muondist_bkg->SetBinContent(60,0.00503168);
  muondist_bkg->SetBinContent(61,0.00439957);
  muondist_bkg->SetBinContent(62,0.00380356);
  muondist_bkg->SetBinContent(63,0.00326258);
  muondist_bkg->SetBinContent(64,0.00278945);
  muondist_bkg->SetBinContent(65,0.00235148);
  muondist_bkg->SetBinContent(66,0.00195177);
  muondist_bkg->SetBinContent(67,0.00160316);
  muondist_bkg->SetBinContent(68,0.00129723);
  muondist_bkg->SetBinContent(69,0.00103706);
  muondist_bkg->SetBinContent(70,0.000813606);
  muondist_bkg->SetBinContent(71,0.000617169);
  muondist_bkg->SetBinContent(72,0.000617169);
  muondist_bkg->SetBinContent(73,0.000617169);
  muondist_bkg->SetBinContent(74,0.000617169);
  muondist_bkg->SetBinContent(75,0.000617169);
  muondist_bkg->SetBinContent(76,0.000617169);
  muondist_bkg->SetBinContent(77,0.000617169);
  muondist_bkg->SetBinContent(78,0.000617169);
  muondist_bkg->SetBinContent(79,0.000617169);
  muondist_bkg->SetBinContent(80,0.000617169);
  muondist_bkg->SetBinContent(81,0.000617169);
  muondist_bkg->SetBinContent(82,0.000617169);
  muondist_bkg->SetBinContent(83,0.000617169);
  muondist_bkg->SetBinContent(84,0.000617169);
  muondist_bkg->SetBinContent(85,0.000617169);
  muondist_bkg->SetBinContent(86,0.000617169);
  muondist_bkg->SetBinContent(87,0.000617169);
  muondist_bkg->SetBinContent(88,0.000617169);
  muondist_bkg->SetBinContent(89,0.000617169);
  muondist_bkg->SetBinContent(90,0.000617169);
  muondist_bkg->SetBinContent(91,0.000617169);
  muondist_bkg->SetBinContent(92,0.000617169);
  muondist_bkg->SetBinContent(93,0.000617169);
  muondist_bkg->SetBinContent(94,0.000617169);
  muondist_bkg->SetBinContent(95,0.000617169);
  muondist_bkg->SetBinContent(96,0.000617169);
  muondist_bkg->SetBinContent(97,0.000617169);
  muondist_bkg->SetBinContent(98,0.000617169);
  muondist_bkg->SetBinContent(99,0.000617169);
  muondist_bkg->SetBinContent(100,0.000617169);
  muondist_bkg->SetBinContent(101,4.61328e-07);
  muondist_bkg->SetBinError(1,3.48987e-06);
  muondist_bkg->SetBinError(2,6.04898e-06);
  muondist_bkg->SetBinError(3,7.81097e-06);
  muondist_bkg->SetBinError(4,9.2328e-06);
  muondist_bkg->SetBinError(5,1.04544e-05);
  muondist_bkg->SetBinError(6,1.15776e-05);
  muondist_bkg->SetBinError(7,1.25666e-05);
  muondist_bkg->SetBinError(8,1.34685e-05);
  muondist_bkg->SetBinError(9,1.43263e-05);
  muondist_bkg->SetBinError(10,1.51322e-05);
  muondist_bkg->SetBinError(11,1.59228e-05);
  muondist_bkg->SetBinError(12,1.66582e-05);
  muondist_bkg->SetBinError(13,1.73464e-05);
  muondist_bkg->SetBinError(14,1.80452e-05);
  muondist_bkg->SetBinError(15,1.86572e-05);
  muondist_bkg->SetBinError(16,1.92668e-05);
  muondist_bkg->SetBinError(17,1.98069e-05);
  muondist_bkg->SetBinError(18,2.03126e-05);
  muondist_bkg->SetBinError(19,2.07861e-05);
  muondist_bkg->SetBinError(20,2.11899e-05);
  muondist_bkg->SetBinError(21,2.15674e-05);
  muondist_bkg->SetBinError(22,2.19375e-05);
  muondist_bkg->SetBinError(23,2.22078e-05);
  muondist_bkg->SetBinError(24,2.24675e-05);
  muondist_bkg->SetBinError(25,2.2702e-05);
  muondist_bkg->SetBinError(26,2.28919e-05);
  muondist_bkg->SetBinError(27,2.30401e-05);
  muondist_bkg->SetBinError(28,2.31377e-05);
  muondist_bkg->SetBinError(29,2.32221e-05);
  muondist_bkg->SetBinError(30,2.3269e-05);
  muondist_bkg->SetBinError(31,2.3286e-05);
  muondist_bkg->SetBinError(32,2.3299e-05);
  muondist_bkg->SetBinError(33,2.32589e-05);
  muondist_bkg->SetBinError(34,2.31789e-05);
  muondist_bkg->SetBinError(35,2.30538e-05);
  muondist_bkg->SetBinError(36,2.28987e-05);
  muondist_bkg->SetBinError(37,2.27077e-05);
  muondist_bkg->SetBinError(38,2.25019e-05);
  muondist_bkg->SetBinError(39,2.22663e-05);
  muondist_bkg->SetBinError(40,2.19939e-05);
  muondist_bkg->SetBinError(41,2.17071e-05);
  muondist_bkg->SetBinError(42,2.13719e-05);
  muondist_bkg->SetBinError(43,2.1027e-05);
  muondist_bkg->SetBinError(44,2.06594e-05);
  muondist_bkg->SetBinError(45,2.02344e-05);
  muondist_bkg->SetBinError(46,1.97911e-05);
  muondist_bkg->SetBinError(47,1.93159e-05);
  muondist_bkg->SetBinError(48,1.88127e-05);
  muondist_bkg->SetBinError(49,1.82553e-05);
  muondist_bkg->SetBinError(50,1.76742e-05);
  muondist_bkg->SetBinError(51,1.70627e-05);
  muondist_bkg->SetBinError(52,1.63859e-05);
  muondist_bkg->SetBinError(53,1.57009e-05);
  muondist_bkg->SetBinError(54,1.49657e-05);
  muondist_bkg->SetBinError(55,1.42197e-05);
  muondist_bkg->SetBinError(56,1.34527e-05);
  muondist_bkg->SetBinError(57,1.26641e-05);
  muondist_bkg->SetBinError(58,1.19367e-05);
  muondist_bkg->SetBinError(59,1.12114e-05);
  muondist_bkg->SetBinError(60,1.05136e-05);
  muondist_bkg->SetBinError(61,9.83106e-06);
  muondist_bkg->SetBinError(62,9.14093e-06);
  muondist_bkg->SetBinError(63,8.46594e-06);
  muondist_bkg->SetBinError(64,7.82807e-06);
  muondist_bkg->SetBinError(65,7.18729e-06);
  muondist_bkg->SetBinError(66,6.54801e-06);
  muondist_bkg->SetBinError(67,5.93449e-06);
  muondist_bkg->SetBinError(68,5.33831e-06);
  muondist_bkg->SetBinError(69,4.77307e-06);
  muondist_bkg->SetBinError(70,4.22768e-06);
  muondist_bkg->SetBinError(71,3.68211e-06);
  muondist_bkg->SetBinError(101,1.0067e-07);
  muondist_bkg->SetEntries(4.47715e+07);
  muondist_bkg->GetXaxis()->SetTitle("distance from mu [mm]");
  muondist_bkg->GetYaxis()->SetTitle("events");
  return true;
}

bool triggermakelikehists = makelikehists();

static float pdf(float x, TH1D * const hist, const bool interpolate)
{
  const float histmax=hist->GetXaxis()->GetXmax();
  const float histmin=hist->GetXaxis()->GetXmin();

  if(x < histmin) x=histmin;
  if(x > histmax) x=histmax;

  if(interpolate){
    return hist->Interpolate(x);
  }
  else{
    const int xbin = hist->FindBin(x);
    return hist->GetBinContent(xbin);
  }
}

// max dt to look for previous mus
// Emily's official Li-9 likelihood time, scaled down to B-12
static const double maxdtmu = 0.7e9 * 20.20 / 178.3;

static float getlike(const float dist, const int nneut)
{
  //get the likelihood for the muon and event pair.

  // don't interpolate with neutron multiplicity because it's an int
  const float signn   = pdf(nneut, numNeutrons_sig, 0);
  const float bkgnn   = pdf(nneut, numNeutrons_bkg, 0);

  const float sigdist = pdf(dist,  muondist_sig, 1);
  const float bkgdist = pdf(dist,  muondist_bkg, 1);

  //value of signal joint-pdf
  float sigtot  = signn*sigdist;
  //value of background joint-pdf
  const float bkgtot  = bkgnn*bkgdist;

  // I think this number is pretty meaningless, since it is tuned for
  // the amount of li-9 and IBD, but maybe as long as the results are
  // reasonable spread between 0 and 1, it is ok.
  const float prior_ratio  = 1/(17.7 * 10.5 * maxdtmu * 1e-9);

  //multiply by prior probabilities
  sigtot *= prior_ratio;

  return (sigtot+bkgtot == 0) ? 0 : sigtot/(sigtot+bkgtot);
}

struct twolike{
  float like, altlike;
};

static twolike lb12like(const vector<track> & ts, const float x,
                        const float y, const float z, const float dtim)
{
  twolike max;
  max.like = max.altlike = 0;
  for(unsigned int i = 0; i < ts.size(); i++){
    const double dt = dtim - ts[i].tim;

    if(dt < 0) fprintf(stderr, "dt = %f\n", dt);

    const float dist = ptol(ts[i], x, y, z);
    const float like = getlike(dist, ts[i].nn);

    // Emily's likelihood ignores time except to exclude
    // muons that are too old
    if(dt < maxdtmu && like > max.like) max.like = like;

    // But what if instead we use the time in the likelihood?
    // There's no log here, right?
    const float altlike = like * exp(-dt/20.20e6);
    if(altlike > max.altlike) max.altlike = altlike;
  }
  return max;
}

double logis(const double x,
             const double p0, const double p1, const double p2)
{
  return p0/(1+ exp(-(x - p1)/p2));
}

// Return a michel/gamma/neutron energy adjusted for the baseline
// shift (presumably) after a muon.
double corrmiche(const double e, const double me)
{
  double pars[3] ={ 0.583, 287., 94.2};
  return e  - logis(me, pars[0], pars[1], pars[2]);
}

bool isnenergy(const double e, const double dt, const double me)
{
  const double corre = dt > 5500? e: corrmiche(e, me);

  return (corre > 1.8 && corre < 2.6) || (corre > 4.0 && corre < 10 );
}

static const char *qrms_name[2]      = { "qrms",      "QRMS"      };
static const char *ctmqtqall_name[2] = { "ctmqtqall", "Qratio"    };
static const char *ctrmsts_name[2]   = { "ctrmsts",   "RMSTstart" };
static const char *qdiff_name[2]     = { "qdiff",     "Qdiff"     };
static const char *ctEvisID_name[2]  = { "ctEvisID",  "EvisID"    };
static const char *ctq_name[2]       = { "ctq",       "ChargeID"  };
static const char *trgtime_name[2]   = { "trgtime",   "TrigTime"  };
static const char *run_name[2]       = { "run",       "RunNumber" };
static const char *coinov_name[2]    = { "coinov",    "NONONONO"  };
static const char *trgId_name[2]     = { "trgId",     "TriggerID" };
static const char *ctX_name[2]       = { "ctX",       "Vtx_BAMA"  };

#define TRANS(old, new) \
static inline void get_##old(TBranch * const br, const int i, \
                             const int whichname, dataparts & bits) \
{ \
  br->GetEntry(i); \
  if(whichname) bits.old = bits.new; \
}

TRANS(ctmqtqall, Qratio[0]); // or 1?
TRANS(ctrmsts, RMSTstart);
TRANS(ctEvisID, EvisID);
TRANS(qrms, QRMS);
TRANS(qdiff, Qdiff);
TRANS(run, RunNumber);

static void get_trgId(TBranch * br, const int i,
                      const int whichname, dataparts & bits,
                      const TTree * const fidotree)
{
  br->GetEntry(i);
  if(whichname) bits.trgId = bits.TriggerID;

  if(bits.trgId >= fidotree->GetEntries()){
    fprintf(stderr, "%s indexed FIDO entry %d, but only have %lld\n",
      whichname==0?"reduced":"JP", bits.trgId, fidotree->GetEntries());
    exit(1);
  }
}

static inline void get_ctX(TBranch * const br, const int i,
                           const int whichname, dataparts & bits)
{
  br->GetEntry(i);
  if(whichname)
    for(int i = 0; i < 3; i++)
      bits.ctX[i] = bits.Vtx_BAMA[i];
}

static int nnaftermu(const unsigned int muoni, dataparts & bits,
                      TTree * const chtree)
{
  const int whichname = !strcmp(chtree->GetName(), "GI");
  TBranch
    * const qrmsbr     = chtree->GetBranch(qrms_name[whichname]),
    * const ctmqtqallbr= chtree->GetBranch(ctmqtqall_name[whichname]),
    * const ctrmstsbr  = chtree->GetBranch(ctrmsts_name[whichname]),
    * const qdiffbr    = chtree->GetBranch(qdiff_name[whichname]),
    * const ctEvisIDbr = chtree->GetBranch(ctEvisID_name[whichname]),
    * const ctqbr      = chtree->GetBranch(ctq_name[whichname]),
    * const trgtimebr  = chtree->GetBranch(trgtime_name[whichname]);

  const double mindtmu      = 1e6;//dt to look for prev neutrons/veto

  trgtimebr->GetEntry(muoni);
  const double mutime = bits.trgtime;

  int found = 0;
  for(unsigned int i = muoni+1; i < chtree->GetEntries(); i++){
    trgtimebr->GetEntry(i);
    const double dt = bits.trgtime - mutime;
    if(dt > mindtmu) break;

    get_ctmqtqall(ctmqtqallbr, i, whichname, bits);
    get_ctrmsts  (ctrmstsbr,   i, whichname, bits);
    get_qdiff    (qdiffbr,     i, whichname, bits);
    get_qrms     (qrmsbr,      i, whichname, bits);

    // pass light noise
    if(lightnoise(bits.qrms, bits.ctmqtqall,
                  bits.ctrmsts, bits.qdiff)) continue;

    get_ctEvisID(ctEvisIDbr, i, whichname, bits);
    if(bits.ctEvisID == 0){
      ctqbr->GetEntry(i);
      bits.ctEvisID = bits.ctq/34e3;
    }

    // right energy for a neutron capture
    if(!isnenergy(bits.ctEvisID, 1e4, 0 /* no correction */)) continue;

    found++;
  }
  return found;
}

static void searchfrommuon(dataparts & bits, TTree * const chtree,
                           TTree * const fitree,
                           const unsigned int muoni,
                           const searchtype search,
                           const double timeleft)
{
  const double mutime = bits.trgtime;
  const double entr_mux = bits.ids_entr_x,
               entr_muy = bits.ids_entr_y,
               entr_muz = bits.ids_entr_z;
  const double mux = near?bits.ids_end_x:fidocorrx(bits.ids_end_x),
               muy = near?bits.ids_end_y:fidocorry(bits.ids_end_y),
               muz = near?bits.ids_end_z:fidocorrz(bits.ids_end_z);
  const float gclen = bits.ids_gclen;
  const float idexitqf = bits.id_idexitqf;
  const float ivqbal = bits.id_ivqbal;
  const int mutrgid = bits.trgId;
  const int murun = bits.run;
  const bool mucoinov = bits.coinov;
  const float mudchi2 = bits.ids_chi2 - bits.id_chi2;
  const float murchi2 =
    bits.ids_chi2/(bits.nidtubes+bits.nivtubes-6);

  // This will get used for the high-purity cut in downstream anslysis.
  // Note that this uses the ids_ variables, which is *different*
  // from the cut used to select in searchforamuon(), i.e. the
  // loose cut. This is aesthetically unpleasing, but perfectly
  // valid statistically. Changing them both to be one or the other
  // (and therefore aesthetically pleasing) and then retuning the
  // high-purity cuts would not give a better result, but it would give
  // a *different* result, and the current methods and results are
  // already blessed.
  const float muivdedx = bits.ids_ivlen?
    bits.fido_qiv/(bits.ids_ivlen-bits.ids_buflen): 0;

  const float mufqid = bits.fido_qid,
              mufqiv = bits.fido_qiv,
              muctqid = bits.ctq,
              muctqiv = bits.ctqIV;

  unsigned int nneutronanydist[2] = {0,0}, ngdneutronanydist[2] = {0,0};
  unsigned int nneutronnear[2] = {0,0}, ngdneutronnear[2] = {0,0};

  double michelt = 0, michelx = 0, michely = 0, michelz = 0,
         michele = 0, michdist = 0;
  bool followingov = false;
  double followingovtime = 0, followingqivtime = 0;
  double firstlatenearneutrontime = 0, firstlatenearneutronenergy = 0;
  double firstneutrontime = 0, firstneutronenergy = 0;
  float followingqiv = 0;

  const int whichname = !strcmp(chtree->GetName(), "GI");
  TBranch
    * const nidtubesbr      = fitree->GetBranch("nidtubes"),
    * const fido_qivbr      = fitree->GetBranch("fido_qiv"),
    * const fido_didfitbr   = fitree->GetBranch("id_didfit"),
    * const fido_entrxbr    = fitree->GetBranch("id_entr_x"),
    * const fido_entrybr    = fitree->GetBranch("id_entr_y"),
    * const fido_entrzbr    = fitree->GetBranch("id_entr_z"),
    * const fido_endxbr     = fitree->GetBranch("id_end_x"),
    * const fido_endybr     = fitree->GetBranch("id_end_y"),
    * const fido_endzbr     = fitree->GetBranch("id_end_z"),

    * const hambr           = chtree->GetBranch("Trk_MuHamID"),
    * const runbr           = chtree->GetBranch(run_name[whichname]),
    * const coinovbr        = chtree->GetBranch(coinov_name[whichname]),
    * const trgIdbr         = chtree->GetBranch(trgId_name[whichname]),
    * const ctXbr           = chtree->GetBranch(ctX_name[whichname]),
    * const qrmsbr          = chtree->GetBranch(qrms_name[whichname]),
    * const ctmqtqallbr     = chtree->GetBranch(ctmqtqall_name[whichname]),
    * const ctrmstsbr       = chtree->GetBranch(ctrmsts_name[whichname]),
    * const qdiffbr         = chtree->GetBranch(qdiff_name[whichname]),
    * const ctEvisIDbr      = chtree->GetBranch(ctEvisID_name[whichname]),
    * const ctqbr           = chtree->GetBranch(ctq_name[whichname]),
    * const trgtimebr       = chtree->GetBranch(trgtime_name[whichname]);

  const double max_micht = near?30000:5500;
  const double max_time_probably_a_mich = 5500;


  vector<double> ntoprint_e, ntoprint_t;
  vector<cart> ntoprint_x;
  vector<int> ntoprint_ov;
  for(unsigned int i = muoni+1; i < chtree->GetEntries(); i++){

    trgtimebr->GetEntry(i);
    if(coinovbr) coinovbr->GetEntry(i);
    const double dt = bits.trgtime - mutime;

    // For any uncut Michels or (hopefully) prompt gammas from muon
    // capture, record the time and energy. The *highest* energy
    // event in the time window is accepted. I'm going to allow IV
    // energy since these may be very close to the muon event. Note that
    // we will often count this Michel as a neutron also for the zeroth
    // element of the count arrays (but not the first), so don't double
    // count by accident.
    if(dt < max_micht && !bits.coinov){
      get_ctEvisID(ctEvisIDbr, i, whichname, bits);
      if(bits.ctEvisID == 0){
        ctqbr->GetEntry(i);
        bits.ctEvisID = bits.ctq/34e3;
      }
      get_ctX(ctXbr, i, whichname, bits);
      if(bits.ctEvisID > michele){
        michele = bits.ctEvisID, michelt = dt;
        michelx = bamacorrxy(bits.ctX[0], bits.ctEvisID);
        michely = bamacorrxy(bits.ctX[1], bits.ctEvisID);
        michelz = bamacorrxy(bits.ctX[2], bits.ctEvisID);
        michdist = sqrt(pow(mux-michelx, 2) +
                        pow(muy-michely, 2) +
                        pow(muz-michelz, 2));
      }
    }

    if(dt < max_micht){
      get_trgId(trgIdbr, i, whichname, bits, fitree);

      fido_qivbr->GetEntry(bits.trgId);
      if(bits.coinov){
         followingov = true;
         followingovtime = dt;
      }
      if(bits.fido_qiv > followingqiv){
        followingqiv = bits.fido_qiv;
        followingqivtime = dt;
      }
    }

    // Skip past retriggers and whatnot, like DC3rdPub
    if(dt < 500) continue;

    // Not more than ~4 nH lifetimes
    if(dt > 800e3) break;

    get_ctmqtqall(ctmqtqallbr, i, whichname, bits);
    get_ctrmsts  (ctrmstsbr,   i, whichname, bits);
    get_qdiff    (qdiffbr,     i, whichname, bits);
    get_qrms     (qrmsbr,      i, whichname, bits);

    // pass light noise
    if(lightnoise(bits.qrms, bits.ctmqtqall,
                  bits.ctrmsts, bits.qdiff)) continue;

    get_ctEvisID(ctEvisIDbr, i, whichname, bits);
    if(bits.ctEvisID == 0){
      ctqbr->GetEntry(i);
      bits.ctEvisID = bits.ctq/34e3;
    }
    // right energy for a neutron capture
    if(!isnenergy(bits.ctEvisID, dt, mufqid/8300)) continue;

    get_ctX(ctXbr, i, whichname, bits);

    // Note, not 1000mm, but 800mm.
    const bool nnear =
      sqrt(pow(mux - bamacorrxy(bits.ctX[0], bits.ctEvisID), 2)+
           pow(muy - bamacorrxy(bits.ctX[1], bits.ctEvisID), 2)+
           pow(muz - bamacorrz( bits.ctX[2], bits.ctEvisID), 2)) < 800;

    const bool gd = bits.ctEvisID > 4.0 && bits.ctEvisID < 10
                 && bits.trgtime - mutime < 150e3;

    if(search == neutron){
      ntoprint_e.push_back(bits.ctEvisID);
      ntoprint_t.push_back(dt);
      ntoprint_ov.push_back(bits.coinov);
      ntoprint_x.push_back(cart(bamacorrxy(bits.ctX[0],bits.ctEvisID),
                                bamacorrxy(bits.ctX[1],bits.ctEvisID),
                                bamacorrz( bits.ctX[2],bits.ctEvisID)));
    }

    const bool alsoamichel = dt < max_time_probably_a_mich;
    nneutronanydist[0]++;
    if(firstneutrontime == 0){
      firstneutrontime = dt;
      firstneutronenergy = bits.ctEvisID;
    }
    if(gd) ngdneutronanydist[0]++;
    if(nnear) nneutronnear[0]++;
    if(gd && nnear) ngdneutronnear[0]++;
    if(!alsoamichel){
      nneutronanydist[1]++;
      if(gd) ngdneutronanydist[1]++;
      if(nnear){
        // special case since these have turned out to be the most
        // useful count of neutrons. Record the time and energy of the
        // first one
        nneutronnear[1]++;
        if(firstlatenearneutrontime == 0){
          firstlatenearneutrontime = dt;
          firstlatenearneutronenergy = bits.ctEvisID;
        }
      }
      if(gd && nnear) ngdneutronnear[1]++;
    }
  }

  for(unsigned int i = 0; i < ntoprint_e.size(); i++)
    printf("%u %d %lf %lf %lf %d %lf %lf %lf %lf\n", i,
           (int)ntoprint_e.size(), mufqid, ntoprint_t[i], ntoprint_e[i],
           ntoprint_ov[i], mufqiv,
           ntoprint_x[i].x, ntoprint_x[i].y, ntoprint_x[i].z);

  unsigned int got = 0, printed = 0;

  // positions of putative isotope decays
  double ix[2], iy[2], iz[2];

  double lastmuontime = mutime, lastgcmuontime = mutime,
         lastvalidtime = mutime;

  vector<track> tmuons;

  for(unsigned int i = muoni+1; i < chtree->GetEntries(); i++){

    trgtimebr->GetEntry(i);
    get_trgId(trgIdbr, i, whichname, bits, fitree);

    const double itime = bits.trgtime;
    const double dt_ms = (itime - mutime)/1e6;
    const double ttlastvalid =(itime-lastvalidtime  )/1e6;
    const double ttlastmuon  =(itime-lastmuontime   )/1e6;
    const double ttlastgcmuon=(itime-lastgcmuontime )/1e6;

     // NOTE-luckplan
     #define LATEFORM \
     "%d %d %d " \
     "%f %f %f %f %f %f " \
     "%d %d %d %d %d %d %d %d " \
     "%lf %.1f %.1f %.1f %.0lf %f %f %f %f %.0f " \
     "%f %f %f %f %f %f %f %f %d %f %f %f %d %f %f %f %f %f %f"

     #define LATEVARS \
     murun, mutrgid, mucoinov, \
     mux, muy, muz, mudchi2, murchi2, muivdedx, \
     ngdneutronnear[0], ngdneutronanydist[0], \
     nneutronnear[0],  nneutronanydist[0], \
     ngdneutronnear[1], ngdneutronanydist[1], \
     nneutronnear[1],  nneutronanydist[1], \
     michele, michelx, michely, michelz, \
     michelt, gclen, entr_mux, entr_muy, entr_muz, \
     michdist, \
     mufqid, mufqiv, muctqid, muctqiv, \
     timeleft, ttlastvalid, ttlastmuon, \
     ttlastgcmuon, followingov, followingovtime, \
     followingqiv, followingqivtime, printed, \
     firstlatenearneutrontime, firstlatenearneutronenergy, \
     firstneutrontime, firstneutronenergy, \
     idexitqf, ivqbal

    if(dt_ms > maxtime){ // stop looking and print muon info
      // This is just overhead for the neutron and be12 searches,
      // so don't do it.
      if(search != neutron && search != be12 && printed == 0)
        printf("0 0 0 0 0 0 0 0 0 0 " LATEFORM "\n", LATEVARS);
      break;
    }

    // Require at least 500us since the last muon so we don't count
    // neutrons as isotope decays
    if(ttlastmuon < 0.5) goto end;

    if(dt_ms < 1) goto end; // Go past all H neutron captures

    get_ctEvisID(ctEvisIDbr, i, whichname, bits);
    if(bits.ctEvisID == 0){
      ctqbr->GetEntry(i);
      bits.ctEvisID = bits.ctq/34e3;
    }

    // Ignore low energy accidentals and H-neutrons
    if(bits.ctEvisID < minenergy) goto end;

    // Ignore events above the end point + res
    if(bits.ctEvisID > maxenergy) goto end;

    fido_qivbr->GetEntry(bits.trgId);
    // No IV, OV energy
    if(bits.fido_qiv > fido_qiv_muon_def) goto end;

    if(coinovbr) coinovbr->GetEntry(i);
    if(bits.coinov) goto end;

    get_ctmqtqall(ctmqtqallbr, i, whichname, bits);
    get_ctrmsts  (ctrmstsbr,   i, whichname, bits);
    get_qdiff    (qdiffbr,     i, whichname, bits);
    get_qrms     (qrmsbr,      i, whichname, bits);

    // XXX disable for now for ND, but probably want to put back later
    if(!near){
      if(lightnoise(bits.qrms, bits.ctmqtqall, bits.ctrmsts, bits.qdiff))
        goto end;
    }

    get_ctX(ctXbr, i, whichname, bits);

    ix[got] = bamacorrxy(bits.ctX[0], bits.ctEvisID),
    iy[got] = bamacorrxy(bits.ctX[1], bits.ctEvisID),
    iz[got] = bamacorrz( bits.ctX[2], bits.ctEvisID);

    {
      // If the be12 search, record distance to previous selected event,
      // not always to muon, otherwise always to the muon.
      const double dist = (got == 0 || search != be12)?
        sqrt(pow(mux  -ix[0],2)+pow(muy  -iy[0],2)+pow(muz  -iz[0],2)):
        sqrt(pow(ix[0]-ix[1],2)+pow(iy[0]-iy[1],2)+pow(iz[0]-iz[1],2));

      // And, optionally, they must be near each other.
      if(distcut != 0 && dist > distcut) goto end;

      get_run(runbr, i, whichname, bits);
      get_trgId(trgIdbr, i, whichname, bits, fitree);

      const twolike b12like =
        lb12like(tmuons, ix[got], iy[got], iz[got], bits.trgtime);
      if(search != neutron)
        // NOTE-luckplan
        printf("%d %lf %f %f %f %f %f %f %f %f " LATEFORM "%c",
               bits.trgId, dt_ms, dist,
               bits.ctEvisID, bits.fido_qiv, ix[got], iy[got], iz[got],
               b12like.like, b12like.altlike,
               LATEVARS,
               search == be12?' ':'\n');
      printed++;

      // If searching for be12, use the ix/iy/iz arrays and stop when we
      // get 2 decays; these are printed on the same line. Otherwise,
      // keep going up to the time limit and put each on its own line.
      if(search == be12) got++;

      if(got >= 2) break;
    }

    end:

    // Note the time of this event if it is a muon.
    if(coinovbr) coinovbr->GetEntry(i);
    fido_qivbr->GetEntry(bits.trgId);
    fido_didfitbr->GetEntry(bits.trgId);
    get_ctEvisID(ctEvisIDbr, i, whichname, bits);
    if(bits.ctEvisID == 0){
      ctqbr->GetEntry(i);
      bits.ctEvisID = bits.ctq/34e3;
    }

    // NOTE-lungbloke: See comment at other.
    if(bits.coinov || bits.fido_qiv > fido_qiv_muon_def
      || bits.ctEvisID > 60)
      lastmuontime = bits.trgtime;

    // Note the time of this even if it is valid, which for me means
    // not light noise and at least 0.4 MeV
    get_ctmqtqall(ctmqtqallbr, i, whichname, bits);
    get_ctrmsts  (ctrmstsbr,   i, whichname, bits);
    get_qdiff    (qdiffbr,     i, whichname, bits);
    get_qrms     (qrmsbr,      i, whichname, bits);

    trgtimebr->GetEntry(i);

    if(!lightnoise(bits.qrms, bits.ctmqtqall, bits.ctrmsts, bits.qdiff)
       && bits.ctEvisID > 0.4)
      lastvalidtime = bits.trgtime;

    fido_entrzbr->GetEntry(bits.trgId);
    fido_endzbr ->GetEntry(bits.trgId);
    nidtubesbr->GetEntry(bits.trgId);
    if(bits.id_didfit && bits.nidtubes > 30 &&
       bits.id_entr_z > bits.id_end_z){
      lastgcmuontime = bits.trgtime;

      fido_entrxbr->GetEntry(bits.trgId);
      fido_entrybr->GetEntry(bits.trgId);
      fido_endxbr ->GetEntry(bits.trgId);
      fido_endybr ->GetEntry(bits.trgId);

      const double mutime = bits.trgtime;

      tmuons.push_back(maketrack(
        bits.id_entr_x, bits.id_entr_y, bits.id_entr_z,
        bits.id_end_x,  bits.id_end_y,  bits.id_end_z,
        mutime, nnaftermu(i, bits, chtree))); // warning: changes bits
    }
    else if(hambr){
      // If FIDO through-going tracks are not available, but Ham is,
      // use it to get the B-12 likelihood.
      hambr->GetEntry(bits.trgId);
      get_ctEvisID(ctEvisIDbr, i, whichname, bits);

      if(bits.Trk_MuHamID[0][0] != 0 && bits.ctEvisID > 60 &&
         bits.Trk_MuHamID[0][2] > bits.Trk_MuHamID[1][2]){
        lastgcmuontime = bits.trgtime;
        const double mutime = bits.trgtime;

        tmuons.push_back(maketrack(
          bits.Trk_MuHamID[0][0], bits.Trk_MuHamID[0][1],
                                  bits.Trk_MuHamID[0][2],
          bits.Trk_MuHamID[1][0], bits.Trk_MuHamID[1][1],
                                  bits.Trk_MuHamID[1][2],
          mutime, nnaftermu(i, bits, chtree))); // warning: changes bits
      }
    }

    // NOTE-luckplan
    // at EOR, be sure to print muon info
    if(search != neutron && i == chtree->GetEntries()-1 && printed == 0)
      printf("0 0 0 0 0 0 0 0 0 0 " LATEFORM "\n", LATEVARS);
  }
  if(got){
    if(search != neutron) printf("\n");
    fflush(stdout);
  }
}

/* Gets the time of the end of the run */
static double geteortime(dataparts & parts, TTree * const chtree,
                         const unsigned int start)
{
  const int whichname = !strcmp(chtree->GetName(), "GI");
  TBranch * const runbr = chtree->GetBranch(run_name[whichname]);
  TBranch * const timbr = chtree->GetBranch(trgtime_name[whichname]);

  get_run(runbr, start, whichname, parts);
  const int run = parts.run;
  double eortime = 0;

  for(int e = start; e < chtree->GetEntries(); e++){
    get_run(runbr, e, whichname, parts);
    if(run == parts.run){
      timbr->GetEntry(e);
      eortime = parts.trgtime;
    }
    else{
      break;
    }
  }

  return eortime;
}

static void searchforamuon(dataparts & parts, TTree * const chtree,
                           TTree * const fitree,
                           const searchtype search)
{
  double lastmuontime = 0;
  double eortime = 0;
  int run = 0;

  const int whichname = !strcmp(chtree->GetName(), "GI");
  TBranch * const runbr     = chtree->GetBranch(run_name[whichname]);
  TBranch * const trgIdbr   = chtree->GetBranch(trgId_name[whichname]);
  TBranch * const trgtimebr = chtree->GetBranch(trgtime_name[whichname]);
  TBranch * const coinovbr  = chtree->GetBranch(coinov_name[whichname]);

  TBranch * const nidtubesbr   = fitree->GetBranch("nidtubes");
  TBranch * const nivtubesbr   = fitree->GetBranch("nivtubes");
  TBranch * const fido_qivbr   = fitree->GetBranch("fido_qiv");
  TBranch * const fido_qidbr   = fitree->GetBranch("fido_qid");
  TBranch * const ids_didfitbr = fitree->GetBranch("ids_didfit");
  TBranch * const id_buflenbr  = fitree->GetBranch("id_buflen");
  TBranch * const id_chi2br    = fitree->GetBranch("id_chi2");
  TBranch * const id_entr_xbr  = fitree->GetBranch("id_entr_x");
  TBranch * const id_entr_ybr  = fitree->GetBranch("id_entr_y");
  TBranch * const id_ivlenbr   = fitree->GetBranch("id_ivlen");
  TBranch * const ids_chi2br   = fitree->GetBranch("ids_chi2");
  TBranch * const ids_end_xbr  = fitree->GetBranch("ids_end_x");
  TBranch * const ids_end_ybr  = fitree->GetBranch("ids_end_y");
  TBranch * const ids_end_zbr  = fitree->GetBranch("ids_end_z");
  TBranch * const ids_entr_xbr = fitree->GetBranch("ids_entr_x");
  TBranch * const ids_entr_ybr = fitree->GetBranch("ids_entr_y");
  TBranch * const ids_entr_zbr = fitree->GetBranch("ids_entr_z");

  for(unsigned int mi = 0; mi < chtree->GetEntries()-1 &&
                           mi < fitree->GetEntries()-1; mi++){
    get_run(runbr, mi, whichname, parts);
    get_trgId(trgIdbr, mi, whichname, parts, fitree);

    if(parts.run != run){
      run = parts.run;
      lastmuontime = 0; // Assume a muon right before the run started.
      eortime = geteortime(parts, chtree, mi);
    }

    ids_didfitbr->GetEntry(parts.trgId);

    // Must be possible that it's a stopper (or not for neutron search)
    if(!((search == neutron) ^ parts.ids_didfit)) goto end;

    trgtimebr->GetEntry(mi);

    // Require at least 500us since the last muon so we don't count
    // neutrons from the previous one as belonging to this one.
    if(parts.trgtime - lastmuontime < 500e3) goto end;

    fido_qivbr->GetEntry(parts.trgId);
    if(parts.fido_qiv < fido_qiv_muon_def) goto end;

    fido_qidbr->GetEntry(parts.trgId);
    if(parts.fido_qid/(near?13500:8300) > 700) goto end;
    if(parts.fido_qid/(near?13500:8300) < (search == neutron?1:60)) goto end;

    nivtubesbr->GetEntry(parts.trgId);
    nidtubesbr->GetEntry(parts.trgId);

    if(parts.nidtubes+parts.nivtubes <= 6) goto end;

    if(search != neutron){
      ids_chi2br->GetEntry(parts.trgId);
      if(!near){
        if(parts.ids_chi2/(parts.nidtubes+parts.nivtubes-6) > 10) goto end;
      }

      // These correctly use the *uncorrected* position
      if(search == buffer){
        // Accept *only* if the muon clearly appears to be exiting
        ids_end_xbr->GetEntry(parts.trgId);
        ids_end_ybr->GetEntry(parts.trgId);
        ids_end_zbr->GetEntry(parts.trgId);
        if(pow(parts.ids_end_x,2)+pow(parts.ids_end_y,2)<pow(1708-35,2)
            &&
           parts.ids_end_z > -1786+(near?40:35)) goto end;
      }
      else{
        // Do *not* accept if the muon appears to be exiting
        ids_end_xbr->GetEntry(parts.trgId);
        ids_end_ybr->GetEntry(parts.trgId);
        if(pow(parts.ids_end_x,2)+pow(parts.ids_end_y,2)>pow(1708-35,2))
          goto end;

        ids_end_zbr->GetEntry(parts.trgId);
        if(parts.ids_end_z < -1786+(near?40:35)) goto end;
      }

      ids_entr_zbr->GetEntry(parts.trgId);
      id_ivlenbr->GetEntry(parts.trgId);
      id_buflenbr->GetEntry(parts.trgId);
      if(parts.ids_entr_z > 11500 -(near?37:62)*
         parts.fido_qiv/(parts.id_ivlen-parts.id_buflen)) goto end;

      // disable for now for ND. Do I want to put it back? Currently, it
      // looks like this has no power to tell a good event from a bad
      // one, but, well, I'm writing it out so I can always look again.
      if(!near){
        id_chi2br->GetEntry(parts.trgId);
        ids_chi2br->GetEntry(parts.trgId);
        if(parts.ids_chi2-parts.id_chi2 > 800) goto end;
      }

      // XXX disable experimentally for ND
      if(!near){
        id_entr_xbr->GetEntry(parts.trgId);
        id_entr_ybr->GetEntry(parts.trgId);
        ids_entr_xbr->GetEntry(parts.trgId);
        ids_entr_ybr->GetEntry(parts.trgId);
        if(pow(parts.id_entr_x, 2)+pow(parts.id_entr_y, 2) < pow(1000, 2) &&
           pow(parts.ids_entr_x,2)+pow(parts.ids_entr_y,2) > pow(2758, 2)) goto end;
      }
    }

    chtree->GetEntry(mi);
    fitree->GetEntry(parts.trgId);

    searchfrommuon(parts, chtree, fitree, mi, search,
    // Record how much time is left to the end of the run so that we can
    // exclude events that are too close when analyzing.
                             (eortime - parts.trgtime)/1e6);
    end:

    // See note above for same code.
    if(coinovbr) coinovbr->GetEntry(mi);
    fido_qivbr->GetEntry(parts.trgId);

    // Note-lungbloke: Notice how this cut (which vetos muons) is
    // slightly different from the other one at lungbloke (which vetos
    // beta decays). In the tech note, the additional cut for > 60 MeV
    // in the ID isn't mentioned. Fortunately, it makes very close to
    // no difference (3e-8 -- 0.04%, depending on run) since almost all
    // events with > 60 MeV in the ID either have an OV coincidenece or
    // some IV energy.
    if(parts.coinov || parts.fido_qiv > fido_qiv_muon_def){
      trgtimebr->GetEntry(mi);
      lastmuontime = parts.trgtime;
    }
  }
}

int main(int argc, char ** argv)
{
  if(argc < 8 || argc%2 != 0){
    fprintf(stderr,
      "b12search (near|far) distcut maxtime[ms] minenergy maxenergy "
      "[reduced or JP file, fido file]*\n\n");

    fprintf(stderr,
      "To run the be12 search, looking for exactly two decays:\n"
      "be12search (near|far) distcut maxtime[ms] minenergy maxenergy "
      "[reduced or JP file, fido file]*\n\n");

    fprintf(stderr,
      "To accept throughgoing muons for neutron studies:\n"
      "neutronsearch (near|far) distcut maxtime[ms] minenergy maxenergy "
      "[reduced or JP file, fido file]*\n\n");

    fprintf(stderr,
      "To search for buffer stopping muons:\n"
      "buffersearch (near|far) distcut maxtime[ms] minenergy maxenergy "
      "[reduced or JP file, fido file]*\n\n");

    fprintf(stderr,
      "To check files only, use\n"
      "checkb12search (near|far) 0 foo foo foo "
      "[reduced or JP file, fido file]*\n");

    fprintf(stderr, "\nA distcut of zero means unlimited\n");

    exit(1);
  }

  if(!strcmp(argv[1], "near")) near = true;
  else if(!strcmp(argv[1], "far")) near = false;
  else{
    fprintf(stderr, "Give \"near\" or \"far\" as the first argument\n");
    exit(1);
  }

  distcut = atof(argv[2]);
  maxtime = atof(argv[3]);
  minenergy = atof(argv[4]);
  maxenergy = atof(argv[5]);

  const searchtype search =
    !strcmp(basename(argv[0]),    "be12search")? be12:
    !strcmp(basename(argv[0]), "neutronsearch")? neutron:
    !strcmp(basename(argv[0]),  "buffersearch")? buffer:
                                                 b12;

  gErrorIgnoreLevel = kFatal;
  int errcode = 0;

  if(search == neutron)
    printf("i/I:nn/I:fq/F:t/F:e/F:ov/I:fqiv/F:x/F:y/F:z/F");
  // NOTE-luckplan
  else printf(
    "trig/I:"
    "dt/F:"
    "dist/F:"
    "e/F:"
    "eiv/F:"
    "dx/F:"
    "dy/F:"
    "dz/F:"
    "b12like/F:"
    "b12altlike/F:"
    "run/I:"
    "mutrig/I:"
    "ovcoin/I:"
    "mx/F:"
    "my/F:"
    "mz/F:"
    "dchi2/F:"
    "rchi2/F:"
    "ivdedx/F:"
    "ngdnear/I:"
    "ngd/I:"
    "nnear/I:"
    "n/I:"
    "latengdnear/I:"
    "latengd/I:"
    "latennear/I:"
    "laten/I:"
    "miche/F:"
    "michx/F:"
    "michy/F:"
    "michz/F:"
    "micht/F:"
    "gclen/F:"
    "fex/F:"
    "fey/F:"
    "fez/F:"
    "michd/F:"
    "fq/F:"
    "fqiv/F:"
    "cq/F:"
    "cqiv/F:"
    "timeleft/F:"
    "ttlastvalid/F:"
    "ttlastmuon/F:"
    "ttlastgcmuon/F:"
    "followingov/O:"
    "followingovtime/F:"
    "followingqiv/F:"
    "followingqivtime/F:"
    "ndecay/I:"
    "firstlatenneart/F:"
    "firstlatenneare/F:"
    "firstnt/F:"
    "firstne/F:"
    "idexitqf/F:"
    "ivqbal/F");
  if(search == be12) // same as above with 2s appended to each name
    printf(          // some are dumb, since, i.e., mutrig === mutrig2
    ":trig2/I:"
    "dt2/F:"
    "dist2/F:"
    "e2/F:"
    "eiv2/F:"
    "dx2/F:"
    "dy2/F:"
    "dz2/F:"
    "b12like2/F:"
    "b12altlike2/F:"
    "run2/I:"
    "mutrig2/I:"
    "ovcoin2/I:"
    "mx2/F:"
    "my2/F:"
    "mz2/F:"
    "dchi2/F:"
    "rchi22/F:"
    "ivdedx2/F:"
    "ngdnear2/I:"
    "ngd2/I:"
    "nnear2/I:"
    "n2/I:"
    "latengdnear2/I:"
    "latengd2/I:"
    "latennear2/I:"
    "laten2/I:"
    "miche2/F:"
    "michx2/F:"
    "michy2/F:"
    "michz2/F:"
    "micht2/F:"
    "gclen2/F:"
    "fex2/F:"
    "fey2/F:"
    "fez2/F:"
    "michd2/F:"
    "fq2/F:"
    "fqiv2/F:"
    "cq2/F:"
    "cqiv2/F:"
    "timeleft2/F:"
    "ttlastvalid2/F:"
    "ttlastmuon2/F:"
    "ttlastgcmuon2/F:"
    "followingov2/O:"
    "followingovtime2/F:"
    "followingqiv2/F:"
    "followingqivtime2/F:"
    "ndecay2/I:"
    "firstlatenneart2/F:"
    "firstlatenneare2/F:"
    "firstnt2/F:"
    "firstne2/F:"
    "idexitqf2/F:"
    "ivqbal2/F");
  printf("\n");

  for(int i = 6; i < argc; i+=2){
    TFile * const chfile = new TFile(argv[i], "read");
    TFile * const fifile = new TFile(argv[i+1], "read");
    dataparts parts;
    TTree * chtree = NULL, * fitree = NULL;

    if(!chfile || chfile->IsZombie()){
      fprintf(stderr, "\nI couldn't read %s\n", argv[i]);
      errcode |= 1;
      goto cleanup;
    }

    chtree = (TTree *)chfile->Get("data");
    if(!chtree) chtree = (TTree *)chfile->Get("GI");
    if(!chtree){
      fprintf(stderr,"\n%s lacks a \"data\" or \"GI\" tree!\n",argv[i]);
      errcode |= 2;
      goto cleanup;
    }

    if(!fifile || fifile->IsZombie()){
      fprintf(stderr, "\nI couldn't read %s\n", argv[i+1]);
      errcode |= 4;
      goto cleanup;
    }

    fitree = (TTree *)fifile->Get("RecoMuonFIDOInfoTree");
    if(!fitree) fitree = (TTree *)fifile->Get("fido");

    if(!fitree){
      fprintf(stderr, "\n%s lacks a fido or RecoMuonFIDOInfoTree!\n",
              argv[i+1]);
      errcode |= 8;
      goto cleanup;
    }

    if(!strcmp(chtree->GetName(), "data")){
      if(chtree->GetEntries() != fitree->GetEntries()){
        fprintf(stderr,
          "\n%s,\n%s:\nreduced has %ld entries, fido %ld (expect same)\n",
          argv[i], argv[i+1],
          long(chtree->GetEntries()), long(fitree->GetEntries()));
        errcode |= 0x10;
      }
    }
    else{ // JP tree is expected to have fewer entries
      if(chtree->GetEntries() > fitree->GetEntries()){
        fprintf(stderr,
          "\n%s,\n%s:\nJP has more entries (%ld) than fido (%ld)\n",
          argv[i], argv[i+1],
          long(chtree->GetEntries()), long(fitree->GetEntries()));
        errcode |= 0x10;
      }
    }

    if(!strcmp(basename(argv[0]), "checkb12search")) goto cleanup;

    fitree->SetMakeClass(1);

    for(unsigned int i = 0; i < noff; i++)
      chtree->SetBranchStatus(turnoff[i], 0);

    fitree->SetBranchStatus("*", 0);

    // Tried different settings, this is the best
    fitree->SetCacheSize(10000000);
    chtree->SetCacheSize(10000000);
    chtree->AddBranchToCache("*");

    memset(&parts, 0, sizeof(parts));

    /******************************************************************/
    /*                   Set up FIDO trees for reading                */
    /******************************************************************/

    #define fSBA(x) fitree->SetBranchStatus(#x, 1); \
                    fitree->SetBranchAddress(#x, &parts.x); \
                    fitree->AddBranchToCache(#x);

    fSBA(ids_didfit);

    fSBA(ids_end_x);
    fSBA(ids_end_y);
    fSBA(ids_end_z);
    fSBA(ids_entr_x);
    fSBA(ids_entr_y);
    fSBA(ids_entr_z);

    fSBA(id_entr_x);
    fSBA(id_entr_y);
    fSBA(id_entr_z);
    fSBA(id_end_x);
    fSBA(id_end_y);
    fSBA(id_end_z);

    fSBA(fido_qid);
    fSBA(fido_qiv);
    fSBA(id_buflen);
    fSBA(id_chi2);
    fSBA(id_didfit);
    fSBA(id_ivlen);
    fSBA(ids_buflen);
    fSBA(ids_chi2);
    fSBA(ids_gclen);
    fSBA(ids_ivlen);
    fSBA(ids_phi);
    fSBA(ids_theta);
    fSBA(nidtubes);
    fSBA(nivtubes);

    fSBA(id_idexitqf);
    fSBA(id_ivqbal);

    /******************************************************************/
    /*   Set up reduced or JP trees for reading, whichever we have    */
    /******************************************************************/

    // if types are compatible, write to same part either way
    #define cSBA(x,y) chtree->SetBranchAddress(\
                        !strcmp(chtree->GetName(), "data")?#x:#y,\
                        &parts.x);

    // Usually have to write to different parts for reduced and JP
    // since they are stored in incompatible types.
    #define cSBA2(x,y) \
      if(!strcmp(chtree->GetName(), "data")) \
        chtree->SetBranchAddress(#x, &parts.x); \
      else \
        chtree->SetBranchAddress(#y, &parts.y);

    // In the JP trees, ChargeID and ChargeIV are arrays of length 3.
    // Arrange for the third element to be in ctq.
    chtree->SetBranchAddress(
      !strcmp(chtree->GetName(), "data")?"ctq":"ChargeID",
      !strcmp(chtree->GetName(), "data")?&parts.ctq:&parts.ctq0);
    chtree->SetBranchAddress(
      !strcmp(chtree->GetName(), "data")?"ctq":"ChargeIV",
      !strcmp(chtree->GetName(), "data")?&parts.ctqIV:&parts.ctqIV0);

    // We might have both this and the FIDO tracking available, but most
    // likely I will not have run FIDO on through-going muons, so Ham is
    // the only through-going information.
    if(strcmp(chtree->GetName(), "data"))
      chtree->SetBranchAddress("Trk_MuHamID", &parts.Trk_MuHamID);

    cSBA(trgtime, TrigTime);

    cSBA2(run, RunNumber);
    cSBA2(trgId, TriggerID);

    // Can't use cSBA2 because one is an array and the other isn't
    if(!strcmp(chtree->GetName(), "data"))
      chtree->SetBranchAddress("ctmqtqall", &parts.ctmqtqall);
    else
      chtree->SetBranchAddress("Qratio", parts.Qratio);

    cSBA2(ctrmsts, RMSTstart);

    // Just not in JP ntuples, but later will be, and will be an
    // unsigned short called OVTrigCoin, as compared to coinov, a bool.
    if(!strcmp(chtree->GetName(), "data"))
      chtree->SetBranchAddress("coinov", &parts.coinov);

    cSBA2(ctEvisID, EvisID);

    cSBA2(qrms, QRMS);
    cSBA2(qdiff, Qdiff);

    if(!strcmp(chtree->GetName(), "data")){
      if(chtree->SetBranchAddress("ctX", parts.ctX) < 0)
        chtree->SetBranchAddress("ctX[3]", parts.ctX);
    }
    else{
      chtree->SetBranchAddress("Vtx_BAMA", parts.Vtx_BAMA);
    }

    searchforamuon(parts, chtree, fitree, search);

    cleanup:

    if(fitree) delete fitree;
    if(fifile) delete fifile;
    if(chtree) delete chtree;
    if(chfile) delete chfile;
  }
  return errcode;
}
