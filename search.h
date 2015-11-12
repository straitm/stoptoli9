bool lightnoise(const float qrms, const float mqtq,
                const float rmsts, const float qdiff);

double fidocorrz(const double z);
double fidocorrx(const double x);
double fidocorry(const double y);

double bamacorrxy(const double xy, const double e);
double bamacorrz (const double  z, const double e);

struct dataparts{
  bool coinov; // not in JP tree at all

  int run, trgId;
  unsigned int RunNumber, TriggerID; // incompatible type in JP tree

  float deltaT;

  float ctmqtqall, ctrmsts;
  double Qratio[2]; // incompatible type in JP tree
  double RMSTstart; // ditto

  float fido_qiv, fido_qid;

  double ctq0, ctq1, ctq; // trickery to match up with JP tree
  double ctqIV0, ctqIV1, ctqIV; // trickery to match up with JP tree

  int nidtubes, nivtubes;

  int ids_didfit;
  float ids_chi2, id_chi2;
  float ids_entr_x, ids_entr_y, ids_entr_z;
  float ids_end_x,  ids_end_y,  ids_end_z;

  float id_entr_x,  id_entr_y,  id_entr_z;
  float id_end_x,   id_end_y,   id_end_z;
  double Trk_MuHamID[2][3]; //incompatible semi-equivalent JP

  float ids_gclen, id_idexitqf, id_ivqbal;
  float ids_ivlen, ids_buflen;
  float ids_theta, ids_phi;
  float id_ivlen, id_buflen;

  bool id_didfit;

  double trgtime; // actually compatible -- only one!

  float ctX[3];
  double Vtx_BAMA[3]; // incompatible type in JP tree
  float ctEvisID;
  double EvisID; // incompatible type in JP tree
  float qrms, qdiff;
  double QRMS, Qdiff; // incompatible type in JP tree

  float pscs, psco; // only in JP
};

static const unsigned int noff = 92;
static const char * turnoff[noff] = {
"lilike",
"ctaplanarity",
"ctFlagMu",
"ctfwhm",
"ctgoodness",
"ctIDMuDeltaT",
"ctIVMuDeltaT",
"ctlightflux",
"ctmqtq",
"ctmqtqflag",
"ctnbadch",
"ctnbadchIV",
"ctnpe",
"ctnpeIV",
"ctnpulse",
"ctnpulseIV",
"ctphi",
"ctqtot",
"ctqtotIV",
"ctR",
"ctrho",
"ctsphericity",
"ctt2tot",
"cttmean",
"cttpeak",
"cttrise",
"ctXmuInGC",
"ctXmuInIV",
"ctXmuOuIV",
"date",
"fido_chi2",
"fido_gclen",
"fido_ivlen",
"fido_minuit_happiness",
"fido_phi",
"fido_stop",
"fido_targlen",
"fido_th",
"fido_used_ov",
"fido_didfit",
"fido_entrx",
"fido_entry",
"fido_entrz",
"fido_endx",
"fido_endy",
"fido_endz",
"hamphi",
"hamth",
"hamx",
"hamxe",
"HEMuDeltaT",
"IVX",
"nev",
"nhit",
"nhitIV",
"novhit",
"novloxy",
"novtrk",
"novupxy",
"ovbadtrk",
"ovloxylike",
"ovloxyx",
"ovloxyy",
"ovloxyz",
"ovtightloxy",
"ovtighttrk",
"ovtightupxy",
"ovtrigid",
"ovtrklike",
"ovtrkphi",
"ovtrkth",
"ovtrkx",
"ovtrky",
"ovupxylike",
"ovupxyx",
"ovupxyy",
"ovupxyz",
"pmtmultpe",
"pmtmultpe_IV",
"timeid",
"timeiv",
"tref",
"trefextIV",
"trefIV",
"trgWord",
"ttovtrig",
"vctnpulse",
"vctnpulseIV",
"vctq",
"vctqIV",
"vcttime",
"vcttimeIV"
};
