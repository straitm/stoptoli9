bool lightnoise(const float qrms, const float mqtq,
                const float rmsts, const float qdiff);

double fidocorrz(const double z);
double fidocorrx(const double x);
double fidocorry(const double y);

double bamacorrxy(const double xy, const double e);
double bamacorrz (const double  z, const double e);

struct dataparts{
  bool coinov;
  int run, trgId;
  float ctmqtqall, ctrmsts;
  float fido_qiv, fido_qid;
  double ctq, ctqIV;
  int nidtubes, nivtubes;

  int ids_didfit;
  float ids_chi2, id_chi2;
  float ids_end_x, ids_end_y, ids_end_z;
  float id_entr_x, id_entr_y;
  float ids_entr_x, ids_entr_y, ids_entr_z;
  float ids_gclen;
  float ids_ivlen, ids_buflen;
  float ids_theta, ids_phi;
  float id_ivlen, id_buflen;

  bool fido_didfit;
  float fido_entrx, fido_entry, fido_entrz;
  float fido_endx, fido_endy, fido_endz;

  double deltaT;
  double trgtime;
  float ctX[3];
  float ctEvisID;
  float qrms, qdiff;
};

static const unsigned int noff = 85;
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
