#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <sstream>
#include <iostream>

#include "TTree.h"
#include "TFile.h"
#include "TError.h"

// having dataparts and outdataparts is dumb.  I should merge them.
struct dataparts{
  bool coinov;
  int run;
  float ctmqtqall, ctrmsts;
  bool fido_stop;
  float fido_qiv, fido_qid;
  float fido_ivlen, fido_buflen, fido_gclen;
  float fido_chi2;
  float fido_endx, fido_endy, fido_endz;
  float fido_entrx, fido_entry, fido_entrz;
  int fido_nidtubes, fido_nivtubes;
  double deltaT;
  double trgtime;
  float ctX[3];
  float ctEvisID;
  float qrms, qdiff;
};

struct outdataparts{
  bool coinov;
  int run;
  int trgId;
  bool fido_stop;
  float fido_qiv, fido_qid;
  int fido_nidtubes, fido_nivtubes;
  float fido_chi2;
  float fido_ivlen, fido_buflen, fido_gclen;
  float fido_endx, fido_endy, fido_endz;
  float fido_entrx, fido_entry, fido_entrz;
  double deltaT;
  double trgtime;
  float ctX[3];
  float ctEvisID;
};

static const unsigned int noff = 85;
static const char * turnoff[noff] = { 
  "nev", "trgId", "trgWord", "date", "tref", "trefIV", "trefextIV",
  "nhitIV", "nhit", "ctq", "ctnpe", "ctqIV", "ctnpeIV", "ctnpulse",
  "ctnbadch", "ctrho", "ctphi", "ctR", "ctmqtq", "ctmqtqflag",
  "ctgoodness", "ctnbadchIV", "ctnpulseIV",
  "vctnpulseIV", "vctqIV", "vcttimeIV", "vctnpulse", "vctq",
  "vcttime", "pmtmultpe", "pmtmultpe_IV", "IVX", 
  "cttrise", "ctfwhm", "ctt2tot", "cttpeak", "cttmean",
  "ctIDMuDeltaT", "ctIVMuDeltaT", "HEMuDeltaT", "ctlightflux",
  "ctFlagMu", "ctXmuInGC", "ctXmuInIV", "ctXmuOuIV", "ctqtot",
  "ctqtotIV", "ctsphericity", "ctaplanarity", "lilike", 
  "timeid", "timeiv", "ovtrigid", "ttovtrig",
  "novhit", "novupxy", "novloxy", "novtrk", "ovupxyx",
  "ovupxyy", "ovupxyz", "ovtightupxy", "ovupxylike", "ovloxyx",
  "ovloxyy", "ovloxyz", "ovtightloxy", "ovloxylike", "ovtrkx",
  "ovtrky", "ovtrkth", "ovtrkphi", "ovtrklike", "ovbadtrk",
  "ovtighttrk", "hamx", "hamxe", "hamth", "hamphi", "fido_didfit", 
  "fido_used_ov", "fido_minuit_happiness",
  "fido_targlen", "fido_th", "fido_phi"
};

int main(int argc, char ** argv)
{
  int errcode = 0;

  if(argc < 3){
    fprintf(stderr, "I need at least two args\n");
    exit(1);
  }

  TFile * outfile = new TFile(argv[1], "recreate");
  TTree * const outtree = new TTree("data", "data");

  outdataparts outparts;

  #define BR(x) outtree->Branch(#x, &outparts.x);
  BR(coinov);
  BR(ctEvisID);
  BR(deltaT);
  BR(fido_buflen)
  BR(fido_chi2);
  BR(fido_endx);
  BR(fido_endy);
  BR(fido_endz);
  BR(fido_entrx);
  BR(fido_entry);
  BR(fido_entrz);
  BR(fido_gclen);
  BR(fido_ivlen);
  BR(fido_nidtubes);
  BR(fido_nivtubes);
  BR(fido_qid);
  BR(fido_qiv);
  BR(fido_stop);
  BR(run);
  BR(trgId);
  BR(trgtime);
  outtree->Branch("ctX[3]", outparts.ctX);

  for(int i = 2; i < argc; i++){
    printf("Reading %s\n", argv[i]);
    TFile * const infile = new TFile(argv[i], "read");
    outfile->cd();

    if(!infile || infile->IsZombie()){
      fprintf(stderr, "I couldn't read %s\n", argv[i]);
      errcode |= 1;
      continue;
    }

    TTree * const data = (TTree *)infile->Get("data");
    if(!data){
      fprintf(stderr, "%s doesn't have a data tree!\n", argv[i]);
      errcode |= 2;
      continue;
    }

    for(unsigned int j = 0; j < noff; j++)
      data->SetBranchStatus(turnoff[j], 0);

    dataparts parts;
    #define SBA(x) data->SetBranchAddress(#x, &parts.x);
    SBA(coinov);
    SBA(ctEvisID);
    SBA(ctmqtqall);
    SBA(ctrmsts);
    SBA(deltaT);
    SBA(fido_buflen);
    SBA(fido_chi2);
    SBA(fido_endx);
    SBA(fido_endy);
    SBA(fido_endz);
    SBA(fido_entrx);
    SBA(fido_entry);
    SBA(fido_entrz);
    SBA(fido_gclen);
    SBA(fido_ivlen);
    SBA(fido_nidtubes);
    SBA(fido_nivtubes);
    SBA(fido_qid);
    SBA(fido_qiv);
    SBA(fido_stop);
    SBA(qdiff);
    SBA(qrms);
    SBA(run);
    SBA(trgtime);
    data->SetBranchAddress("ctX", parts.ctX);

    TBranch * const ebranch = data->GetBranch("ctEvisID");

    TBranch * const sbranch = data->GetBranch("fido_stop");
    TBranch * const lnbr[4] = {
      data->GetBranch("ctrmsts"),
      data->GetBranch("qdiff"),
      data->GetBranch("qrms"),
      data->GetBranch("ctmqtqall") };

    for(int e = 0; e < data->GetEntries(); e++){
      ebranch->GetEntry(e);
      // Don't copy events with no energy in the ID, i.e.
      // failing the "valid trigger" cut.
      if(parts.ctEvisID < 0.4) continue;

      // don't copy light noise, but don't consider something light
      // noise if fido considers it a stopping muon.
      sbranch->GetEntry(e);
      if(!parts.fido_stop){
        lnbr[3]->GetEntry(e);
        if(parts.ctmqtqall > 0.12 || parts.ctmqtqall < 0) continue;

        lnbr[1]->GetEntry(e);
        if(parts.qdiff > 30e3) continue;

        lnbr[0]->GetEntry(e);
        lnbr[2]->GetEntry(e);
        if(parts.ctrmsts >= 36 && parts.qrms >= 464 - 8*parts.ctrmsts)
          continue;
      }

      data->GetEntry(e);

      #define TRANS(x) outparts.x = parts.x;
      TRANS(fido_stop);
      TRANS(fido_endx);
      TRANS(fido_endy);
      TRANS(fido_endz);
      TRANS(fido_entrx);
      TRANS(fido_entry);
      TRANS(fido_entrz);
      TRANS(fido_qiv);
      TRANS(fido_qid);
      TRANS(fido_chi2);
      TRANS(deltaT);
      TRANS(trgtime);
      TRANS(run);
      TRANS(coinov);
      TRANS(ctEvisID);
      TRANS(fido_nidtubes);
      TRANS(fido_nivtubes);
      TRANS(fido_ivlen);
      TRANS(fido_buflen);
      TRANS(fido_gclen);
      memcpy(outparts.ctX, parts.ctX, 3*sizeof(float));
      outparts.trgId = e;

      outtree->Fill();
    }
  
    delete data;
    delete infile;
  }

  outtree->Write();
  outfile->Close();
  return errcode;
}
