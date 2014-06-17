#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <sstream>
#include <iostream>

#include "TTree.h"
#include "TFile.h"
#include "TError.h"

struct dataparts{
  bool coinov;
  int run, trgId;
  float ctmqtqall, ctrmsts;
  bool fido_stop;
  float fido_qiv, fido_qid;
  float fido_chi2;
  float fido_endx, fido_endy, fido_endz;
  int fido_nidtubes, fido_nivtubes;
  float fido_ivlen, fido_buflen;
  double deltaT;
  double trgtime;
  float ctX[3];
  float ctEvisID;
  float qrms, qdiff;
};

static const unsigned int noff = 88;
static const char * turnoff[noff] = { 
  "nev", "trgWord", "date", "tref", "trefIV", "trefextIV",
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
  "fido_gclen", "fido_targlen", "fido_entrx",
  "fido_entry", "fido_entrz", "fido_th", "fido_phi"
};


static double maxtime = 1000;
static double minenergy = 4;
static double maxenergy = 14;

// True if this is light noise. If the variables aren't filled, i.e. all
// zeros, then this will return false.
static bool lightnoise(const float qrms, const float mqtq,
                       const float rmsts, const float qdiff)
{
  if(mqtq > 0.12 || mqtq < 0) return true;
  if(qdiff > 30e3) return true;
  if(rmsts >= 36 && qrms >= 464 - 8*rmsts) return true;
  return false;
}

static void searchfrommuon(dataparts & parts, TTree * const data,
                           const unsigned int muoni)
{
  const double mutime = parts.trgtime;
  const double mux = parts.fido_endx,
               muy = parts.fido_endy,
               muz = parts.fido_endz;
  const int mutrgid = parts.trgId;
  const int murun = parts.run;
  const bool mucoinov = parts.coinov;
  const float murchi2 =
    parts.fido_chi2/(parts.fido_nidtubes+parts.fido_nivtubes-6);
  const float muivdedx =
    parts.fido_qiv/(parts.fido_ivlen-parts.fido_buflen);

  unsigned int nneutron = 0;
  for(unsigned int i = 0; i < data->GetEntries(); i++){
    data->GetEntry(muoni+i);

    // Not more than ~4 nH lifetimes
    if(parts.trgtime - mutime > 800e3) break;

    if(parts.run != murun) break; // Stop at run boundaries

    // pass light noise
    if(lightnoise(parts.qrms, parts.ctmqtqall,
                  parts.ctrmsts, parts.qdiff)) continue;

    // right energy for a neutron capture
    if(!((parts.ctEvisID > 1.8 && parts.ctEvisID < 2.6) ||
         (parts.ctEvisID > 4.0 && parts.ctEvisID < 10 )))
      continue;
    
    // Do *not* require any particular position.

    nneutron++;
  }

  unsigned int got = 0;

  double b12x[2], b12y[2], b12z[2];

  for(unsigned int i = muoni+1; i < data->GetEntries(); i++){
    data->GetEntry(i);  

    if(parts.run != murun) break; // Stop at run boundaries

    const double b12time = parts.trgtime; 
    const double dt_ms = (b12time - mutime)/1e6;
    
    if(dt_ms < 1) continue; // Go past all H neutron captures

    if(dt_ms > maxtime) break; // Stop looking

    // Ignore low energy accidentals and H-neutrons
    if(parts.ctEvisID < minenergy) continue;

    // Ignore events above the end point + res
    if(parts.ctEvisID > maxenergy) continue;

    // No IV, OV energy
    if(parts.fido_qiv > 1000) continue;
    if(parts.coinov) continue;

    // pass light noise
    if(lightnoise(parts.qrms, parts.ctmqtqall,
                  parts.ctrmsts, parts.qdiff)) continue;

    b12x[got] = parts.ctX[0]*(0.970+0.013*parts.ctEvisID),
    b12y[got] = parts.ctX[1]*(0.970+0.013*parts.ctEvisID),
    b12z[got] = parts.ctX[2]*(0.985+0.005*parts.ctEvisID);

    const double dist = got == 0? sqrt(pow(mux-b12x[0], 2)
                                      +pow(muy-b12y[0], 2)
                                      +pow(muz-b12z[0], 2)):
                                  sqrt(pow(b12x[1]-b12x[0], 2)
                                      +pow(b12y[1]-b12y[0], 2)
                                      +pow(b12z[1]-b12z[0], 2));

    // And they must be near each other.
    if(dist > 800) continue;

    printf("B-12_for %d %d %d is %d dt %lf dist %f energy %f "
           "b12_at %f %f %f chi2,ivdedx: %f %f neutrons %d ",
           parts.run, mutrgid, mucoinov, parts.trgId, dt_ms, dist,
           parts.ctEvisID, b12x[got], b12y[got], b12z[got], murchi2,
           muivdedx, nneutron);
    got++;

    if(got >= 2) break;
  }
  if(got){
    printf("\n");
    fflush(stdout);
  }
}

static void search(dataparts & parts, TTree * const data)
{
  data->GetBranch("trgtime")->GetEntry(data->GetEntries()-1);
  const double eortime = parts.trgtime;

  double lastmuontime = 0; 

  for(unsigned int mi = 0; mi < data->GetEntries()-1; mi++){
    data->GetBranch("fido_stop")->GetEntry(mi);

    // Must be reasonably sure that it's a stopper
    if(!parts.fido_stop) goto end;

    data->GetBranch("trgtime")->GetEntry(mi);

    // Don't use a muon that is within our search window length from the
    // end of the run
    if(eortime - parts.trgtime < maxtime) break;

    // Require at least 500us since the last muon so we don't count
    // neutrons from the previous one as belonging to this one.
    if(parts.trgtime - lastmuontime < 500e3) goto end;

    data->GetEntry(mi);

    if(parts.fido_qiv < 5000) goto end;
    if(parts.fido_qid/8300 > 700) goto end;
    if(parts.fido_chi2/(parts.fido_nidtubes+parts.fido_nivtubes-6) > 10)
      goto end;

    if(parts.fido_nidtubes+parts.fido_nivtubes < 30) goto end;

    data->GetEntry(mi+1);

    // And must not have a Michel, loose cuts
    if(parts.deltaT < 5500 && parts.ctEvisID > 12) goto end;

    data->GetEntry(mi);
    searchfrommuon(parts, data, mi);

    end:

    data->GetBranch("coinov")->GetEntry(mi);
    data->GetBranch("fido_qiv")->GetEntry(mi);

    if(parts.coinov || parts.fido_qiv > 5000){
      data->GetBranch("trgtime")->GetEntry(mi);
      lastmuontime = parts.trgtime;
    }
  }
}

int main(int argc, char ** argv)
{
  if(argc < 5){
    fprintf(stderr, "b12search maxtime minenergy maxenergy files...\n");
    exit(1);
  }

  maxtime = atof(argv[1]);
  minenergy = atof(argv[2]);
  maxenergy = atof(argv[3]);

  gErrorIgnoreLevel = kFatal;
  int errcode = 0;
  for(int i = 4; i < argc; i++){
    TFile * const infile = new TFile(argv[i], "read");

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

    for(unsigned int i = 0; i < noff; i++)
      data->SetBranchStatus(turnoff[i], 0);

    dataparts parts;
    memset(&parts, 0, sizeof(parts));
    #define SBA(x) data->SetBranchAddress(#x, &parts.x);
    SBA(fido_stop);
    SBA(fido_endx);
    SBA(fido_endy);
    SBA(fido_endz);
    SBA(fido_qiv);
    SBA(fido_qid);
    SBA(fido_chi2);
    SBA(fido_nidtubes);
    SBA(fido_nivtubes);
    SBA(fido_ivlen);
    SBA(fido_buflen);
    SBA(deltaT);
    SBA(trgtime);
    SBA(run);
    SBA(trgId);
    SBA(ctmqtqall);
    SBA(ctrmsts);
    SBA(coinov);
    SBA(ctEvisID);
    SBA(qrms);
    SBA(qdiff);
    if(data->SetBranchAddress("ctX", parts.ctX) < 0)
      data->SetBranchAddress("ctX[3]", parts.ctX);

    search(parts, data);
  
    delete data;
    delete infile;
  }
  return errcode;
}
