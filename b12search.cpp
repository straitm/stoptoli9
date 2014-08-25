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
  float fido_qiv, fido_qid;
  int fido_nidtubes, fido_nivtubes;

  int ids_didfit;
  float ids_chi2;
  float ids_end_x, ids_end_y, ids_end_z;
  float ids_entr_x, ids_entr_y, ids_entr_z;
  float ids_gclen;
  float ids_ivlen, ids_buflen;

  double deltaT;
  double trgtime;
  float ctX[3];
  float ctEvisID;
  float qrms, qdiff;
};

static const unsigned int noff = 95;
static const char * turnoff[noff] = { 
  "fido_stop", "fido_chi2", "fido_endx", "fido_endy", "fido_endz",
  "fido_entrx", "fido_entry", "fido_entrz", "fido_gclen",
  "fido_ivlen", "fido_gclen",
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
  "fido_used_ov", "fido_minuit_happiness", "fido_targlen",
  "fido_th", "fido_phi"
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

static void searchfrommuon(dataparts & parts, TTree * const chtree,
                           TTree * const fitree,
                           const unsigned int muoni)
{
  const double mutime = parts.trgtime;
  const double entr_mux = parts.ids_entr_x,
               entr_muy = parts.ids_entr_y,
               entr_muz = parts.ids_entr_z;
  const double mux = parts.ids_end_x,
               muy = parts.ids_end_y,
               muz = parts.ids_end_z-170./3400.*(1700-parts.ids_end_z);
  const float gclen = parts.ids_gclen;
  const int mutrgid = parts.trgId;
  const int murun = parts.run;
  const bool mucoinov = parts.coinov;
  const float murchi2 =
    parts.ids_chi2/(parts.fido_nidtubes+parts.fido_nivtubes-6);
  const float muivdedx =
    parts.fido_qiv/(parts.ids_ivlen-parts.ids_buflen);

  unsigned int nneutronanydist = 0, ngdneutronanydist = 0;
  unsigned int nneutronnear = 0, ngdneutronnear = 0;

  double michelt = 0, michele = 0;

  for(unsigned int i = 1; i < chtree->GetEntries(); i++){
    chtree->GetEntry(muoni+i);
    fitree->GetEntry(muoni+i);

    if(parts.run != murun) break; // Stop at run boundaries

    const double dt = parts.trgtime - mutime;

    // For any uncut Michels or (hopefully) prompt gammas from muon
    // capture, record the time and energy. Note that if we are
    // operating on my microdsts, this will not pick up triggers with E
    // < 0.4MeV in the ID or light noise. I'm going to allow IV energy
    // since these may be very close to the muon event. Note that we 
    // will often count this Michel as a neutron also, so don't double
    // count by accident.
    if(dt < 5500 && !parts.coinov && michelt == 0 && michele == 0)
      michele = parts.ctEvisID, michelt = dt;

    // Skip past retriggers and whatnot, like DC3rdPub
    if(dt < 500) continue;

    // Not more than ~4 nH lifetimes
    if(dt > 800e3) break;

    // pass light noise
    if(lightnoise(parts.qrms, parts.ctmqtqall,
                  parts.ctrmsts, parts.qdiff)) continue;

    // right energy for a neutron capture
    if(!((parts.ctEvisID > 1.8 && parts.ctEvisID < 2.6) ||
         (parts.ctEvisID > 4.0 && parts.ctEvisID < 10 )))
      continue;
    
    const bool near = sqrt(pow(mux - parts.ctX[0], 2)+
                           pow(muy - parts.ctX[1], 2)+
                           pow(muz - parts.ctX[2], 2)) < 800; 

    const bool gd = parts.ctEvisID > 4.0 && parts.ctEvisID < 10
                 && parts.trgtime - mutime < 150e3;

    nneutronanydist++;
    if(gd) ngdneutronanydist++;
    if(near) nneutronnear++;
    if(gd && near) ngdneutronnear++;
  }

  unsigned int got = 0;

  // positions of putative isotope decays
  double ix[2], iy[2], iz[2];

  double lastmuontime = 0; 
  for(unsigned int i = muoni+1; i < chtree->GetEntries(); i++){
    chtree->GetEntry(i);  
    fitree->GetEntry(i);  

    if(parts.run != murun) break; // Stop at run boundaries

    const double itime = parts.trgtime; 
    const double dt_ms = (itime - mutime)/1e6;
    
    // Require at least 500us since the last muon so we don't count
    // neutrons as isotope decays
    if(itime - lastmuontime < 500e3) goto end;

    if(dt_ms < 1) goto end; // Go past all H neutron captures

    if(dt_ms > maxtime) break; // Stop looking

    // Ignore low energy accidentals and H-neutrons
    if(parts.ctEvisID < minenergy) goto end;

    // Ignore events above the end point + res
    if(parts.ctEvisID > maxenergy) goto end;

    // No IV, OV energy
    if(parts.fido_qiv > 1000) goto end;
    if(parts.coinov) goto end;

    // pass light noise
    if(lightnoise(parts.qrms, parts.ctmqtqall,
                  parts.ctrmsts, parts.qdiff)) goto end;

    ix[got] = parts.ctX[0]*(0.970+0.013*parts.ctEvisID),
    iy[got] = parts.ctX[1]*(0.970+0.013*parts.ctEvisID),
    iz[got] = parts.ctX[2]*(0.985+0.005*parts.ctEvisID);

    {
      // Record distance to previous selected event, not always to muon
      const double dist = got == 0?
        sqrt(pow(mux  -ix[0],2)+pow(muy  -iy[0],2)+pow(muz  -iz[0],2)):
        sqrt(pow(ix[0]-ix[1],2)+pow(iy[0]-iy[1],2)+pow(iz[0]-iz[1],2));

      // And they must be near each other.
      if(dist > 800) goto end;

      // run:itrig:coinov:mutrig:dt:dist:e:x:y:z:foo:chi2:ivdedx:
      // ngdnear:ngd:nnear:n:miche:micht
      printf("iso_for %d %d %d is %d dt %lf dist %f energy %f decay_at"
             " %f %f %f muon_at %f %f %f chi2,ivdedx: %f %f "
             "neutrons %d %d %d %d michel "
             "%lf %lf morefido %f %f %f %f ",
             murun, mutrgid, mucoinov, parts.trgId, dt_ms, dist,
             parts.ctEvisID, ix[got], iy[got], iz[got], mux, muy, muz,
             murchi2,
             muivdedx, ngdneutronnear, ngdneutronanydist,
                         nneutronnear,  nneutronanydist,
             michele, michelt, gclen, entr_mux, entr_muy, entr_muz);
      got++;

      if(got >= 2) break;
    }

    end:

    // Note the time of this event if it is a muon. Note that if using
    // my microdsts, only events with ID energy are in the input files,
    // so this only selects muons that cross the ID, which I think is
    // fine.
    chtree->GetBranch("coinov")->GetEntry(i);
    chtree->GetBranch("fido_qiv")->GetEntry(i);

    if(parts.coinov || parts.fido_qiv > 5000){
      chtree->GetBranch("trgtime")->GetEntry(i);
      lastmuontime = parts.trgtime;
    }

  }
  if(got){
    printf("\n");
    fflush(stdout);
  }
}

static double geteortime(dataparts & parts, TTree * const chtree,
                         const unsigned int start)
{
  TBranch * const runbr = chtree->GetBranch("run");
  TBranch * const timbr = chtree->GetBranch("trgtime");

  runbr->GetEntry(start);
  const int run = parts.run;
  double eortime = 0;

  for(int e = start; e < chtree->GetEntries(); e++){
    runbr->GetEntry(e);
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

static void search(dataparts & parts, TTree * const chtree,
                   TTree * const fitree)
{
  double lastmuontime = 0; 
  double eortime = 0;

  int run = 0;

  TBranch * const runbr = chtree->GetBranch("run");
  TBranch * const timbr = chtree->GetBranch("trgtime");
  TBranch * const stpbr = fitree->GetBranch("ids_didfit");
  
  for(unsigned int mi = 0; mi < chtree->GetEntries()-1; mi++){
    runbr->GetEntry(mi);

    if(parts.run != run){
      run = parts.run;
      lastmuontime = 0; // Assume a muon right before the run started.
      eortime = geteortime(parts, chtree, mi);
    }

    stpbr->GetEntry(mi);

    // Must be possible that it's a stopper
    if(!parts.ids_didfit) goto end;

    timbr->GetEntry(mi);

    // Don't use a muon that is within our search window length from the
    // end of the run
    if((eortime - parts.trgtime)/1e6 < maxtime) continue;

    // Require at least 500us since the last muon so we don't count
    // neutrons from the previous one as belonging to this one.
    if(parts.trgtime - lastmuontime < 500e3) goto end;

    chtree->GetEntry(mi);

    // XXX are these biasing me?
    //if(parts.fido_qiv < 5000) goto end;
    //if(parts.fido_qid/8300 > 700) goto end;
    //if(parts.fido_chi2/(parts.fido_nidtubes+parts.fido_nivtubes-6) > 10)
    //  goto end;

    if(parts.fido_nidtubes+parts.fido_nivtubes < 30) goto end;

    fitree->GetEntry(mi);
    searchfrommuon(parts, chtree, fitree, mi);

    end:

    // See note above for same code.
    chtree->GetBranch("coinov")->GetEntry(mi);
    chtree->GetBranch("fido_qiv")->GetEntry(mi);

    if(parts.coinov || parts.fido_qiv > 5000){
      chtree->GetBranch("trgtime")->GetEntry(mi);
      lastmuontime = parts.trgtime;
    }
  }
}

int main(int argc, char ** argv)
{
  if(argc < 6 && argc%2 != 0){
    fprintf(stderr, "b12search maxtime[ms] minenergy maxenergy "
                    "[cheetah file, fido file]*\n");
    exit(1);
  }

  maxtime = atof(argv[1]);
  minenergy = atof(argv[2]);
  maxenergy = atof(argv[3]);

  gErrorIgnoreLevel = kFatal;
  int errcode = 0;
  for(int i = 4; i < argc; i+=2){
    fputs(".", stderr);
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
    if(!chtree){
      fprintf(stderr, "\n%s lacks a \"data\" tree!\n", argv[i]);
      errcode |= 2;
      goto cleanup;
    }

    if(!fifile || fifile->IsZombie()){
      fprintf(stderr, "\nI couldn't read %s\n", argv[i+1]);
      errcode |= 4;
      goto cleanup;
    }

    fitree = (TTree *)fifile->Get("RecoMuonFIDOInfoTree");
    if(!fitree){
      fprintf(stderr, "\n%s lacks a RecoMuonFIDOInfoTree!\n", argv[i+1]);
      errcode |= 8;
      goto cleanup;
    }

    if(chtree->GetEntries() != fitree->GetEntries()){
      fprintf(stderr,
              "\n%s,\n%s:\ncheetah has %ld entries, but fido has %ld\n",
              argv[i], argv[i+1],
              long(chtree->GetEntries()), long(fitree->GetEntries()));
      errcode |= 0x10;
      goto cleanup;
    }

    fitree->SetMakeClass(1);

    for(unsigned int i = 0; i < noff; i++)
      chtree->SetBranchStatus(turnoff[i], 0);

    fitree->SetBranchStatus("*", 0);

    memset(&parts, 0, sizeof(parts));
    #define cSBA(x) chtree->SetBranchAddress(#x, &parts.x);
    #define fSBA(x) fitree->SetBranchStatus(#x, 1); \
                    fitree->SetBranchAddress(#x, &parts.x);
    fSBA(ids_didfit);
    fSBA(ids_end_x);
    fSBA(ids_end_y);
    fSBA(ids_end_z);
    fSBA(ids_entr_x);
    fSBA(ids_entr_y);
    fSBA(ids_entr_z);
    fSBA(ids_gclen);
    fSBA(ids_chi2);
    fSBA(ids_ivlen);
    fSBA(ids_buflen);

    cSBA(fido_nidtubes);
    cSBA(fido_nivtubes);
    cSBA(fido_qiv);
    cSBA(fido_qid);
    cSBA(deltaT);
    cSBA(trgtime);
    cSBA(run);
    cSBA(trgId);
    cSBA(ctmqtqall);
    cSBA(ctrmsts);
    cSBA(coinov);
    cSBA(ctEvisID);
    cSBA(qrms);
    cSBA(qdiff);
    if(chtree->SetBranchAddress("ctX", parts.ctX) < 0)
      chtree->SetBranchAddress("ctX[3]", parts.ctX);

    search(parts, chtree, fitree);

    cleanup:
  
    if(fitree) delete fitree;
    if(fifile) delete fifile;
    if(chtree) delete chtree;
    if(chfile) delete chfile;
  }
  return errcode;
}
