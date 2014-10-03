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

#include "search.h"

static double maxtime = 1000;
static double minenergy = 4;
static double maxenergy = 14;

struct cart{
  double x, y, z;
};


struct track{
  float x0, y0, z0, x1, y1, z1;
  double tim;
};

static track maketrack(const float x0, const float y0, const float z0,
                       const float x1, const float y1, const float z1,
                       const float tim)
{
  track t;
  t.x0 = x0;
  t.y0 = y0;
  t.z0 = z0;
  t.x1 = x1;
  t.y1 = y1;
  t.z1 = z1;
  t.tim = tim;
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

static float lb12like(const vector<track> & ts, const float x,
                      const float y, const float z, const float dtim)
{
  float llike = -1000;
  for(unsigned int i = 0; i < ts.size(); i++){
    const float dist = ptol(ts[i], x, y, z);
    const double dt_ns = dtim - ts[i].tim;
    const float thisllike = (-dist/690.) /* econover thesis */
                          + (-dt_ns*log(2)/20.20e6);
    if(thisllike > llike) llike = thisllike;
  }
  return llike;
}

static void searchfrommuon(dataparts & bits, TTree * const chtree,
                           const unsigned int muoni,
                           const bool is_be12search,
                           const double timeleft)
{
  const double mutime = bits.trgtime;
  const double entr_mux = bits.ids_entr_x,
               entr_muy = bits.ids_entr_y,
               entr_muz = bits.ids_entr_z;
  const double mux = fidocorrxy(bits.ids_end_x),
               muy = fidocorrxy(bits.ids_end_y),
               muz = fidocorrz(bits.ids_end_z);
  const float gclen = bits.ids_gclen;
  const int mutrgid = bits.trgId;
  const int murun = bits.run;
  const bool mucoinov = bits.coinov;
  const float murchi2 =
    bits.ids_chi2/(bits.fido_nidtubes+bits.fido_nivtubes-6);
  const float muivdedx =
    bits.fido_qiv/(bits.ids_ivlen-bits.ids_buflen);

  const float mufqid = bits.fido_qid,
              mufqiv = bits.fido_qiv,
              muctqid = bits.ctq,
              muctqiv = bits.ctqIV;

  unsigned int nneutronanydist[2] = {0,0}, ngdneutronanydist[2] = {0,0};
  unsigned int nneutronnear[2] = {0,0}, ngdneutronnear[2] = {0,0};

  double michelt = 0, michele = 0, michdist = 0;

  TBranch 
    * const runbr      = chtree->GetBranch("run"),
    * const coinovbr   = chtree->GetBranch("coinov"),
    * const trgIdbr    = chtree->GetBranch("trgId"),
    * const ctXbr      = chtree->GetBranch("ctX"),
    * const qrmsbr     = chtree->GetBranch("qrms"),
    * const ctmqtqallbr= chtree->GetBranch("ctmqtqall"),
    * const ctrmstsbr  = chtree->GetBranch("ctrmsts"),
    * const qdiffbr    = chtree->GetBranch("qdiff"),
    * const fido_qivbr = chtree->GetBranch("fido_qiv"),
    * const fido_didfitbr=chtree->GetBranch("fido_didfit"),
    * const fido_entrxbr=chtree->GetBranch("fido_entrx"),
    * const fido_entrybr=chtree->GetBranch("fido_entry"),
    * const fido_entrzbr=chtree->GetBranch("fido_entrz"),
    * const fido_endxbr= chtree->GetBranch("fido_endx"),
    * const fido_endybr= chtree->GetBranch("fido_endy"),
    * const fido_endzbr= chtree->GetBranch("fido_endz"),
    * const ctEvisIDbr = chtree->GetBranch("ctEvisID"),
    * const trgtimebr  = chtree->GetBranch("trgtime");
  
  double deadtime = 0, nondeadenergy = 0;
  for(unsigned int i = muoni+1; i < chtree->GetEntries(); i++){
    trgtimebr->GetEntry(i);
    ctEvisIDbr->GetEntry(i);
    if(bits.trgtime-mutime < 6000) continue;
    deadtime = bits.trgtime-mutime;
    nondeadenergy = bits.ctEvisID;
    break;
  }


  for(unsigned int i = muoni+1; i < chtree->GetEntries(); i++){

#ifdef MULTIRUNFILES
    runbr->GetEntry(i);
    if(bits.run != murun) break; // Stop at run boundaries
#endif

    trgtimebr->GetEntry(i);
    coinovbr->GetEntry(i);
    const double dt = bits.trgtime - mutime;

    // For any uncut Michels or (hopefully) prompt gammas from muon
    // capture, record the time and energy. Note that if we are
    // operating on my microdsts, this will not pick up triggers with E
    // < 0.4MeV in the ID or light noise. I'm going to allow IV energy
    // since these may be very close to the muon event. Note that we
    // will often count this Michel as a neutron also for the zeroth
    // element of the count arrays (but not the first), so don't double
    // count by accident.  
    if(dt < 5500 && !bits.coinov && michelt == 0 && michele == 0){
      ctEvisIDbr->GetEntry(i);
      ctXbr->GetEntry(i);
      michele = bits.ctEvisID, michelt = dt;
      michdist =
        sqrt(pow(mux-bamacorrxy(bits.ctX[0], bits.ctEvisID), 2) +
             pow(muy-bamacorrxy(bits.ctX[1], bits.ctEvisID), 2) +
             pow(muz-bamacorrz( bits.ctX[2], bits.ctEvisID), 2));
    }

    // Skip past retriggers and whatnot, like DC3rdPub
    if(dt < 500) continue;

    // Not more than ~4 nH lifetimes
    if(dt > 800e3) break;

    qrmsbr->GetEntry(i);
    ctmqtqallbr->GetEntry(i);
    ctrmstsbr->GetEntry(i);
    qdiffbr->GetEntry(i);

    // pass light noise
    if(lightnoise(bits.qrms, bits.ctmqtqall,
                  bits.ctrmsts, bits.qdiff)) continue;

    ctEvisIDbr->GetEntry(i);
    // right energy for a neutron capture
    if(!((bits.ctEvisID > 1.8 && bits.ctEvisID < 2.6) ||
         (bits.ctEvisID > 4.0 && bits.ctEvisID < 10 )))
      continue;

    ctXbr->GetEntry(i);
    const bool near =
      sqrt(pow(mux - bamacorrxy(bits.ctX[0], bits.ctEvisID), 2)+
           pow(muy - bamacorrxy(bits.ctX[1], bits.ctEvisID), 2)+
           pow(muz - bamacorrz( bits.ctX[2], bits.ctEvisID), 2)) < 800;

    const bool gd = bits.ctEvisID > 4.0 && bits.ctEvisID < 10
                 && bits.trgtime - mutime < 150e3;

    const bool alsoamichel = dt < 5500;
    nneutronanydist[0]++;
    if(gd) ngdneutronanydist[0]++;
    if(near) nneutronnear[0]++;
    if(gd && near) ngdneutronnear[0]++;
    if(!alsoamichel){
      nneutronanydist[1]++;
      if(gd) ngdneutronanydist[1]++;
      if(near) nneutronnear[1]++;
      if(gd && near) ngdneutronnear[1]++;
    }
  }

  unsigned int got = 0, printed = 0;

  // positions of putative isotope decays
  double ix[2], iy[2], iz[2];

  double lastmuontime = mutime, lastgcmuontime = mutime,
         lastvalidtime = mutime;

  vector<track> tmuons;

  for(unsigned int i = muoni+1; i < chtree->GetEntries(); i++){

#ifdef MULTIRUNFILES
    runbr->GetEntry(i);
    if(bits.run != murun) break; // Stop at run boundaries
#endif

    trgtimebr->GetEntry(i);

    const double itime = bits.trgtime;
    const double dt_ms = (itime - mutime)/1e6;
    const double ttlastvalid =(itime-lastvalidtime  )/1e6;
    const double ttlastmuon  =(itime-lastmuontime   )/1e6;
    const double ttlastgcmuon=(itime-lastgcmuontime )/1e6;

    // Require at least 500us since the last muon so we don't count
    // neutrons as isotope decays
    if(ttlastmuon < 0.5) goto end;

    if(dt_ms < 1) goto end; // Go past all H neutron captures

    if(dt_ms > maxtime){ // stop looking
      if(printed == 0)
        // NOTE-luckplan
        printf("0 0 0 0 0 0 0 0 "
               #define LATEFORM \
               "%d %d %d " \
               "%f %f %f %f %f " \
               "%d %d %d %d %d %d %d %d " \
               "%lf %.0lf %f %f %f %f %.0f %f %f " \
               "%f %f %f %f %f %f %f %f"
               LATEFORM
               "\n",
               /* 0, 0, 0, 0, 0, 0, 0, */
               #define LATEVARS \
               murun, mutrgid, mucoinov, \
               mux, muy, muz, murchi2, muivdedx, \
               ngdneutronnear[0], ngdneutronanydist[0], \
               nneutronnear[0],  nneutronanydist[0], \
               ngdneutronnear[1], ngdneutronanydist[1], \
               nneutronnear[1],  nneutronanydist[1], \
               michele, michelt, gclen, entr_mux, entr_muy, entr_muz, \
               deadtime, nondeadenergy, michdist, \
               mufqid, mufqiv, muctqid, muctqiv, \
               timeleft, ttlastvalid, ttlastmuon, \
               ttlastgcmuon
              LATEVARS);
      break;
    }

    ctEvisIDbr->GetEntry(i);

    // Ignore low energy accidentals and H-neutrons
    if(bits.ctEvisID < minenergy) goto end;

    // Ignore events above the end point + res
    if(bits.ctEvisID > maxenergy) goto end;
    
    fido_qivbr->GetEntry(i);
    // No IV, OV energy
    if(bits.fido_qiv > 1000) goto end;

    coinovbr->GetEntry(i);
    if(bits.coinov) goto end;

    qrmsbr->GetEntry(i);
    ctmqtqallbr->GetEntry(i);
    ctrmstsbr->GetEntry(i);
    qdiffbr->GetEntry(i);

    // pass light noise
    if(lightnoise(bits.qrms, bits.ctmqtqall, bits.ctrmsts, bits.qdiff))
      goto end;

    ctXbr->GetEntry(i);

    ix[got] = bamacorrxy(bits.ctX[0], bits.ctEvisID),
    iy[got] = bamacorrxy(bits.ctX[1], bits.ctEvisID),
    iz[got] = bamacorrz( bits.ctX[2], bits.ctEvisID);

    {
      // If the be12 search, record distance to previous selected event,
      // not always to muon, otherwise always to the muon.
      const double dist = got == 0 || !is_be12search?
        sqrt(pow(mux  -ix[0],2)+pow(muy  -iy[0],2)+pow(muz  -iz[0],2)):
        sqrt(pow(ix[0]-ix[1],2)+pow(iy[0]-iy[1],2)+pow(iz[0]-iz[1],2));

      // And they must be near each other.
      if(dist > 800) goto end;

      trgIdbr->GetEntry(i);
      runbr->GetEntry(i);

      // NOTE-luckplan
      printf("%d %lf %f %f %f %f %f %f " LATEFORM "%c",
             bits.trgId, dt_ms, dist,
             bits.ctEvisID, ix[got], iy[got], iz[got],
             lb12like(tmuons, ix[got], iy[got], iz[got], bits.trgtime),
             LATEVARS,
             is_be12search?' ':'\n');
      printed++;
      
      // If searching for be12, use the ix/iy/iz arrays and stop when we
      // get 2 decays; these are printed on the same line. Otherwise,
      // keep going up to the time limit and put each on its own line.
      if(is_be12search) got++;

      if(got >= 2) break;
    }

    end:

    // Note the time of this event if it is a muon. Note that if using
    // my microdsts, only events with ID energy are in the input files,
    // so this only selects muons that cross the ID, which I think is
    // fine.
    coinovbr->GetEntry(i);
    fido_qivbr->GetEntry(i);
    fido_didfitbr->GetEntry(i);
    ctEvisIDbr->GetEntry(i);

    if(bits.coinov || bits.fido_qiv > 5000 || bits.ctEvisID > 60)
      lastmuontime = bits.trgtime;

    if(bits.fido_qiv > 5000 && bits.ctEvisID > 60)
      lastgcmuontime = bits.trgtime;

    if(bits.fido_didfit){
      fido_entrxbr->GetEntry(i);
      fido_entrybr->GetEntry(i);
      fido_entrzbr->GetEntry(i);
      fido_endxbr ->GetEntry(i);
      fido_endybr ->GetEntry(i);
      fido_endzbr ->GetEntry(i);

      tmuons.push_back(maketrack(
        bits.fido_entrx, bits.fido_entry, bits.fido_entrz,
        bits.fido_endx, bits.fido_endy, bits.fido_endz,
        bits.trgtime));
    }

    // Note the time of this even if it is valid, which for me means
    // not light noise and at least 0.4 MeV
    qrmsbr->GetEntry(i);
    ctmqtqallbr->GetEntry(i);
    ctrmstsbr->GetEntry(i);
    qdiffbr->GetEntry(i);
    ctEvisIDbr->GetEntry(i);

    if(!lightnoise(bits.qrms, bits.ctmqtqall, bits.ctrmsts, bits.qdiff)
       && bits.ctEvisID > 0.4)
      lastvalidtime = bits.trgtime;
  }
  if(got){
    printf("\n");
    fflush(stdout);
  }
}

/* Gets the time of the end of the run */
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
                   TTree * const fitree, const bool is_be12search)
{
  double lastmuontime = 0;
  double eortime = 0;
  int run = 0;

  TBranch * const runbr           = chtree->GetBranch("run");
  TBranch * const trgtimebr       = chtree->GetBranch("trgtime");
  TBranch * const coinovbr        = chtree->GetBranch("coinov");
  TBranch * const fido_qivbr      = chtree->GetBranch("fido_qiv");
  TBranch * const fido_qidbr      = chtree->GetBranch("fido_qid");
  TBranch * const fido_nidtubesbr = chtree->GetBranch("fido_nidtubes");
  TBranch * const fido_nivtubesbr = chtree->GetBranch("fido_nivtubes");

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

  for(unsigned int mi = 0; mi < chtree->GetEntries()-1; mi++){
    runbr->GetEntry(mi);

    if(parts.run != run){
      run = parts.run;
      lastmuontime = 0; // Assume a muon right before the run started.
      eortime = geteortime(parts, chtree, mi);
    }

    ids_didfitbr->GetEntry(mi);

    // Must be possible that it's a stopper
    if(!parts.ids_didfit) goto end;

    trgtimebr->GetEntry(mi);

    // Require at least 500us since the last muon so we don't count
    // neutrons from the previous one as belonging to this one.
    if(parts.trgtime - lastmuontime < 500e3) goto end;

    fido_qivbr->GetEntry(mi);

    if(parts.fido_qiv < 5000) goto end;

    fido_qidbr->GetEntry(mi);

    if(parts.fido_qid/8300 > 700) goto end;

    ids_chi2br->GetEntry(mi);
    fido_nivtubesbr->GetEntry(mi);
    fido_nidtubesbr->GetEntry(mi);

    if(parts.ids_chi2/(parts.fido_nidtubes+parts.fido_nivtubes-6) > 10)
      goto end;

    ids_end_xbr->GetEntry(mi);
    ids_end_ybr->GetEntry(mi);

    // This correctly uses the *uncorrected* position
    if(pow(parts.ids_end_x, 2)+pow(parts.ids_end_y, 2) > pow(1708-35,2))
      goto end;

    ids_end_zbr->GetEntry(mi);

    // This correctly uses the *uncorrected* position
    if(parts.ids_end_z < -1786+35) goto end;

    ids_entr_zbr->GetEntry(mi);
    id_ivlenbr->GetEntry(mi);
    id_buflenbr->GetEntry(mi);

    if(parts.ids_entr_z > 11500 -
       62*parts.fido_qiv/(parts.id_ivlen-parts.id_buflen)) goto end;

    id_chi2br->GetEntry(mi);
    ids_chi2br->GetEntry(mi);

    if(parts.ids_chi2-parts.id_chi2 > 800) goto end;

    id_entr_xbr->GetEntry(mi);
    id_entr_ybr->GetEntry(mi);
    ids_entr_xbr->GetEntry(mi);
    ids_entr_ybr->GetEntry(mi);

    if(pow(parts.id_entr_x, 2)+pow(parts.id_entr_y, 2) < pow(1000, 2) &&
       pow(parts.ids_entr_x,2)+pow(parts.ids_entr_y,2) > pow(2758, 2))
      goto end;

    if(parts.fido_nidtubes+parts.fido_nivtubes < 6) goto end;
 
    chtree->GetEntry(mi);
    fitree->GetEntry(mi);

    searchfrommuon(parts, chtree, mi, is_be12search, 
    // Record how much time is left to the end of the run so that we can
    // exclude events that are too close when analyzing.
                             (eortime - parts.trgtime)/1e6);
    end:

    // See note above for same code.
    coinovbr->GetEntry(mi);
    fido_qivbr->GetEntry(mi);

    if(parts.coinov || parts.fido_qiv > 5000){
      trgtimebr->GetEntry(mi);
      lastmuontime = parts.trgtime;
    }
  }
}

int main(int argc, char ** argv)
{
  if(argc < 6 && argc%2 != 0){
    fprintf(stderr,
            "b12search maxtime[ms] minenergy maxenergy "
            "[cheetah file, fido file]*\n\n");

    fprintf(stderr,
            "To run the be12 search, looking for exactly two decays:\n"
            "be12search maxtime[ms] minenergy maxenergy "
            "[cheetah file, fido file]*\n\n");

    fprintf(stderr,
            "To check files only, use\n"
            "checkb12search foo foo foo "
            "[cheetah file, fido file]*\n");

    exit(1);
  }

  maxtime = atof(argv[1]);
  minenergy = atof(argv[2]);
  maxenergy = atof(argv[3]);

  const bool is_be12search = !strcmp(basename(argv[0]), "be12search");

  gErrorIgnoreLevel = kFatal;
  int errcode = 0;

  fprintf(stderr, "Processing %d runs\n", (argc-4)/2);

  // NOTE-luckplan
  printf("trig/I:dt/F:dist/F:e/F:dx/F:dy/F:dz/F:b12like/F:"
         "run/I:mutrig/I:ovcoin/I:mx/F:my/F:mz/F:"
         "chi2/F:ivdedx/F:ngdnear/I:ngd/I:nnear/I:n/I:latengdnear/I:"
         "latengd/I:latennear/I:laten/I:miche/F:micht/F:gclen/F:"
         "fex/F:fey/F:fez/F:deadt/F:"
         "deade/F:michd/F:fq/F:fqiv/F:cq/F:cqiv/F:timeleft/F:"
         "ttlastvalid/F:ttlastmuon/F:ttlastgcmuon/F"
         "\n");

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
      fprintf(stderr,"\n%s lacks a RecoMuonFIDOInfoTree!\n",argv[i+1]);
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
    #define cSBA(x) chtree->SetBranchAddress(#x, &parts.x);
    #define fSBA(x) fitree->SetBranchStatus(#x, 1); \
                    fitree->SetBranchAddress(#x, &parts.x); \
                    fitree->AddBranchToCache(#x);
 
    fSBA(ids_didfit);
    fSBA(ids_end_x);
    fSBA(ids_end_y);
    fSBA(ids_end_z);
    fSBA(id_entr_x);
    fSBA(id_entr_y);
    fSBA(ids_entr_x);
    fSBA(ids_entr_y);
    fSBA(ids_entr_z);
    fSBA(ids_gclen);
    fSBA(ids_chi2);
    fSBA(ids_ivlen);
    fSBA(ids_buflen);
    fSBA(id_chi2);
    fSBA(id_ivlen);
    fSBA(id_buflen);

    cSBA(fido_nidtubes);
    cSBA(fido_nivtubes);
    cSBA(fido_qiv);
    cSBA(fido_qid);
    cSBA(fido_didfit);
    cSBA(fido_entrx);
    cSBA(fido_entry);
    cSBA(fido_entrz);
    cSBA(fido_endx);
    cSBA(fido_endy);
    cSBA(fido_endz);
    cSBA(ctq);
    cSBA(ctqIV);
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

    search(parts, chtree, fitree, is_be12search);

    cleanup:

    if(fitree) delete fitree;
    if(fifile) delete fifile;
    if(chtree) delete chtree;
    if(chfile) delete chfile;
  }
  return errcode;
}
