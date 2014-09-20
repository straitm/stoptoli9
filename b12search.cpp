#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <sstream>
#include <iostream>

#include "TTree.h"
#include "TFile.h"
#include "TError.h"

#include "search.h"

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
                           const unsigned int muoni,
                           const bool is_be12search,
                           const double timeleft)
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

  const float mufqid = parts.fido_qid,
              mufqiv = parts.fido_qiv,
              muctqid = parts.ctq,
              muctqiv = parts.ctqIV;

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
    * const ctEvisIDbr = chtree->GetBranch("ctEvisID"),
    * const trgtimebr  = chtree->GetBranch("trgtime");
    
  
  double deadtime = 0, nondeadenergy = 0;
  for(unsigned int i = muoni+1; i < chtree->GetEntries(); i++){
    trgtimebr->GetEntry(i);
    ctEvisIDbr->GetEntry(i);
    if(parts.trgtime-mutime < 6000) continue;
    deadtime = parts.trgtime-mutime;
    nondeadenergy = parts.ctEvisID;
    break;
  }


  for(unsigned int i = muoni+1; i < chtree->GetEntries(); i++){

#ifdef MULTIRUNFILES
    runbr->GetEntry(i);
    if(parts.run != murun) break; // Stop at run boundaries
#endif

    trgtimebr->GetEntry(i);
    coinovbr->GetEntry(i);
    const double dt = parts.trgtime - mutime;

    // For any uncut Michels or (hopefully) prompt gammas from muon
    // capture, record the time and energy. Note that if we are
    // operating on my microdsts, this will not pick up triggers with E
    // < 0.4MeV in the ID or light noise. I'm going to allow IV energy
    // since these may be very close to the muon event. Note that we
    // will often count this Michel as a neutron also for the zeroth
    // element of the count arrays (but not the first), so don't double
    // count by accident.  
    if(dt < 5500 && !parts.coinov && michelt == 0 && michele == 0){
      ctEvisIDbr->GetEntry(i);
      ctXbr->GetEntry(i);
      michele = parts.ctEvisID, michelt = dt;
      michdist = sqrt(pow(mux-parts.ctX[0], 2) +
                      pow(muy-parts.ctX[1], 2) +
                      pow(muz-parts.ctX[2], 2));
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
    if(lightnoise(parts.qrms, parts.ctmqtqall,
                  parts.ctrmsts, parts.qdiff)) continue;

    ctEvisIDbr->GetEntry(i);
    // right energy for a neutron capture
    if(!((parts.ctEvisID > 1.8 && parts.ctEvisID < 2.6) ||
         (parts.ctEvisID > 4.0 && parts.ctEvisID < 10 )))
      continue;

    ctXbr->GetEntry(i);
    const bool near = sqrt(pow(mux - parts.ctX[0], 2)+
                           pow(muy - parts.ctX[1], 2)+
                           pow(muz - parts.ctX[2], 2)) < 800;

    const bool gd = parts.ctEvisID > 4.0 && parts.ctEvisID < 10
                 && parts.trgtime - mutime < 150e3;

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

  double lastmuontime = 0;
  for(unsigned int i = muoni+1; i < chtree->GetEntries(); i++){

#ifdef MULTIRUNFILES
    runbr->GetEntry(i);
    if(parts.run != murun) break; // Stop at run boundaries
#endif

    trgtimebr->GetEntry(i);

    const double itime = parts.trgtime;
    const double dt_ms = (itime - mutime)/1e6;

    // Require at least 500us since the last muon so we don't count
    // neutrons as isotope decays
    if(itime - lastmuontime < 500e3) goto end;

    if(dt_ms < 1) goto end; // Go past all H neutron captures

    if(dt_ms > maxtime){ // stop looking
      if(printed == 0)
        printf("%d %d %d 0 0 0 0 "
               "0 0 0 %f %f %f %f %f "
               "%d %d %d %d %d %d %d %d "
               "%lf %.0lf %f %f %f %f %.0f %f %f "
               "%f %f %f %f %f\n",
               murun, mutrgid, mucoinov, mux, muy, muz,
               murchi2, muivdedx,
               ngdneutronnear[0], ngdneutronanydist[0],
               nneutronnear[0],  nneutronanydist[0],
               ngdneutronnear[1], ngdneutronanydist[1],
               nneutronnear[1],  nneutronanydist[1],
               michele, michelt, gclen, entr_mux, entr_muy, entr_muz,
               deadtime, nondeadenergy, michdist,
               mufqid, mufqiv, muctqid, muctqiv, timeleft);
      break;
    }

    ctEvisIDbr->GetEntry(i);

    // Ignore low energy accidentals and H-neutrons
    if(parts.ctEvisID < minenergy) goto end;

    // Ignore events above the end point + res
    if(parts.ctEvisID > maxenergy) goto end;
    
    fido_qivbr->GetEntry(i);
    // No IV, OV energy
    if(parts.fido_qiv > 1000) goto end;

    coinovbr->GetEntry(i);
    if(parts.coinov) goto end;

    qrmsbr->GetEntry(i);
    ctmqtqallbr->GetEntry(i);
    ctrmstsbr->GetEntry(i);
    qdiffbr->GetEntry(i);

    // pass light noise
    if(lightnoise(parts.qrms, parts.ctmqtqall,
                  parts.ctrmsts, parts.qdiff)) goto end;

    ctXbr->GetEntry(i);

    ix[got] = parts.ctX[0]*(0.970+0.013*parts.ctEvisID),
    iy[got] = parts.ctX[1]*(0.970+0.013*parts.ctEvisID),
    iz[got] = parts.ctX[2]*(0.985+0.005*parts.ctEvisID);

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

      printf("%d %d %d %d %lf %f %f "
             "%f %f %f %f %f %f %f %f "
             "%d %d %d %d %d %d %d %d "
             "%lf %.0lf %f %f %f %f %.0f %f %f "
             "%f %f %f %f %f%c",
             murun, mutrgid, mucoinov, parts.trgId, dt_ms, dist,
             parts.ctEvisID, ix[got], iy[got], iz[got], mux, muy, muz,
             murchi2, muivdedx,
             ngdneutronnear[0], ngdneutronanydist[0],
             nneutronnear[0],  nneutronanydist[0],
             ngdneutronnear[1], ngdneutronanydist[1],
             nneutronnear[1],  nneutronanydist[1],
             michele, michelt, gclen, entr_mux, entr_muy, entr_muz,
             deadtime, nondeadenergy, michdist,
             mufqid, mufqiv, muctqid, muctqiv,
             timeleft,
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

    if(parts.coinov || parts.fido_qiv > 5000)
      lastmuontime = parts.trgtime;

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

    if(pow(parts.ids_end_x, 2)+pow(parts.ids_end_y, 2) > pow(1708-35,2))
      goto end;

    ids_end_zbr->GetEntry(mi);

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

  printf("run:mutrig:ovcoin:trig:dt:dist:e:dx:dy:dz:mx:my:mz:chi2:"
         "ivdedx:ngdnear:ngd:nnear:n:latengdnear:latengd:latennear:"
         "laten:miche:micht:gclen:fex:fey:fez:deadt:deade:michd:fq:"
         "fqiv:cq:cqiv:timeleft\n");

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

    memset(&parts, 0, sizeof(parts));
    #define cSBA(x) chtree->SetBranchAddress(#x, &parts.x);
    #define fSBA(x) fitree->SetBranchStatus(#x, 1); \
                    fitree->SetBranchAddress(#x, &parts.x);
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
