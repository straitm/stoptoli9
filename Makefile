all: stoptoli9 b12search

# No code depends on these results.  They are just to illustrate examples 
# in the technote.
technote: \
  distcuteff_byenergy_finalfit.out.technote \
  distcuteff_topbottom_finalfit.out.technote

analysis: \
  technote \
  b12groundstate_finalfit.out \
  c9_finalfit.out \
  be12_finalfit.out \
  b8_finalfit.out \
  n12_finalfit.out \
  b14_finalfit.out \
  n16_finalfit.out \
  he6_finalfit_0.out \
  he6_finalfit_1.out \
  he6_finalfit_2.out \
  he6_finalfit_3.out \
  li9_finalfit_0.out \
  li9_finalfit_1.out

li9_finalfit_C.so: \
  consts.h \
  distcuteff_targetgc_finalfit.out.h \
  totallivetime_finalfit.out.h \
  noncarbondenominators_finalfit.out.h \
  carbondenominators_finalfit.out.h \
  reactorpowerbin.C \
  li9_finalfit.C
	root -n -l -b -q li9_finalfit.C++'(-2)'

li9_finalfit_-1.out.h: \
  li9_finalfit_C.so \
  run_analysis_li9.sh
	./run_analysis_li9.sh -1

li9_finalfit_0.out: \
  li9_finalfit_C.so \
  run_analysis_li9.sh
	./run_analysis_li9.sh 0

li9_finalfit_1.out: \
  li9_finalfit_C.so \
  li9_finalfit_-1.out.h \
  run_analysis_li9.sh
	./run_analysis_li9.sh 1

n16_finalfit.out: \
  consts.h \
  sub_muon_eff.out.h \
  distcuteff_wholeloose_finalfit.out.h \
  distcuteff_b12like_finalfit.out.h \
  totallivetime_finalfit.out.h \
  noncarbondenominators_finalfit.out.h \
  carbondenominators_finalfit.out.h \
  run_analysis_n16.sh \
  n16_finalfit.C
	./run_analysis_n16.sh

b14_finalfit.out: \
  consts.h \
  sub_muon_eff.out.h \
  distcuteff_b12like_finalfit.out.h \
  distcuteff_wholeloose_finalfit.out.h \
  totallivetime_finalfit.out.h \
  noncarbondenominators_finalfit.out.h \
  run_analysis_b14.sh \
  b14_finalfit.C
	./run_analysis_b14.sh

n12_finalfit.out: \
  consts.h \
  sub_muon_eff.out.h \
  n12cutefficiency_finalfit.out.h \
  distcuteff_wholeloose_finalfit.out.h \
  totallivetime_finalfit.out.h \
  noncarbondenominators_finalfit.out.h \
  run_analysis_n12.sh \
  n12_finalfit.C
	./run_analysis_n12.sh

c9_finalfit.out: \
  consts.h \
  sub_muon_eff.out.h \
  distcuteff_wholeloose_finalfit.out.h \
  totallivetime_finalfit.out.h \
  noncarbondenominators_finalfit.out.h \
  run_analysis_c9.sh \
  c9_finalfit.C
	./run_analysis_c9.sh

he6_finalfit_C.so: \
  he6_finalfit.C \
  sub_muon_eff.out.h \
  consts.h \
  distcuteff_b12like_finalfit.out.h \
  distcuteff_he6_finalfit.out.h \
  totallivetime_finalfit.out.h \
  carbondenominators_finalfit.out.h
	root -n -l -b -q he6_finalfit.C++'(-1)'

he6_finalfit_0.out: \
  he6_finalfit_C.so \
  run_analysis_he6.sh
	./run_analysis_he6.sh 0

he6_finalfit_1.out: \
  he6_finalfit_C.so \
  run_analysis_he6.sh
	./run_analysis_he6.sh 1

he6_finalfit_2.out: \
  he6_finalfit_C.so \
  run_analysis_he6.sh
	./run_analysis_he6.sh 2

he6_finalfit_3.out: \
  he6_finalfit_C.so \
  run_analysis_he6.sh
	./run_analysis_he6.sh 3

b8_finalfit.out: \
  consts.h \
  sub_muon_eff.out.h \
  b8cutefficiency.out.h \
  distcuteff_wholeloose_finalfit.out.h \
  totallivetime_finalfit.out.h \
  carbondenominators_finalfit.out.h \
  run_analysis_b8.sh \
  b8_finalfit.C
	./run_analysis_b8.sh

be12_finalfit.out: \
  consts.h \
  sub_muon_eff.out.h \
  distcuteff_wholeloose_finalfit.out.h \
  b12cutefficiency_finalfit.out.h \
  totallivetime_finalfit.out.h \
  carbondenominators_finalfit.out.h \
  run_analysis_be12.sh \
  be12_finalfit.C
	./run_analysis_be12.sh

li8_finalfit.out.h: \
  consts.h \
  sub_muon_eff.out.h \
  distcuteff_b12like_finalfit.out.h \
  distcuteff_wholeloose_finalfit.out.h \
  carbondenominators_finalfit.out.h \
  totallivetime_finalfit.out.h \
  li8cutefficiency_finalfit.out.h \
  b12cutefficiency_finalfit.out.h \
  run_analysis_li8.sh \
  li8_finalfit.C
	./run_analysis_li8.sh

b12groundstate_finalfit.out: \
  consts.h \
  fullb12_finalfit.out.h \
  totallivetime_finalfit.out.h \
  b12cutefficiency_finalfit.out.h \
  carbondenominators_finalfit.out.h \
  b12gamma_finalfit_0.out.h \
  b12gamma_finalfit_1.out.h \
  run_analysis_b12groundstate.sh \
  b12groundstate_finalfit.C
	./run_analysis_b12groundstate.sh

b12gamma_finalfit_C.so: \
  b12gamma_finalfit.C \
  consts.h \
  sub_muon_eff.out.h \
  distcuteff_wholeloose_finalfit.out.h \
  li8cutefficiency_finalfit.out.h \
  b12cutefficiency_finalfit.out.h \
  totallivetime_finalfit.out.h \
  li8_finalfit.out.h \
  fullb12_finalfit.out.h \
  carbondenominators_finalfit.out.h
	root -n -l -b -q b12gamma_finalfit.C++'(-1)'

b12gamma_finalfit_0.out.h: \
  b12gamma_finalfit_C.so \
  run_analysis_b12gamma.sh
	./run_analysis_b12gamma.sh 0

b12gamma_finalfit_1.out.h: \
  b12gamma_finalfit_C.so \
  run_analysis_b12gamma.sh
	./run_analysis_b12gamma.sh 1

fullb12_finalfit.out.h: \
  consts.h \
  sub_muon_eff.out.h \
  b12cutefficiency_finalfit.out.h \
  carbondenominators_finalfit.out.h \
  run_analysis_fullb12.sh \
  fullb12_finalfit.C
	./run_analysis_fullb12.sh

n12cutefficiency_finalfit.out.h: \
  n12cutefficiency_finalfit.C \
	./run_analysis_n12cutefficiency.sh

b8cutefficiency_finalfit.out.h: \
  b8cutefficiency_finalfit.C \
	./run_analysis_b8cutefficiency.sh

li8cutefficiency_finalfit.out.h: \
  li8cutefficiency_finalfit.C \
	./run_analysis_li8cutefficiency.sh

2cutefficiency_finalfit.out.h: \
  consts.h \
  b12cutefficiency_finalfit.C \
  b12spectrum.C
	./run_analysis_b12cutefficiency.sh

distcuteff_byenergy_finalfit.out.technote: \
  distcuteff_byenergy_finalfit.C run_analysis_distcuteff.sh \
  distcuteff_include.C b12like_finalfit.C reactorpowerbin.C \
  b12cutefficiency_finalfit.out.h totallivetime_finalfit.out.h
	./run_analysis_distcuteff.sh byenergy

distcuteff_b12like_finalfit.out.h: \
  distcuteff_b12like_finalfit.C run_analysis_distcuteff.sh \
  distcuteff_include.C b12like_finalfit.C reactorpowerbin.C \
  b12cutefficiency_finalfit.out.h totallivetime_finalfit.out.h
	./run_analysis_distcuteff.sh b12like

distcuteff_he6_finalfit.out.h: \
  distcuteff_he6_finalfit.C run_analysis_distcuteff.sh \
  distcuteff_include.C b12like_finalfit.C reactorpowerbin.C \
  b12cutefficiency_finalfit.out.h totallivetime_finalfit.out.h
	./run_analysis_distcuteff.sh he6

distcuteff_targetgc_finalfit.out.h: \
  distcuteff_targetgc_finalfit.C run_analysis_distcuteff.sh \
  distcuteff_include.C b12like_finalfit.C reactorpowerbin.C \
  b12cutefficiency_finalfit.out.h totallivetime_finalfit.out.h
	./run_analysis_distcuteff.sh targetgc

distcuteff_topbottom_finalfit.out.technote: \
  distcuteff_topbottom_finalfit.C run_analysis_distcuteff.sh \
  distcuteff_include.C b12like_finalfit.C reactorpowerbin.C \
  b12cutefficiency_finalfit.out.h totallivetime_finalfit.out.h
	./run_analysis_distcuteff.sh topbottom

distcuteff_wholeloose_finalfit.out.h: \
  distcuteff_wholeloose_finalfit.C run_analysis_distcuteff.sh \
  distcuteff_include.C b12like_finalfit.C reactorpowerbin.C \
  b12cutefficiency_finalfit.out.h totallivetime_finalfit.out.h
	./run_analysis_distcuteff.sh wholeloose

sub_muon_eff.out.h: \
  getmuondeadtime.sh \
  getliveanddeadtime.sh \
  run_analysis_sub_muon_eff.sh
	./run_analysis_sub_muon_eff.sh

totallivetime_finalfit.out.h: \
  reactorpowerbin.C \
  totallivetime_finalfit.C \
  run_analysis_totallivetime.sh \
  runlist
	./run_analysis_totallivetime.sh

musicchargeratio_finalfit.out.h: \
  musicchargeratio_finalfit.C
	./run_analysis_musicchargeratio.sh

carbondenominators_finalfit.out.h: \
  consts.h \
  sub_muon_eff.out.h \
  totallivetime_finalfit.out.h \
  li8cutefficiency_finalfit.out.h \
  b12cutefficiency_finalfit.out.h \
  run_analysis_carbondenominators.sh \
  musicchargeratio_finalfit.out.h \
  mucount_finalfit.C \
  carbondenominators_finalfit.C
	./run_analysis_carbondenominators.sh

noncarbondenominators_finalfit.out.h: \
  consts.h \
  carbondenominators_finalfit.out.h \
  run_analysis_noncarbondenominators.sh \
  noncarbondenominators_finalfit.C
	./run_analysis_noncarbondenominators.sh

search.o: search.cpp
	g++ -O2 -c search.cpp -Wall -Wextra -lm

b12search.o: b12search.cpp search.h
	g++ -ffast-math -O3 -c b12search.cpp `root-config --cflags` -Wall -Wextra

stoptoli9.o: stoptoli9.cpp search.h
	g++ -O2 -c stoptoli9.cpp `root-config --cflags` -Wall -Wextra

b12search: b12search.o search.o
	g++ -o b12search b12search.o search.o `root-config --libs`

stoptoli9: stoptoli9.o search.o
	g++ -o stoptoli9 stoptoli9.o search.o `root-config --libs`

clean:
	rm -f stoptoli9 b12search stoptoli9.o b12search.o search.o \
      *_C.d *_C.so  AutoDict*cxx*

analysisclean:
	rm -f *_C.d *_C.so AutoDict*cxx* *.out *.out.h *.technote *.out.fail
