all: stoptoli9 b12search

analysis: fullb12_finalfit.out c9_finalfit.out li8_finalfit.out

c9_finalfit.out: c9finalfit.C loosecaptures_finalfit_out.h dcfluids_finalfit_out.h
	./run_analysis_c9.sh

li8_finalfit.out: li8finalfit.C loosecaptures_finalfit_out.h
	./run_analysis_li8.sh

fullb12_finalfit.out: fullb12finalfit.C loosecaptures_finalfit_out.h
	./run_analysis_fullb12.sh

loosecaptures_finalfit_out.h: consts.h mucountfinalfit.C b12like_finalfit.C
	./run_analysis_b12like.sh

dcfluids_finalfit_out.h: dcfluids_finalfit.C consts.h loosecaptures_finalfit_out.h
	./run_dcfluids.sh

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
	rm *.out *_out.h *.technote
