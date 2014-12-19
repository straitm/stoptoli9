all: stoptoli9 b12search nsearch fissionsearch microdst

search.o: search.cpp
	g++ -O2 -c search.cpp -Wall -Wextra -lm

microdst: microdst.cpp
	g++ -O2 -o microdst microdst.cpp `root-config --libs --cflags` -Wall -Wextra

fissionsearch: fissionsearch.cpp
	g++ -O2 -o fissionsearch fissionsearch.cpp `root-config --libs --cflags` -Wall -Wextra

b12search.o: b12search.cpp search.h
	g++ -ffast-math -O3 -c b12search.cpp `root-config --cflags` -Wall -Wextra

stoptoli9.o: stoptoli9.cpp search.h
	g++ -O2 -c stoptoli9.cpp `root-config --cflags` -Wall -Wextra

b12search: b12search.o search.o
	g++ -o b12search b12search.o search.o `root-config --libs`

stoptoli9: stoptoli9.o search.o
	g++ -o stoptoli9 stoptoli9.o search.o `root-config --libs`

clean: 
	rm -f stoptoli9 b12search nsearch fissionsearch microdst stoptoli9.o b12search.o search.o
