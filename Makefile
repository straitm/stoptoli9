all: stoptoli9 b12search nsearch fissionsearch microdst

microdst: microdst.cpp
	g++ -O2 -o microdst microdst.cpp `root-config --libs --cflags` -Wall -Wextra

fissionsearch: fissionsearch.cpp
	g++ -O2 -o fissionsearch fissionsearch.cpp `root-config --libs --cflags` -Wall -Wextra

nsearch: nsearch.cpp
	g++ -O2 -o nsearch nsearch.cpp `root-config --libs --cflags` -Wall -Wextra

b12search: b12search.cpp search.h
	g++ -ffast-math -O3 -o b12search b12search.cpp `root-config --libs --cflags` -Wall -Wextra

stoptoli9: stoptoli9.cpp search.h
	g++ -O2 -o stoptoli9 stoptoli9.cpp `root-config --libs --cflags` -Wall -Wextra

clean: 
	rm stoptoli9 b12search nsearch fissionsearch microdst
