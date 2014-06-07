all: stoptoli9 b12search

b12search: b12search.cpp
	g++ -O2 -o b12search b12search.cpp `root-config --libs --cflags` -Wall -Wextra

stoptoli9: stoptoli9.cpp
	g++ -O2 -o stoptoli9 stoptoli9.cpp `root-config --libs --cflags` -Wall -Wextra
