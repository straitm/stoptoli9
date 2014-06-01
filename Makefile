all: stoptoli9

stoptoli9: stoptoli9.cpp
	g++ -O2 -o stoptoli9 stoptoli9.cpp `root-config --libs --cflags` -Wall -Wextra
