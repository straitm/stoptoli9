/*
  Takes in a b12search output file and makes all possible pairs of
  events following each muon in the style of be12search.
*/
#include <stdio.h>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
using namespace std;

const double maxtime = 10000; // in ms

void printandclearbuf(vector<string> & buf)
{
  for(unsigned int i = 0; i < buf.size(); i++)
    for(unsigned int j = i+1; j < buf.size(); j++)
      printf("%s %s\n", buf[i].c_str(), buf[j].c_str());
  buf.resize(0);
}

int main()
{
  string ins, oldrun, oldtrig;
  double dt;
  vector<string> buf;
  while(getline(cin, ins)){
    string murun, mutrig;
    stringstream ss(ins);

    // Throw out first field
    for(int i = 0; i < 1; i++) ss >> murun;
    ss >> dt;
    if(dt > maxtime) continue;
    // throw out the next seven fields
    for(int i = 0; i < 7; i++) ss >> murun;

    if(!(ss >> murun >> mutrig)){
      cerr << "Got line with fewer than eleven fields, ignoring\n";
      continue;
    }

    if(murun != oldrun || mutrig != oldtrig) printandclearbuf(buf);
    oldrun = murun;
    oldtrig = mutrig;
    buf.push_back(ins);
  }
  printandclearbuf(buf);
}
