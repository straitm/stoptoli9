/*
 * Read in energy/time pairs and output columns as follows:
 *
 * Muon energy (muon = over 60MeV)
 * Time since muon up to 5.5mus
 * event energy
 *
 */

#include <stdio.h>
#include <iostream>
#include <vector>
using std::vector;
using std::cin;

struct in{
  float energy;
  double time;

  in(float e, double t){ energy = e; time = t; }
};

int main()
{
  /* quick and dirty solution: read in the whole damn thing so I can
  easily scan around. Of course, can't run this on huge files. */
  vector<in> indata;
  float ein = 0; double tin = 0;
  while(cin >> ein >> tin) indata.push_back(in(ein, tin));

  for(unsigned int i = 0; i < indata.size(); i++){
    if(indata[i].energy < 60) continue;
    for(unsigned int j = i+1; j < indata.size(); j++){
      if(indata[j].time > indata[i].time + 5500) break;
      if(indata[j].time - indata[i].time > 0) 
        printf("%f %f %f\n", indata[i].energy,
                             indata[j].time - indata[i].time,
                             indata[j].energy);
    }
  }
}
