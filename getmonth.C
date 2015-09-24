#include <fstream>

int getmonth(int run)
{
  ifstream infile("b12months");
  int frun, month;
  while(infile >> frun >> month)
    if(frun == run){
      fprintf(stderr, "%d\n", month);
      return month;
    }
}
