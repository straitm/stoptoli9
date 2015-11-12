#include "TTree.h"

int b12ntuple2root(const char * fn)
{
  TTree t("t", "t");
  if(!t.ReadFile(fn)) return 1;
  t.SaveAs(Form("%s.root", fn));
  return 0;
}
