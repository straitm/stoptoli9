#include "TTree.h"

void b12ntuple2root(const char * fn, bool hasheader)
{
  TTree t("t", "t");
  t.ReadFile(fn, hasheader?"":"trig/I:dt/F:dist/F:e/F:dx/F:dy/F:dz/F:"
         "run/I:mutrig/I:ovcoin/I:mx/F:my/F:mz/F:"
         "chi2/F:ivdedx/F:ngdnear/I:ngd/I:nnear/I:n/I:latengdnear/I:"
         "latengd/I:latennear/I:laten/I:miche/F:micht/F:gclen/F:"
         "fex/F:fey/F:fez/F:deadt/F:"
         "deade/F:michd/F:fq/F:fqiv/F:cq/F:cqiv/F:timeleft/F:"
         "ttlastvalid/F:ttlastmuon/F");
  t.SaveAs(Form("%s.root", fn));
}
