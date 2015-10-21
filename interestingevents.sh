#!/bin/bash

(root -b /home/cp/strait/stoptoli9/seq1.ntuple.root << EOF

t->Scan("mutrig", "dt > 0 && dt < 30 && e > 4 && e < 16 && fqiv < 200e3 && miche < 12 && abs(dz) > 10 && abs(dz) < 1700 && dx**2+dy**2 < 1700**2 && run == $1");

EOF
) | grep '^\* ' | awk '{print $4}' | grep '[1-9]'
