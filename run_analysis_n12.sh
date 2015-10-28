#!/bin/bash

out=n12_finalfit.out 
set -o pipefail

root -b -q n12finalfit.C | tee /tmp/$$.$out &&
mv -f /tmp/$$.$out $out &&
grep TECHNOTE $out > n12_finalfit_out.technote
