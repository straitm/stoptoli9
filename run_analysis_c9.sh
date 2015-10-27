#!/bin/bash

out=c9_finalfit.out 
set -o pipefail

root -b -q c9finalfit.C"('o')" | tee /tmp/$$.$out &&
root -b -q c9finalfit.C"('n')" | tee -a /tmp/$$.$out &&
mv -f /tmp/$$.$out $out &&
grep TECHNOTE $out > c9_finalfit_out.technote
