#!/bin/bash

out=n16_finalfit.out 
set -o pipefail

root -b -q n16finalfit.C+O | tee /tmp/$$.$out &&
mv -f /tmp/$$.$out $out &&
grep TECHNOTE $out > n16_finalfit_out.technote
