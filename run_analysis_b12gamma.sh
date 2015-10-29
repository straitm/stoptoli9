#!/bin/bash

out=b12gamma_finalfit.out 
set -o pipefail

root -b -q b12gammafinalfit.C+ | tee /tmp/$$.$out && 
mv -f /tmp/$$.$out $out
grep TECHNOTE $out > b12gamma_finalfit_out.technote
