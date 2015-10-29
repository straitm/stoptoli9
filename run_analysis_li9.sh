#!/bin/bash

out=li9_finalfit.out 
set -o pipefail

root -b -q li9finalfit.C+O'(-1)' | tee /tmp/$$.$out &&
root -b -q li9finalfit.C+O'(1)'  | tee -a /tmp/$$.$out &&
root -b -q li9finalfit.C+O'(0)'  | tee -a /tmp/$$.$out &&
mv -f /tmp/$$.$out $out &&
grep TECHNOTE $out > li9_finalfit_out.technote
