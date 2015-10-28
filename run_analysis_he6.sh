#!/bin/bash

out=he6_finalfit.out 
set -o pipefail

root -b -q he6finalfit.C+'(0)' | tee /tmp/$$.$out &&
root -b -q he6finalfit.C+'(1)' | tee -a /tmp/$$.$out &&
root -b -q he6finalfit.C+'(2)' | tee -a /tmp/$$.$out &&
root -b -q he6finalfit.C+'(3)' | tee -a /tmp/$$.$out &&
mv -f /tmp/$$.$out $out &&
grep TECHNOTE $out > he6_finalfit_out.technote
