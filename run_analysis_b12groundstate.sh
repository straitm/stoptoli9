#!/bin/bash

out=b12groundstate_finalfit.out 
set -o pipefail

root -b -q b12groundstate_finalfit.C | tee /tmp/$$.$out &&
mv -f /tmp/$$.$out $out &&
grep TECHNOTE $out > b12groundstate_finalfit_out.technote
