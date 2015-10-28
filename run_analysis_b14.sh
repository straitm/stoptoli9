#!/bin/bash

out=b14_finalfit.out 
set -o pipefail

root -b -q b14finalfit.C | tee /tmp/$$.$out &&
mv -f /tmp/$$.$out $out &&
grep TECHNOTE $out > b14_finalfit_out.technote
