#!/bin/bash

out=fullb12_finalfit.out 
set -o pipefail

root -b -q fullb12finalfit.C+ | tee /tmp/$$.$out && 
mv -f /tmp/$$.$out $out
grep TECHNOTE $out > fullb12_finalfit_out.technote
