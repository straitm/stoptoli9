#!/bin/bash

out=be12_finalfit.out 
set -o pipefail

root -b -q be12finalfit.C | tee /tmp/$$.$out &&
mv -f /tmp/$$.$out $out &&
grep TECHNOTE $out > be12_finalfit_out.technote
