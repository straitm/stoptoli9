#!/bin/bash

out=b8_finalfit.out 
set -o pipefail

root -b -q b8finalfit.C | tee /tmp/$$.$out &&
mv -f /tmp/$$.$out $out &&
grep TECHNOTE $out > b8_finalfit_out.technote
