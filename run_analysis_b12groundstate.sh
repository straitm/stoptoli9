#!/bin/bash

out=b12groundstate_finalfit.out 
set -o pipefail

if ! root -b -q b12groundstate_finalfit.C | tee /tmp/$$.$out; then
  mv -f /tmp/$$.$out $out.fail
  exit 1
fi

mv -f /tmp/$$.$out $out
grep TECHNOTE $out > b12groundstate_finalfit_out.technote
