#!/bin/bash

out=n16_finalfit.out 
set -o pipefail

if ! root -b -q n16finalfit.C+O | tee /tmp/$$.$out; then
  mv -f /tmp/$$.$out $out.fail
  exit 1
fi

mv -f /tmp/$$.$out $out
grep TECHNOTE $out > n16_finalfit_out.technote
