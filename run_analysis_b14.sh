#!/bin/bash

out=b14_finalfit.out 
set -o pipefail

if ! root -b -q b14finalfit.C | tee /tmp/$$.$out; then
  mv -f /tmp/$$.$out $out.fail
  exit 1
fi

mv -f /tmp/$$.$out $out
grep TECHNOTE $out > b14_finalfit_out.technote
