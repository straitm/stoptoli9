#!/bin/bash

out=b8_finalfit.out 
set -o pipefail

if ! root -b -q b8finalfit.C | tee /tmp/$$.$out; then
  mv -f /tmp/$$.$out $out.fail
  exit 1
fi

mv -f /tmp/$$.$out $out
grep TECHNOTE $out > b8_finalfit_out.technote
