#!/bin/bash

out=be12_finalfit.out 
set -o pipefail

if ! root -b -q be12finalfit.C | tee /tmp/$$.$out; then
  mv -f /tmp/$$.$out $out.fail
  exit 1
fi

mv -f /tmp/$$.$out $out
grep TECHNOTE $out > be12_finalfit_out.technote
