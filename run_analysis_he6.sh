#!/bin/bash

out=he6_finalfit.out 
set -o pipefail

if ! root -b -q he6finalfit.C+O'(0)' | tee /tmp/$$.$out ||
   ! root -b -q he6finalfit.C+O'(1)' | tee -a /tmp/$$.$out ||
   ! root -b -q he6finalfit.C+O'(2)' | tee -a /tmp/$$.$out ||
   ! root -b -q he6finalfit.C+O'(3)' | tee -a /tmp/$$.$out; then
  mv -f /tmp/$$.$out $out.fail
  exit 1
fi

mv -f /tmp/$$.$out $out
grep TECHNOTE $out > he6_finalfit_out.technote
