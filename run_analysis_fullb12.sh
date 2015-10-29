#!/bin/bash

out=fullb12_finalfit.out 
set -o pipefail

if ! root -b -q fullb12finalfit.C+O | tee /tmp/$$.$out; then
  mv -f /tmp/$$.$out $out.fail
  exit 1
fi

mv -f /tmp/$$.$out $out

printf '/* Automatically generated by fullb12_finalfit.C */\n' \
  > fullb12_finalfit_out.h
grep '^const ' $out >> fullb12_finalfit_out.h

grep TECHNOTE $out > fullb12_finalfit_out.technote
