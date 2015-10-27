#!/bin/bash

out=li8_finalfit.out 

root -b -q li8finalfit.C+ | tee /tmp/$$.$out &&
mv -f /tmp/$$.$out $out
grep TECHNOTE $out > li8_finalfit_out.technote
