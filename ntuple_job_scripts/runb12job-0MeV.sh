#!/bin/bash

n=$1
shift
. sr

out=/cp/s4/strait/fullfido-100s-0-25MeV-20151117.d/$n

mkdir -p $(dirname $out)

if [ $(hostname) == cps4 ]; then
  B12SEARCH=../b12search
else
  B12SEARCH=../b12search-SL5
fi

cd ~strait/stoptoli9/ntuple_job_scripts/

if ! [ -e $out ] || file $out | grep -q empty; then
  $B12SEARCH far 400 100000 0 25 $(cat b12args$n) > $out && root -b -q ../b12ntuple2root.C+'("'$out'")' && gzip $out
fi
