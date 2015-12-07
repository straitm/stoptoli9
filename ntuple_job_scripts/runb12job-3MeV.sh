#!/bin/bash

n=$1
shift
. sr

out=/cp/s4/strait/fullfido-300s-3-25MeV-20151203.d/$n

mkdir -p $(dirname $out)

if [ $(hostname) == cps4 ]; then
  B12SEARCH=../b12search
  TOROOT=/home/cp/strait/stoptoli9/b12ntuple2root.C
else
  B12SEARCH=../b12search-SL5
  TOROOT=/home/cp/strait/stoptoli9/b12ntuple2root_SL5.C
fi

cd ~strait/stoptoli9/ntuple_job_scripts/

if ! [ -e $out ] || file $out | grep -q empty; then
  $B12SEARCH far 0 300000 3 25 $(cat b12args$n) > $out && root -b -q ${TOROOT}+'("'$out'")' && gzip $out
fi
