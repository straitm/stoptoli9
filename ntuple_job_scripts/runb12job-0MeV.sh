#!/bin/bash

n=$1

out=/cp/s4/strait/fullfido-100s-0-25MeV-20151117.d/$n

mkdir -p $(dirname $out)

if ! [ -e $out ] || file $out | grep -q empty; then
  ../b12search far 400 100000 0 25 $(cat b12args$n) > $out && root -b -q ../b12ntuple2root.C+'("'$out'")' && gzip $out
fi
