#!/bin/bash

n=$1

out=/cp/s4/strait/fullfido-300s-3-25MeV-20151117.d/$n

mkdir -p $(dirname $out)

if ! [ -e $out ] || file $out | grep -q empty; then
  ../b12search far 0 300000 3 25 $(cat b12args$n) > $out && root -b -q ../b12ntuple2root.C+'("'$out'")' && xz -1 $out
fi
