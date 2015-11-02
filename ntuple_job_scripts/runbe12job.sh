#!/bin/bash

n=$1

out=/cp/s4/strait/be12-20151102.d/$n

mkdir -p $(dirname $out)

if ! [ -e $out ] || file $out | grep -q empty; then
  ../be12search far 400 3000 3 25 $(cat b12args$n) > $out && root -b -q ../b12ntuple2root.C+'("'$out'", true)' && xz $out
fi
