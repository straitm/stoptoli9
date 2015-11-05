#!/bin/bash

n=$1

out=/cp/s4/strait/fullfido-muinfo-20151105.d/$n

mkdir -p $(dirname $out)

if ! [ -e $out ] || file $out | grep -q empty; then
  ../b12search far 0 0 0 0 $(cat b12args$n) > $out && root -b -q ../b12ntuple2root.C+'("'$out'", true)' && xz $out
fi
