#!/bin/bash

run=$1

out=/cp/s4/strait/fullfido-near-neutron.d/neutron-$run

mkdir -p $(dirname $out)

while ! [ -e $out ] || file $out | grep -q empty; do
  ../neutronsearch near 0 0 3 25 /cp/s4/strait/jplighttree/SEQ13/data.$run.root /cp/s4/strait/ndfido/$run/$run.seq5.fido.root  > $out && root -b -q ../b12ntuple2root.C+'("'$out'")' && xz $out
done
