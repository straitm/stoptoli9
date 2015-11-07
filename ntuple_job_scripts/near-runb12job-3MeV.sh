#!/bin/bash

run=$1

out=/cp/s4/strait/fullfido-near.d/$run

mkdir -p $(dirname $out)

if ( ! [ -e $out ] || file $out | grep -q empty) && ! [ -e $out.xz ]; then
  ../b12search near 800 10000 3 25 /cp/s4/strait/jplighttree/SEQ12/data.$run.root /cp/s4/strait/ndfido_seq2/$run/$run.seq4.fido.root  > $out && root -b -q ../b12ntuple2root.C+'("'$out'", true)' && xz $out
fi
