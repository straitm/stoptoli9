#!/bin/bash

n=$1

if [ "$2" == post3rdpub ]; then
  base=shortb12args-post3rdpub
  suffix=post3rdpub
else
  base=shortb12args
  suffix=
fi

out=/cp/s4/strait/fullfido-neutron-20150727.d/$n.$suffix

mkdir -p $(dirname $out)

while ! [ -e $out ] || file $out | grep -q empty; do
  ../neutronsearch far 0 0 3 25 $(cat $base$n) > $out
done
