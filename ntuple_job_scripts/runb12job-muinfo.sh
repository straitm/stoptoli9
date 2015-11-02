#!/bin/bash

n=$1

out=/cp/s4/strait/fullfido-511keV-20150911.d/$n

mkdir -p $(dirname $out)

while ! [ -e $out ] || file $out | grep -q empty; do
  ../buffersearch far 0 0 3 25 $(cat shortb12args$n) > $out
done

