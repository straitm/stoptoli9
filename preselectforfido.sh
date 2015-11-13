#!/bin/bash

run=$1
listout=$2
sli9dir=/home/cp/strait/stoptoli9/
datadir=/cp/s4/strait/

root -l -b -n -q \
$sli9dir/preselectforfido.C'("'$datadir'/jplighttree/SEQ13/data.'$run'.root")' \
| grep -vE '^$|Processing' > $listout
