#!/bin/bash

run=$1
listout=$2

if [ -e $listout ]; then
  echo Not clobbering existing $listout
  exit 1
fi

sli9dir=/home/cp/strait/stoptoli9/
datadir=/cp/s4/strait/

root -l -b -n -q \
$sli9dir/preselectforfido.C'("'$datadir'/jplighttree/SEQ13/data.'$run'.root")' \
| grep -vE '^$|Processing' > $listout
