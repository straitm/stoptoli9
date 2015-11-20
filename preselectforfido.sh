#!/bin/bash

run=$1
listout=$2

if [ -e $listout ]; then
  echo Not clobbering existing $listout
  exit 1
fi

sli9dir=/home/cp/strait/stoptoli9/
datadir=/cp/s4/strait/

JP=$datadir/jplighttree/SEQ13/data.$run.root
if ! [ -e $JP ]; then
  JP=$datadir/jplighttree/SEQ12/data.$run.root
fi
if ! [ -e $JP ]; then
  echo Cannot find a JP file for this run, exiting
  exit 1
fi


root -l -b -n -q $sli9dir/preselectforfido.C'("'$JP'")' \
| grep -vE '^$|Processing' > $listout
