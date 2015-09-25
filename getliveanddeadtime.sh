#!/bin/bash

run=$1

live=$(root -b -q getrunlivetime.C'('$run')' 2> /dev/null | grep -vE '[*A-Z]|^$')
dead=$(./getmuondeadtime.sh $run)
halfdead=$(echo $dead | awk '{print $1}')
onedead=$(echo $dead | awk '{print $2}')

echo $run $live $halfdead $onedead
