#!/bin/bash

run=$1

live=$(root -b -q getrunlivetime.C'('$run')' 2> /dev/null | grep -vE '[*A-Z]|^$')
dead=$(./getmuondeadtime.sh $run)

echo $run $live $dead $(dc <<< "9k$live $dead -p")
