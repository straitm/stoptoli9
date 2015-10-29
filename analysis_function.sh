#!/bin/bash

set -o pipefail
tmp=/tmp/$$.${name}_finalfit.out 
macro=${name}_finalfit.C

finish(){
  name=$1
  out=${name}_finalfit.out 
  headerout=${name}_finalfit_out.h

  mv -f /tmp/$$.$out $out
  grep TECHNOTE $out > $technoteout

  if grep -q '^const ' $out; then
    printf '/* Automatically generated by %s analysis */\n' $name > $headerout
    grep '^const ' $out >> $headerout
  fi

  true # Always suceed
}

fail(){
  name=$1
  out=${name}_finalfit.out 
  mv -f /tmp/$$.$out $out.fail
  exit 1
}
