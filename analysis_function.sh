#!/bin/bash

set -o pipefail

if [ $macroparameter ]; then
  extrafilenamepart=_$macroparameter
fi
tmp=/tmp/$$.${name}_finalfit$extrafilenamepart.out 
macro=${name}_finalfit.C

finish(){
  name=$1
  out=${name}_finalfit$extrafilenamepart.out 
  headerout=${name}_finalfit$extrafilenamepart.out.h
  technoteout=${name}_finalfit$extrafilenamepart.out.technote

  mv -f $tmp $out
  grep TECHNOTE $out > $technoteout

  if grep -q '^const ' $out; then
    printf '/* Automatically generated by %s analysis */\n' $name > $headerout
    grep '^const ' $out >> $headerout
  fi

  echo $name $macroparameter analysis finished sucessfully

  true # Always suceed
}

fail(){
  name=$1
  out=${name}_finalfit.out 
  mv -f $tmp $out.fail
  cat $out.fail
  exit 1
}
