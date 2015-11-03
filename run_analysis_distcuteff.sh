#!/bin/bash

name=distcuteff_$1
. analysis_function.sh

if ! root -l -n -b -q ${macro}+O &> $tmp; then
  fail $name
else
  if grep -q "segmentation violation" $tmp; then
    fail $name
  else
    finish $name
  fi
fi
