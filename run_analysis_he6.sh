#!/bin/bash

name=he6
macroparameter=$1
. analysis_function.sh

if ! root -l -n -b -q ${macro}+O'('$macroparameter')' &> $tmp; then
  fail $name
else
  finish $name
fi
