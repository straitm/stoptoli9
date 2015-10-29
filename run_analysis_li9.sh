#!/bin/bash

name=li9
macroparameter=$1
. analysis_function.sh

if ! root -b -q ${macro}+O'('$macroparameter')' | tee $tmp; then
  fail $name
else
  finish $name
fi
