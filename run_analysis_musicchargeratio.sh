#!/bin/bash

name=musicchargeratio
macroparameter=$1
. analysis_function.sh

#                                    true = far
if ! root -l -n -b -q ${macro}+O'("'$macroparameter'")' &> $tmp; then
  fail $name
else
  finish $name
fi
