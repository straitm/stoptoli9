#!/bin/bash

name=b12gamma
macroparameter=$1
. analysis_function.sh

if ! root -b -q ${macro}+O'('$macroparameter')' &> $tmp; then
  fail $name
else
  finish $name
fi