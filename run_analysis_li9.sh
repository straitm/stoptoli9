#!/bin/bash

name=li9
. analysis_function.sh

if ! root -b -q ${macro}+O'(-1)' | tee $tmp ||
   ! root -b -q ${macro}+O'(1)'  | tee -a $tmp ||
   ! root -b -q ${macro}+O'(0)'  | tee -a $tmp; then
  fail $name
else
  finish $name
fi
