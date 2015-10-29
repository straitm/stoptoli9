#!/bin/bash

name=he6
. analysis_function.sh

if ! root -b -q ${macro}+O'(0)' | tee $tmp ||
   ! root -b -q ${macro}+O'(1)' | tee -a $tmp ||
   ! root -b -q ${macro}+O'(2)' | tee -a $tmp ||
   ! root -b -q ${macro}+O'(3)' | tee -a $tmp; then
  fail $name
else
  finish $name
fi
