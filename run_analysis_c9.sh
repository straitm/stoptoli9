#!/bin/bash

name=c9
. analysis_function.sh

if ! root -b -q ${macro}"('o')" | tee $tmp ||
   ! root -b -q ${macro}"('n')" | tee -a $tmp; then
  fail $name
else
  finish $name
fi
