#!/bin/bash

name=c9
. analysis_function.sh

if ! root -l -n -b -q ${macro}"('o')" &> $tmp ||
   ! root -l -n -b -q ${macro}"('n')" &>> $tmp; then
  fail $name
else
  finish $name
fi
