#!/bin/bash

name=c9
. analysis_function.sh

if ! root -b -q ${macro}"('o')" &> $tmp ||
   ! root -b -q ${macro}"('n')" &>> $tmp; then
  fail $name
else
  finish $name
fi
