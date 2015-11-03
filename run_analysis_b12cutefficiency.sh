#!/bin/bash

name=b12cutefficiency
. analysis_function.sh

if ! root -l -n -b -q ${macro}+O'(4)' &> $tmp ||
   ! root -l -n -b -q ${macro}+O'(3)' &>> $tmp; then
  fail $name
else
  finish $name
fi
