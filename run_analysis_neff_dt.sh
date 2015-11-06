#!/bin/bash

name=neff_dt
. analysis_function.sh

if ! root -l -n -b -q ${macro}+O'(1)' &>  $tmp ||
   ! root -l -n -b -q ${macro}+O'(2)' &>> $tmp ||
   ! root -l -n -b -q ${macro}+O'(3)' &>> $tmp ||
   ! root -l -n -b -q ${macro}+O'(4)' &>> $tmp; then
  fail $name
else
  finish $name
fi
