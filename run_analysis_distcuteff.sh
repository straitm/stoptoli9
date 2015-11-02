#!/bin/bash

name=distcuteff
. analysis_function.sh

if ! root -l -n -b -q ${macro}+O &> $tmp; then
  fail $name
else
  finish $name
fi
