#!/bin/bash

name=distcuteff
. analysis_function.sh

if ! root -l -n -b -q ${macro}+ &> $tmp; then
  fail $name
else
  finish $name
fi
