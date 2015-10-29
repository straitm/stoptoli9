#!/bin/bash

name=be12
. analysis_function.sh

if ! root -b -q $macro | tee $tmp; then
  fail $name
else
  finish $name
fi
