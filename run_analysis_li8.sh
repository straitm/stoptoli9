#!/bin/bash

name=li8
. analysis_function.sh

if ! root -b -q ${macro}+O | tee $tmp; then
  fail $name
else
  finish $name
fi
