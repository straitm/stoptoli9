#!/bin/bash

name=fullb12
. analysis_function.sh

if ! root -b -q ${macro}+O | tee $tmp; then
  fail $name
else
  finish $name
fi
