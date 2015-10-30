#!/bin/bash

name=b12groundstate
. analysis_function.sh

if ! root -b -q $macro &> $tmp; then
  fail $name
else
  finish $name
fi
