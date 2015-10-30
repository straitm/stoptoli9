#!/bin/bash

name=b14
. analysis_function.sh

if ! root -b -q $macro &> $tmp; then
  fail $name
else
  finish $name
fi
