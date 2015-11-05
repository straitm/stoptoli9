#!/bin/bash

name=b8cutefficiency
. analysis_function.sh

if ! root -l -n -b -q ${macro}+O'(4)' &> $tmp; then
  fail $name
else
  finish $name
fi
