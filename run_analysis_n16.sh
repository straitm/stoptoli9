#!/bin/bash

name=n16
. analysis_function.sh

if ! root -b -q ${macro}+O &> $tmp; then
  fail $name
else
  finish $name
fi
