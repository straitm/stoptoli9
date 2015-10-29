#!/bin/bash
# Runs the B-12-like analysis on the high-purity sample, giving its
# results. Runs it for the anti-high-purity sample, and uses the ratio
# to give the number of captures in the loose sample.

name=carbondenominators
. analysis_function.sh

if ! root -b -q ${macro}+O | tee $tmp; then
  fail $name
else
  finish $name
fi
