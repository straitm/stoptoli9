#!/bin/bash

# Runs the B-12-like analysis on the high-purity sample, giving its
# results. Runs it for the anti-high-purity sample, and uses the ratio
# to give the number of captures in the loose sample. Sends output for
# tech note to loosecaptures_finalfit.out and to be used downstream to
# loosecaptures_finalfit_out.h

root -b -q loosecaptures_finalfit.C+ | \
  tee loosecaptures_finalfit.out | \
  tee /dev/stderr | \
  grep 'const double' > loosecaptures_finalfit_out.h
