#!/bin/bash

(root -b -q c9finalfit.C"('o')";
root -b -q c9finalfit.C"('n')") | tee c9_finalfit.out

grep TECHNOTE c9_finalfit.out > c9_finalfit_out.technote
