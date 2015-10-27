#!/bin/bash

root -b -q fullb12finalfit.C+

grep TECHNOTE fullb12_finalfit.out > fullb12_finalfit_out.technote
