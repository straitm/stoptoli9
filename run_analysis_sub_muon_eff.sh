#!/bin/bash

name=sub_muon_eff
tmp=/tmp/$$.$name.out
tmp1=/tmp/$$.$name.out1
tmp2=/tmp/$$.$name.out2

# This is at least mostly CPU-bound
for n in $(cat runlist); do
  echo ./getliveanddeadtime.sh $n
done | parallel -j 10 | tee $tmp1 | awk '{\
  tot+=$2;\
  dead05+=$3;\
  dead10+=$4\
}\
END{\
  printf("TECHNOTE 4.1.2: subsequent muon efficiency for 0.5ms veto %.2f%%\n", 100*(1-dead05/tot));\
  printf("TECHNOTE 4.1.2: subsequent muon efficiency for 1.0ms veto %.2f%%\n", 100*(1-dead10/tot));\
  printf("const double sub_muon_eff05 = %.9f;\n", 1-dead05/tot);\
  printf("const double sub_muon_eff10 = %.9f;\n", 1-dead10/tot);\
}' &> $tmp2

if [ $? -eq 0 ]; then
  grep '^const ' $tmp2 > $name.out.h
  grep '^TECHNOTE' $tmp2 > $name.out.technote
  cat $tmp1 $tmp2 > $tmp
  rm -f $tmp1 $tmp2
  mv $tmp $name.out
else
  cat $tmp1 $tmp2 > $tmp
  rm -f $tmp1 $tmp2
  mv $tmp $name.fail
fi
