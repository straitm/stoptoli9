#!/bin/bash

root -b -q b12like_finalfit.C++O'(bad)' &> /dev/null

cat run2month | awk '{print $2}' | sort | uniq | \
while read month; do
  printf "echo month: $month; root -b -q b12like_finalfit.C+'(\"run >= %s && run <= %s && e > 4 && e < 14 && dt < %%f\", %f, %f)'\n" \
   $(grep $month run2month | head -n 1 | awk '{print $1}') \
   $(grep $month run2month | tail -n 1 | awk '{print $1}') \
   $(grep $month month_livetime_efficiency | awk '{print $3, $2}')
done | parallel -j 1 2> /dev/null | grep -E "month:|per live" -A 2
