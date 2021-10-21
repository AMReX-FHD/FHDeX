#!/bin/bash

NRUN=100

rm -f res.frac

for ((i=1;i<=$NRUN;i++))
do
  grep -A1 frac R$i/log.spparks --no-group-separator | tail -1 >> res.frac
done

python analyze_frac.py | tee res.frac_stat
