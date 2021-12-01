#!/bin/bash

NRUN=100

for ((i=1;i<=$NRUN;i++))
do
  echo R$i
  log=R$i/log.spparks
  res=R$i/res.coverage

  grep -B1 "Loop" $log --no-group-separator | awk 'NR % 2 == 1{printf("%e\t%d\t%d\t%e\n",$1,$6,$7,$7/($6+$7))}' > $res

  python3 coverage.py $res > ${res}_stat
done
