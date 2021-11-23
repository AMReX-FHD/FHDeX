#!/bin/bash

log=log.spparks
res=res.coverage
scr=coverage_stat.py

if [ ! -f $log ]
then
    echo "ERROR: $log does not exist"
    exit
fi

grep -B1 "Loop" $log --no-group-separator | awk 'NR % 2 == 1{printf("%e\t%d\t%d\t%e\n",$1,$6,$7,$7/($6+$7))}' > $res

gnuplot -persist -e 'plot "res.coverage" u 1:4 w l'

python $scr $res
