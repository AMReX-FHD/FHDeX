#!/bin/bash

log=log.spparks
res=res.coverage

if [ ! -f $log ]
then
    echo "ERROR: $log does not exist"
    exit
fi

grep -B1 "Loop" $log --no-group-separator | awk 'NR % 2 == 1{printf("%e\t%d\t%d\n",$1,$6,$7,$8)}' > $res

gnuplot -persist -e 'plot "res.coverage" u 1:($3/2400/2400) w l'
