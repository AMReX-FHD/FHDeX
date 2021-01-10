#!/bin/bash

log=log.spparks
res=res.coverage

if [ ! -f $log ]
then
    echo "ERROR: $log does not exist"
    exit
fi

grep -B1 "Loop" $log --no-group-separator | awk 'NR % 2 == 1{printf("%e\t%d\t%d\n",$1,$6,$7,$8)}' > $res

gnuplot -persist -e 'plot "res.coverage" u 1:3 w l'
#gnuplot -persist -e 'plot [0:1e-7] "res.coverage" u 1:3 w l,1.519198e-01*(1-exp(-3.451071e+07*x))*1200*1200'
