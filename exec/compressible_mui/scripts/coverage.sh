#!/bin/bash

log=log.spparks
res=res.coverage

if [ ! -f $log ]
then
    echo "ERROR: $log does not exist"
    exit
fi

grep -B1 "Loop" $log --no-group-separator | awk 'NR % 2 == 1{printf("%e\t%d\t%d\t%d\t%d\n",$1,$7,$8,$9,$10)}' > $res
