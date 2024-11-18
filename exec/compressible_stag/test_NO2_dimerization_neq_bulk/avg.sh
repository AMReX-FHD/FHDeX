#!/bin/bash

exec=../../../../FBoxLib/Tools/Postprocessing/F_Src/faverage.Linux.gfortran.exe

outfile=res.avg

if [ -f $outfile ]
then
    rm $outfile
fi

pltfiles=`ls -d plt[0-9]*`

for pltfile in $pltfiles
do
    resfile=${outfile}_${pltfile}
    $exec -p $pltfile -v 3 tInstant rhoYkInstant_0 rhoYkInstant_1 -o $resfile

    tail $resfile -n+3 | awk '{ sum1+=$2; sum1sq+=$2*$2; sum2+=$4; sum2sq+=$4*$4; sum3+=$6; sum3sq+=$6*$6 } END {print sum1/NR, sqrt((sum1sq-sum1*sum1/NR)/NR/(NR-1)), sum2/NR, sqrt((sum2sq-sum2*sum2/NR)/NR/(NR-1)), sum3/NR, sqrt((sum3sq-sum3*sum3/NR)/NR/(NR-1)) }' >> $outfile
done

gnuplot -persist <<-EOFMarker
    plot "res.avg" u 0:1 w lp
EOFMarker
