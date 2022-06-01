#!/bin/bash

exec=../../../../FBoxLib/Tools/Postprocessing/F_Src/faverage.Linux.gfortran.exe
outfile=res.ads
tmpfile=tmp.ads
plotfile=ads.png
pyscr=ads.py
plotscr=ads.plt

pltfiles=`ls -d plt0*`

if [ -f $outfile ]
then
    rm $outfile
fi

for pltfile in $pltfiles
do
    echo $pltfile >> $outfile
    $exec -p $pltfile -o $tmpfile -v 4 surfcovMean_0 surfcovMean_1 surfcovVar_0 surfcovVar_1
    python $pyscr $tmpfile >> $outfile
done

rm $tmpfile

grep "E\[theta1\]" $outfile  > ${outfile}_mean1
grep "Var\[theta1\]" $outfile  > ${outfile}_var1
grep "E\[theta2\]" $outfile  > ${outfile}_mean2
grep "Var\[theta2\]" $outfile  > ${outfile}_var2

gnuplot $plotscr
eog $plotfile
