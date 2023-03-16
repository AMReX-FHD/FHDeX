#!/bin/bash

exec=../../../../../FBoxLib/Tools/Postprocessing/F_Src/faverage.Linux.gfortran.exe
outfile=res.ads
tmpfile=tmp.ads
plotfile=ads.png
pyscr=ads.py

pltfiles=`ls -d plt0*`

if [ -f $outfile ]
then
    rm $outfile
fi

for pltfile in $pltfiles
do
    echo $pltfile >> $outfile
    $exec -p $pltfile -o $tmpfile -v 2 surfcovMean_0 surfcovVar_0
    python $pyscr $tmpfile >> $outfile
done

rm $tmpfile

grep "E\[theta\]" $outfile  > ${outfile}_mean
grep "Var\[theta\]" $outfile  > ${outfile}_var
