#!/bin/bash

exec=../../../../FBoxLib/Tools/Postprocessing/F_Src/faverage.Linux.gfortran.exe
outfile=res.ads
tmpfile=tmp.ads
pyscr=ads.py

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
