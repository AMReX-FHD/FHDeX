#!/bin/bash

exec=../../../../FBoxLib/Tools/Postprocessing/F_Src/faverage.Linux.gfortran.exe
outfile=res.mass_cons
tmpfile=tmp.mass_cons
plotfile=mass_cons.png
pyscr=mass_cons.py

pltfiles=`ls -d plt0*`

if [ -f $outfile ]
then
    rm $outfile
fi

for pltfile in $pltfiles
do
    $exec -p $pltfile -o $tmpfile -v 2 surfcov_0 rhoYkInstant_0
    python $pyscr $tmpfile >> $outfile
done

rm $tmpfile

cat $outfile
