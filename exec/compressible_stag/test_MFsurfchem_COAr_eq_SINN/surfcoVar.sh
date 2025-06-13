#!/bin/bash

exec=/home/chemtae/GIT/FBoxLib/Tools/Postprocessing/F_Src/faverage.Linux.gfortran.exe

outfile=res.surfcoVar
tmpfile=tmp.surfcoVar
plotfile=surfcoVar.png
pyscr=surfcoVar.py
plotscr=surfcoVar.plt

pltfiles=`ls -d plt*00000`

if [ -f $outfile ]
then
    rm $outfile
fi

for pltfile in $pltfiles
do
    echo $pltfile >> $outfile
    $exec -p $pltfile -o $tmpfile -v 6 rhoYk-T theta-rhoYk theta-vx theta-vy theta-vz theta-T
    python $pyscr $tmpfile >> $outfile
done

rm $tmpfile

grep "Var\[rhoYk-T\]" $outfile     > ${outfile}_rhoYk-T
grep "Var\[theta-rhoYk\]" $outfile > ${outfile}_theta-rhoYk
grep "Var\[theta-vx\]" $outfile    > ${outfile}_theta-vx
grep "Var\[theta-vy\]" $outfile    > ${outfile}_theta-vy
grep "Var\[theta-vz\]" $outfile    > ${outfile}_theta-vz
grep "Var\[theta-T\]" $outfile     > ${outfile}_theta-T

gnuplot $plotscr
eog $plotfile
