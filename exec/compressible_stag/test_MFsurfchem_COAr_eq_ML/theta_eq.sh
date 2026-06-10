#!/bin/bash

exec=../../../../FBoxLib/Tools/Postprocessing/F_Src/faverage.Linux.gfortran.exe
outfile=res.theta_eq
tmpfile=tmp.theta_eq
pyscr=theta_eq.py

if [ -f $outfile ]
then
    rm $outfile
fi

$exec -p plt000000000 -o $tmpfile -v 1 surfcov_0
python $pyscr $tmpfile >> $outfile

rm $tmpfile

cat $outfile
