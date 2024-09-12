#!/bin/bash

exec=../../../../FBoxLib/Tools/Postprocessing/F_Src/faverage.Linux.gfortran.exe

if [ "$#" -eq 0 ]
then 
    pltfiles=(`ls -d plt0*`)
    pltfile=${pltfiles[-1]}
else
    pltfile=$1
fi

$exec -p $pltfile -o res.zavg3 -v 10 rhoMean rhoEMean rhoYkMean_0 rhoYkMean_1 tMean pMean YkMean_0 YkMean_1 XkMean_0 XkMean_1

python entropy_check.py
