#!/bin/bash

exec=../../../../FBoxLib/Tools/Postprocessing/F_Src/faverage.Linux.gfortran.exe

if [ "$#" -eq 0 ]
then
    pltfiles=(`ls -d plt0*`)
    pltfile=${pltfiles[-1]}
else
    pltfile=$1
fi

$exec -p $pltfile -o res.zavg -v 8 tMean rhoMean rhoYkMean_0 rhoYkMean_1 TVarDirect rhoVar rhoYkVar_0 rhoYkVar_1

gnuplot zavg2.plt

mv res.zavg res.zavg2_${pltfile}
mv zavg.png zavg2_${pltfile}.png
eog zavg2_${pltfile}.png
