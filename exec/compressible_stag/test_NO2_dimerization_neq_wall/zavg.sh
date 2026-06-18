#!/bin/bash

exec=../../../../FBoxLib/Tools/Postprocessing/F_Src/faverage.Linux.gfortran.exe

if [ "$#" -eq 0 ]
then
    pltfiles=(`ls -d plt0*`)
    pltfile=${pltfiles[-1]}
else
    pltfile=$1
fi

$exec -p $pltfile -o res.zavg -v 4 tInstant rhoInstant rhoYkInstant_0 rhoYkInstant_1

gnuplot zavg.plt

mv res.zavg res.zavg_${pltfile}
mv zavg.png zavg_${pltfile}.png
eog zavg_${pltfile}.png
