#!/bin/bash

exec=../../../../FBoxLib/Tools/Postprocessing/F_Src/faverage.Linux.gfortran.exe

pltfiles=(`ls -d plt0*`)
pltfile=${pltfiles[-1]}

outfile=res.zavg_rhoY1mean

$exec -p $pltfile -o $outfile -v 1 rhoYkMean_0

gnuplot zavg_rhoY1mean.plt
eog zavg_rhoY1mean.png
