#!/bin/bash

exec=../../../../FBoxLib/Tools/Postprocessing/F_Src/faverage.Linux.gfortran.exe

pltfiles=(`ls -d plt0*`)
pltfile=${pltfiles[-1]}

outfile=res.zavg

$exec -p $pltfile -o $outfile -v 9 rhoVar rhoEVar TVar jxVarFACE jyVarFACE jzVarFACE jxVarCC jyVarCC jzVarCC

gnuplot zavg.plt
eog zavg.png
