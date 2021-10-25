#!/bin/bash

exec=../../../../FBoxLib/Tools/Postprocessing/F_Src/faverage.Linux.gfortran.exe

pltfiles=(`ls -d plt0*`)
pltfile=${pltfiles[-1]}

outfile=res.zavg

$exec -p $pltfile -o $outfile -v 6 rhoVar rhoEVar tVar jxVar jyVar jzVar
