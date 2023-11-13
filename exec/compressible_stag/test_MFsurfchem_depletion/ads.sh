#!/bin/bash

exec=../../../../FBoxLib/Tools/Postprocessing/F_Src/faverage.Linux.gfortran.exe
outfile=res.ads
tmpfile=tmp.ads

pltfiles=`ls -d plt0*`

if [ -f $outfile ]
then
    rm $outfile
fi

for pltfile in $pltfiles
do
    $exec -p $pltfile -o $tmpfile -v 1 surfcov_0
    head -3 $tmpfile | tail -1 >> $outfile
done

rm $tmpfile

gnuplot -persist <<-EOFMarker
    plot [:][0:0.25] "$outfile" u 0:2 w l, 0.205
EOFMarker
