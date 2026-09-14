#!/bin/bash

SPKSCR=in.kmc_nonmui
NRUN=100

# check kmc executable
exec=../SPPARKS_MUI/spk_nonmui
if [ ! -f $exec1 ]
then
  echo "ERROR: kmc executable $exec not found"
  exit
fi

for ((i=1;i<=$NRUN;i++))
do
    echo "***** R$i *****"
    mkdir R$i
    mpirun -np 4 $exec -var SEED $i -log R$i/log.spparks < $SPKSCR
done
