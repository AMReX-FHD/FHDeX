#!/bin/bash

SPKSCR=in.kmc_nonmui

# check kmc executable
exec=../SPPARKS_MUI/spk_nonmui
if [ ! -f $exec1 ]
then
  echo "ERROR: kmc executable $exec not found"
  exit
fi

echo "mpirun -np 16 $exec -var SEED 100 < $SPKSCR"
mpirun -np 16 $exec -var SEED 100 -screen none < $SPKSCR &

echo "try \"tail -f log.spparks\""
echo "try \"./coverage.sh"
