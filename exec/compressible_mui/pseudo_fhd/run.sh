#!/bin/bash

# kmc executable
exec1=../SPPARKS_MUI/spk_mui
if [ ! -f $exec1 ]
then
  echo "ERROR: kmc executable $exec1 not found"
  exit
fi

# pseudo-fhd executable
exec2=./pseudo_fhd.bin
if [ ! -f $exec2 ]
then
  echo "ERROR: pseudo-fhd executable $exec2 not found"
  echo "ERROR: run ./compile.sh"
  exit
fi

# run the two executables simultaneously
spkscr=in.4spec
echo "** running kmc and pseudo-fhd"
time mpirun -np 1 $exec1 -var SEED 100 < $spkscr : -np 1 $exec2

../coverage.sh
