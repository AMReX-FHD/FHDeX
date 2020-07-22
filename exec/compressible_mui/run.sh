#!/bin/bash

# kmc executable
exec1=SPPARKS_MUI/spk_mui
if [ ! -f $exec1 ]
then
  echo "ERROR: spparks executable $exec1 not found"
  exit
fi

# fhd executable 
exec2=./main3d.gnu.DEBUG.MPI.ex
if [ ! -f $exec2 ]
then
  echo "ERROR: fhd executable $exec2 not found"
  exit
fi

# run the two executables simultaneously
spkscr=in.4spec
fhdscr=inputs_equil_3d
echo "** running kmc and fhd"
mpirun -np 1 $exec1 < $spkscr : -np 1 $exec2 $fhdscr
