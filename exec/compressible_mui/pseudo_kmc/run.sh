#!/bin/bash

# pseudo-kmc executable
exec1=./pseudo_kmc.bin
if [ ! -f $exec1 ]
then
  echo "ERROR: pseudo-kmc executable $exec1 not found"
  echo "ERROR: run ./compile.sh"
  exit
fi

# fhd executable
exec2=../main3d.gnu.DEBUG.MPI.ex
if [ ! -f $exec2 ]
then
  echo "ERROR: fhd executable $exec2 not found"
  exit
fi

# run the two executables simultaneously
fhdscr=inputs_equil_3d
echo "** running pseudo-kmc and fhd"
mpirun -np 1 $exec1 : -np 1 $exec2 $fhdscr
