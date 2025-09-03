#!/bin/bash

SPKSCR=in.kmc_amrex
FHDSCR=inputs_fhd_stag_amrex
SPKSEED=$RANDOM

# check kmc executable
exec1=../spk3d.gnu.MPI.ex
if [ ! -f $exec1 ]
then
  echo "ERROR: kmc executable $exec1 not found"
  exit
fi

# check fhd executable
exec2=../fhd3d.gnu.MPI.ex
if [ ! -f $exec2 ]
then
  echo "ERROR: fhd executable $exec2 not found"
  exit
fi

# check number of steps
N1=`grep nstep $SPKSCR | head -1 | awk '{print $4}'`
N2=`grep max_step $FHDSCR | awk '{print $3}'`
if [ "$N1" != "$N2" ]
then
  echo "ERROR: nstep = $N1 (kmc) and max_step = $N2 (fhd) do not match"
  exit
fi

# check timestep size
N1=`grep Trun $SPKSCR | head -1 | awk '{print $4}'`
N2=`grep fixed_dt $FHDSCR | awk '{print $3}'`
if [ "$N1" != "$N2" ]
then
  echo "ERROR: Trun = $N1 (kmc) and fixed_dt = $N2 (fhd) do not match"
  exit
fi

# run the two executables simultaneously
echo "mpirun -np 2 $exec1 -var SEED $SPKSEED -screen none < $SPKSCR : -np 2 $exec2 $FHDSCR > log.fhd &"
mpirun -np 2 $exec1 -var SEED $SPKSEED -screen none < $SPKSCR : -np 2 $exec2 $FHDSCR > log.fhd &

echo "try \"tail -f log.fhd\""
