#!/bin/bash

RUNDIR=RUN
SPKSCR=in.kmc
FHDSCR=inputs_fhd

# check kmc executable
exec1=SPPARKS_MUI/spk_mui
if [ ! -f $exec1 ]
then
  echo "ERROR: kmc executable $exec1 not found"
  exit
fi

# check fhd executable 
exec2=./main3d.gnu.DEBUG.MPI.ex
if [ ! -f $exec2 ]
then
  echo "ERROR: fhd executable $exec2 not found"
  exit
fi

# check RUNDIR and create RUNDIR
if [ -d $RUNDIR ]
then
  echo "ERROR: $RUNDIR already exists"
  exit
fi

# copy scripts
mkdir $RUNDIR
cp $SPKSCR $RUNDIR
cp $FHDSCR $RUNDIR
cd $RUNDIR

# run the two executables simultaneously
echo "** running kmc and fhd"
time mpirun -np 1 ../$exec1 -var SEED 100 -screen none < $SPKSCR : -np 1 ../$exec2 $FHDSCR | tee log.fhd

../scripts/coverage.sh
