#!/bin/bash

RUNDIR=RUN_SF
FHDSCR=inputs_fhd_SF

# check fhd executable
exec2=./main3d.gnu.MPI.ex
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
cp $0 $RUNDIR
cp $FHDSCR $RUNDIR
cp main_driver.cpp $RUNDIR
cd $RUNDIR

# run
echo "** running fhd"
time mpirun -np 4 ../$exec2 $FHDSCR | tee log.fhd
