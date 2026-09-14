#!/bin/bash

RUNDIR=RUN
SPKSCR=in.kmc
FHDSCR=inputs_fhd_coverage

# check fhd executable
exec1=../main3d.gnu.MPI.ex
if [ ! -f $exec1 ]
then
  echo "ERROR: fhd executable $exec1 not found"
  exit
fi

# check kmc executable
exec2=../SPPARKS_MUI/spk_mui
if [ ! -f $exec2 ]
then
  echo "ERROR: kmc executable $exec2 not found"
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
cp $SPKSCR $RUNDIR
cp $FHDSCR $RUNDIR
cp ../main_driver.cpp $RUNDIR
cd $RUNDIR

# check number of steps
N1=`grep nstep $SPKSCR | head -1 | awk '{print $4}'`
N2=`grep max_step $FHDSCR | awk '{print $3}'`
if [ "$N1" != "$N2" ]
then
  echo "ERROR: nstep = $N1 (kmc) and max_step = $N2 (fhd) do not match"
  exit
fi

# run the two executables simultaneously
echo "** running kmc and fhd"
time mpirun -np 4 ../$exec1 $FHDSCR : -np 1 ../$exec2 -var SEED 100 -screen none -in $SPKSCR | tee log.fhd

../../scripts/coverage.sh
