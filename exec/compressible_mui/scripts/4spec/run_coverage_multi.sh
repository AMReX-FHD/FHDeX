#!/bin/bash

RUNDIR=RUN
SPKSCR=in.kmc
FHDSCR=inputs_fhd_coverage
NSAMPLE=64

# check kmc executable
exec1=SPPARKS_MUI/spk_mui
if [ ! -f $exec1 ]
then
  echo "ERROR: kmc executable $exec1 not found"
  exit
fi

# check fhd executable
exec2=./main3d.gnu.MPI.ex
if [ ! -f $exec2 ]
then
  echo "ERROR: fhd executable $exec2 not found"
  exit
fi

for ((i=1;i<=$NSAMPLE;i++))
do
  RUNDIRi=$RUNDIR$i
  SEED=$(( $i * 100 ))

  # check RUNDIRi and create RUNDIRi
  if [ -d $RUNDIRi ]
  then
    echo "ERROR: $RUNDIRi already exists"
  exit
  fi

  # copy scripts
  mkdir $RUNDIRi
  cp $0 $RUNDIRi
  cp $SPKSCR $RUNDIRi
  cp $FHDSCR $RUNDIRi
  cp main_driver.cpp $RUNDIRi
  cd $RUNDIRi

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
  time mpirun -np 1 ../$exec1 -var SEED $SEED -screen none < $SPKSCR : -np 1 ../$exec2 $FHDSCR | tee log.fhd

  ../scripts/coverage.sh

  cd ..
done
