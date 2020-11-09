#!/bin/bash
#SBATCH --qos=regular
#SBATCH --time=180
#SBATCH --nodes=3
#SBATCH --tasks-per-node=32
#SBATCH --constraint=haswell

RUNDIR=RUN
KMCSCR=in.kmc
FHDSCR=inputs_fhd

if [ -d $RUNDIR ]
then
  echo "ERROR: $RUNDIR already exists"
  exit
fi

mkdir $RUNDIR
cp $KMCSCR $RUNDIR
cp $FHDSCR $RUNDIR
cp $0 $RUNDIR
cp mpmd.conf $RUNDIR
cd $RUNDIR

echo "*** START: `date`"

srun -n 80 -l --multi-prog ../mpmd.conf

echo "*** FINISH: `date`"
