#!/bin/bash
#SBATCH --qos=debug
#SBATCH --time=30
#SBATCH --nodes=1
#SBATCH --tasks-per-node=32
#SBATCH --constraint=haswell

RUNDIR=RUN_SF

if [ -d $RUNDIR ]
then
  echo "ERROR: $RUNDIR already exists"
  exit
fi

mkdir $RUNDIR
cd $RUNDIR

srun -n5 -l --multi-prog ../mpmd.conf 
