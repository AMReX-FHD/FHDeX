#!/bin/bash
#SBATCH --partition test
#SBATCH --time=0-00:01:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --export=ALL

#  type 'man sbatch' for more information and options

srun -n 6 --multi-prog mpmd.conf
