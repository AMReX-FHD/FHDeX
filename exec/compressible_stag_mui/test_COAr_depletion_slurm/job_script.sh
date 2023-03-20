#! /bin/bash -l
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=56
#SBATCH --partition dept.appliedmath
#SBATCH --time=140:00:00

# COMMANDS HERE

srun -n 112 -l --multi-prog mpmd.conf
