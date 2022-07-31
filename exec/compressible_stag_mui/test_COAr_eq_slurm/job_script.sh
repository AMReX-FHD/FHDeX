#! /bin/bash -l
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=55
#SBATCH --partition dept.appliedmath
#SBATCH --time=90:00:00

# COMMANDS HERE

srun -n 110 -l --multi-prog mpmd.conf
