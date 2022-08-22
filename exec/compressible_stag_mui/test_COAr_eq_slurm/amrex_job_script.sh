#! /bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=56
#SBATCH --partition test
#SBATCH --time=0-00:30:00

# COMMANDS HERE

srun -n 56 -l --multi-prog amrex_mpmd.conf
