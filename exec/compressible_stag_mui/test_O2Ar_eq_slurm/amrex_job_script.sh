#! /bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=41
#SBATCH --partition pi.ckim103
#SBATCH --time=0-03:00:00

# COMMANDS HERE

srun -n 41 -l --multi-prog amrex_mpmd.conf
