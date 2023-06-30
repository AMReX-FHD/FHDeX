#! /bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=56
#SBATCH --partition dept.appliedmath

# COMMANDS HERE

SECONDS=0
srun -n 56 -l --multi-prog amrex_mpmd.conf
echo $SECONDS >> cpu_time
