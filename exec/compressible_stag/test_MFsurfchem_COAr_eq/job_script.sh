#! /bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --partition pi.ckim103
#SBATCH --time=0-15:00:00

# COMMANDS HERE

srun -n 64 ../main3d.gnu.MPI.ex inputs_fhd_stag
