#! /bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --partition test
#SBATCH --time=0-00:30:00

# COMMANDS HERE

srun -n 16 ../main2d.gnu.MPI.ex inputs_Schlogl_2d
