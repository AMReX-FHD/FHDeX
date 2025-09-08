#! /bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --partition pi.ckim103
#SBATCH --nodelist=gnode011
#SBATCH --time=0-03:00:00

# COMMANDS HERE

srun -n 32 ../main3d.gnu.MPI.ex inputs_NO2_dimerization_eq
