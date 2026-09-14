#! /bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --partition test
#SBATCH --time=0-00:30:00

# COMMANDS HERE

srun -n 32 ../main3d.gnu.MPI.ex inputs_NO2_dimerization_neq_bulk0

srun echo "*** equilibration is done!"
srun echo "*** waiting 30 sec..."
srun sleep 30

srun -n 32 ../main3d.gnu.MPI.ex inputs_NO2_dimerization_neq_bulk1
