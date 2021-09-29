#!/bin/bash
#SBATCH -N 1
#SBATCH -C haswell
#SBATCH -q regular
#SBATCH -J 273_396_0999
#SBATCH -A m3767
#SBATCH --mail-user=ahong@lbl.gov
#SBATCH --mail-type=ALL
#SBATCH -t 48:00:00

#OpenMP settings:
export OMP_NUM_THREADS=1
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

#run the application
srun -n 16 -c 4 --cpu_bind=cores ./main3d.gnu.haswell.MPI.ex input
