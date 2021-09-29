#!/bin/bash
#SBATCH -N 1
#SBATCH -C haswell
#SBATCH -q regular
#SBATCH -J 273_519_0990_half
#SBATCH -A m3767
#SBATCH --mail-user=ahong@lbl.gov
#SBATCH --mail-type=ALL
#SBATCH -t 48:00:00

#OpenMP settings:
export OMP_NUM_THREADS=1
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

#run the application
srun -n 8 -c 8 --cpu_bind=cores ./main3d.gnu.haswell.MPI.ex input
