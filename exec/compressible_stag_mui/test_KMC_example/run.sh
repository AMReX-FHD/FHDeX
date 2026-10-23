#!/bin/bash
#SBATCH -J test_O2ads
#SBATCH -N 1
#SBATCH -n 32
#SBATCH -p pi.ckim103
#SBATCH -t 0-03:00:00

#run the application: 
srun ../SPPARKS_MUI/spk_nonmui -in in.kmc_O2 -log log.spparks
