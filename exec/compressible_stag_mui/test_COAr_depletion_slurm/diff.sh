#!/bin/bash

DIR1=./
DIR2=../test_COAr_depletion_slurm
FILES="submit_job.sh job_script.sh mpmd.conf inputs_fhd_stag in.kmc zavg.sh ads.sh mass_cons.py mass_cons.sh init_rho.py clean.sh diff.sh"

for file in $FILES
do
  echo "** $file"
  diff $DIR1/$file $DIR2/$file
done
