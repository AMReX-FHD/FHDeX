#!/bin/bash

DIR1=./
DIR2=../test_O2Ar_eq_slurm
FILES="submit_job.sh amrex_job_script.sh inputs_fhd_stag params_O2Ar.py zavg.sh zavg.plt ads.sh ads.py ads.plt clean.sh diff.sh README"

for file in $FILES
do
  echo "** $file"
  diff $DIR1/$file $DIR2/$file
done
