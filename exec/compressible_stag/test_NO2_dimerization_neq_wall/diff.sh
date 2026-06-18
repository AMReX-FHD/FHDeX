#!/bin/bash

DIR1=./
DIR2=../test_NO2_dimerization_neq_wall
FILES="submit_job.sh job_script.sh inputs_NO2_dimerization_neq_wall zavg.sh zavg2.sh zavg3.sh zavg.plt zavg2.plt entropy_check.py entropy_check_funcs.py clean.sh diff.sh README"

for file in $FILES
do
  echo "** $file"
  diff $DIR1/$file $DIR2/$file
done
