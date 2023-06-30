#!/bin/bash

DIR1=./
DIR2=../test_NO2_dimerization_neq
FILES="submit_job.sh job_script.sh inputs_NO2_dimerization_neq params_NO2_dimerization_rate_const.py avg.sh clean.sh diff.sh README"

for file in $FILES
do
  echo "** $file"
  diff $DIR1/$file $DIR2/$file
done
