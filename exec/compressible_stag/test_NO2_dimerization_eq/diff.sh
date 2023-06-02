#!/bin/bash

DIR1=./
DIR2=../test_NO2_dimerization_eq
FILES="submit_job.sh job_script.sh inputs_NO2_dimerization_eq params_NO2_dimerization_eq.py params_NO2_N2O4_thermochem.py clean.sh diff.sh README"

for file in $FILES
do
  echo "** $file"
  diff $DIR1/$file $DIR2/$file
done
