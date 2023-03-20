#!/bin/bash

DIR1=./
DIR2=../test_MFsurfchem_depletion
FILES="ads.sh clean.sh diff.sh inputs_fhd_stag run.sh zavg.sh"

for file in $FILES
do
  echo "** $file"
  diff $DIR1/$file $DIR2/$file
done
