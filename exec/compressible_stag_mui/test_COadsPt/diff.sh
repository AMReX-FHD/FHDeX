#!/bin/bash

DIR1=./
DIR2=../test_COadsPt
FILES="clean.sh diff.sh in.kmc inputs_fhd run.sh"

for file in $FILES
do
  echo "** $file"
  diff $DIR1/$file $DIR2/$file
done
