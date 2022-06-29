#!/bin/bash

DIR1=./
DIR2=../test_COadsPt
FILES="inputs_fhd_stag run.sh zavg.sh zavg.plt ads.sh ads.py ads.plt clean.sh"

for file in $FILES
do
  echo "** $file"
  diff $DIR1/$file $DIR2/$file
done
