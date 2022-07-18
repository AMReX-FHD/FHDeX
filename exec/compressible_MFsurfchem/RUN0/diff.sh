#!/bin/bash

DIR1=./
DIR2=../RUN0
FILES="params_CO_Ar_eq_SF.py inputs_fhd run.sh zavg.sh zavg.plt clean.sh"

for file in $FILES
do
  echo "** $file"
  diff $DIR1/$file $DIR2/$file
done
