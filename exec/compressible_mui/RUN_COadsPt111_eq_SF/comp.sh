#!/bin/bash

COMPDIR=../RUN

FILES=`ls`

for file in $FILES
do
  echo "****** $file ******"
  diff $file $COMPDIR/$file
  echo
done
