#!/bin/bash

COMPDIR=../RUN
#COMPDIR=../RUN_COadsPt111_coverage

FILES=`ls`

for file in $FILES
do
  echo "****** $file ******"
  diff $file $COMPDIR/$file
  echo
done
