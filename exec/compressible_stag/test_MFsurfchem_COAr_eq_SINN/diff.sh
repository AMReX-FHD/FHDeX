#!/bin/bash

DIR1=./
DIR2=../test_MFsurfchem_COAr_eq
FILES="run.sh inputs_fhd_stag params_COAr_eq.py zavg.sh zavg.plt ads.sh ads.py ads.plt mass_cons.py mass_cons.sh theta_eq.py theta_eq.sh zavg_rhoY1mean.sh zavg_rhoY1mean.plt clean.sh diff.sh surfcoVar.sh"

for file in $FILES
do
  echo "** $file"
  diff $DIR1/$file $DIR2/$file
done
