#!/bin/bash

NRUN=100

rm -f res.cov_var
rm -f res.cov_var2

for ((i=1;i<=$NRUN;i++))
do
  grep "var(theta)" R$i/res.coverage_stat >> res.cov_var
  grep "var2(theta)" R$i/res.coverage_stat >> res.cov_var2
done

python analyze_cov_var.py res.cov_var res.cov_var2 | tee res.cov_var_stat

