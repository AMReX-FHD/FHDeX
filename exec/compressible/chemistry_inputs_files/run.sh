#!/bin/bash

inputs_file=inputs_NO2_N2O4_eq_SF

mpirun -n 8 ../main3d.gnu.MPI.ex $inputs_file > scrout & 
