#!/bin/bash

inputs_file=inputs_fhd_stag

mpirun -n 8 ../main3d.gnu.MPI.ex $inputs_file > scrout &
