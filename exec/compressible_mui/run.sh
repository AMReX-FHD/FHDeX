#!/bin/bash

mpirun -np 1 ./main3d.gnu.DEBUG.MPI.ex inputs_equil_3d : -np 1 ./pseudo_kmc.bin
