#!/bin/bash

echo "mpirun -n 4 ./prog1.ex : -n 2 ./prog1.ex"
mpirun -n 4 ./prog1.ex : -n 2 ./prog1.ex
