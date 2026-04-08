#!/usr/bin/bash

## load necessary modules
module load PrgEnv-cray
module load cray-mpich
module load cce

# compiler environment hints
export CC=$(which craycc)
export CXX=$(which craycc)
export FC=$(which crayftn)

make -j10 USE_CUDA=FALSE USE_HIP=FALSE USE_ASSERTION=TRUE
