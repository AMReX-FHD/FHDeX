#!/usr/bin/bash

# required dependencies
module load PrgEnv-gnu
module load craype
module load craype-x86-milan

module list

# optimize CPU microarchitecture for AMD EPYC 3rd Gen (Milan/Zen3)
# note: the cc/CC/ftn wrappers below add those
export CXXFLAGS="-march=znver3"
export CFLAGS="-march=znver3"

# compiler environment hints
export CC=cc
export CXX=CC
export FC=ftn

make -j10 USE_CUDA=FALSE MAX_SPEC=2 USE_ASSERTION=TRUE DEBUG=FALSE
