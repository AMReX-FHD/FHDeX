#!/usr/bin/bash

# required dependencies
module load cray-fftw
module load cmake
module load cudatoolkit

module list

# necessary to use CUDA-Aware MPI and run a job
export CRAY_ACCEL_TARGET=nvidia80

export MPICH_GPU_SUPPORT_ENABLED=1

# optimize CUDA compilation for A100
export AMREX_CUDA_ARCH=8.0

# optimize CPU microarchitecture for AMD EPYC 3rd Gen (Milan/Zen3)
# note: the cc/CC/ftn wrappers below add those
export CXXFLAGS="-march=znver3"
export CFLAGS="-march=znver3"

# compiler environment hints
export CC=cc
export CXX=CC
export FC=ftn
export CUDACXX=$(which nvcc)
export CUDAHOSTCXX=CC

make -j10 USE_CUDA=TRUE USE_HEFFTE_CUFFT=TRUE USE_ASSERTION=TRUE MAX_SPEC=2
