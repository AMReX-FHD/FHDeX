#!/usr/bin/bash

## load necessary modules 
module load cmake/3.23.2
module load craype-accel-amd-gfx90a
module load rocm/5.2.0  # waiting for 5.6 for next bump
module load cray-mpich
module load cce/15.0.0  # must be loaded after rocm

# GPU-aware MPI
export MPICH_GPU_SUPPORT_ENABLED=1

# optimize CUDA compilation for MI250X
export AMREX_AMD_ARCH=gfx90a

# compiler environment hints
export CC=$(which hipcc)
export CXX=$(which hipcc)
export FC=$(which ftn)
export CFLAGS="-I${ROCM_PATH}/include"
export CXXFLAGS="-I${ROCM_PATH}/include -Wno-pass-failed"
export LDFLAGS="-L${ROCM_PATH}/lib -lamdhip64 ${PE_MPICH_GTL_DIR_amd_gfx90a} -lmpi_gtl_hsa"

make -j18 USE_CUDA=FALSE USE_HIP=TRUE DO_TURB=TRUE MAX_SPEC=2
