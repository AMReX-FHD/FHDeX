#!/usr/bin/env fish


# Parse input arguments
argparse "d/delete" -- $argv

if test -n "$_flag_d"
    rm -r amrex
end


# Assume location of the AMReX Executable
if ! set -q AMREX_SOURCE
    set -gx AMREX_SOURCE (realpath  ../../../amrex)
end


# Default to clang compiler
if ! set -q CC && ! set -q CXX
    set -gx CC (which clang)
    set -gx CXX (which clang++)
end


mkdir -p amrex
pushd amrex
cmake -DAMReX_PARTICLES=ON -DAMReX_AMRDATA=ON -DAMReX_EB=ON -DAMReX_OMP=ON \
      -DBUILD_SHARED_LIBS=ON -DAMReX_PIC=ON -DAMReX_FORTRAN_INTERFACES=ON \
      -DAMReX_FORTRAN=ON \
      -DCMAKE_INSTALL_PREFIX="." \
      $AMREX_SOURCE
make -j (math (nproc) / 2) VERBOSE=TRUE
make install
popd
