#!/usr/bin/env fish

mkdir -p dist
pushd dist
cmake .. -DAMReX_ROOT=/Users/blaschke/Science/amrex/dist \
         -DCMAKE_C_COMPILER=(which gcc-9) \
         -DCMAKE_CXX_COMPILER=(which g++-9) \
         -DCMAKE_Fortran_COMPILER=(which gfortran-9)
make VERBOSE=TRUE -j
popd
