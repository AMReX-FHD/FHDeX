#!/usr/bin/env fish

mkdir -p dist
pushd dist
cmake .. -DAMReX_ROOT=/Users/johannesblaschke/Science/amrex/dist \
         -DCMAKE_INSTALL_PREFIX=. \
         -DCMAKE_C_COMPILER=(which gcc-9) \
         -DCMAKE_CXX_COMPILER=(which g++-9) \
         -DCMAKE_Fortran_COMPILER=(which gfortran-9)
make VERBOSE=TRUE -j
make install
popd
