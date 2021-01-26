#!/usr/bin/env fish


# Parse input arguments
argparse "d/delete" -- $argv

if test -n "$_flag_d"
    rm -r dist
end


# Assume location of the AMReX Executable
if ! set -q AMREX_ROOT
    set -gx AMREX_ROOT (realpath  ../amrex)
end

# Default to clang compiler
if ! set -q CC && ! set -q CXX
    set -gx CC (which clang)
    set -gx CXX (which clang++)
end

# Build!
mkdir -p dist
pushd dist
cmake .. -DAMReX_ROOT=$AMREX_ROOT -DCMAKE_INSTALL_PREFIX="."
make VERBOSE=TRUE -j
make install
popd
