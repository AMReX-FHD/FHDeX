#ifndef TAGGING_H
#define TAGGING_H

#include <AMReX_Array4.H>

// Tag the cells where phi is below phierr (already scaled to a density):
// the low-density cells become the particle region.
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
state_error (int i, int j, int k,
             amrex::Array4<char> const& tag,
             amrex::Array4<amrex::Real const> const& state,
             amrex::Real phierr, char tagval)
{
    if (state(i,j,k) < phierr) {
        tag(i,j,k) = tagval;
    }
}

#endif
