#ifndef TAGGING_H
#define TAGGING_H

#include <AMReX_Array4.H>

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
state_error (int i, int j, int k,
             amrex::Array4<char> const& tag,
             amrex::Array4<amrex::Real const> const& state,
             amrex::Real phierr, char tagval)
{
#if 1
//    int hack = 0;
    if (state(i,j,k) < phierr) {
        tag(i,j,k) = tagval;
//        hack = 1;
    }
//    if ( j == 127 ) {
//        amrex::Print() << " tag at " << i << " " << hack << std::endl;
//    }
#else
    // Specification of two boxes that only touch at the corner
    if ( (i > 72 and i < 79 and j > 72 and j < 79) or
         (i > 80 and i < 87 and j > 80 and j < 87) ) {
        tag(i,j,k) = tagval;
    }
#endif
}

#endif
