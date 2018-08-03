#ifndef _StochMFlux_H
#define _StochMFlux_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Vector.H>

#include <AMReX_MultiFabUtil.H>

#include <common_namespace.H>

using namespace amrex;
using namespace common;

class StochMFlux {
    
int n_rngs;

// mflux cell centred
MultiFab mflux;

// mflux on nodes in 2d
// mflux on edges in 3d
Vector<std::array< MultiFab, NUM_EDGE > >  mflux_ed;

public:

StochMFlux(int,BoxArray,DistributionMapping);
// void fun() {}

};

StochMFlux::StochMFlux(int nrngs, BoxArray ba, DistributionMapping dmap) {  
    n_rngs = nrngs;
    mflux_ed.resize(n_rngs);

    mflux.define(ba, dmap, 1, std::max(1,filtering_width));

#if (AMREX_SPACEDIM == 2)
    for (int i=0; i<n_rngs; ++i) {
      mflux_ed[i][0].define(convert(ba,nodal_flag), dmap, 2, filtering_width);
    }
#elif (AMREX_SPACEDIM == 3)
    for (int i=0; i<n_rngs; ++i) {
      mflux_ed[i][0].define(convert(ba,nodal_flag_xy), dmap, 2, filtering_width);
      mflux_ed[i][1].define(convert(ba,nodal_flag_xz), dmap, 2, filtering_width);
      mflux_ed[i][2].define(convert(ba,nodal_flag_yz), dmap, 2, filtering_width);
    }
#endif
} 

#endif
