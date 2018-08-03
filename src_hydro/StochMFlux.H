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

#if (AMREX_SPACEDIM == 2)
  const int ncomp_cc = 2;
#elif (AMREX_SPACEDIM == 3)
  const int ncomp_cc = 3;
#endif
  const int ncomp_ed = 2;
  
  int n_rngs;
  int stoch_stress_form;

  Geometry geom;
  BoxArray ba;
  DistributionMapping dmap;

  // mflux cell centred
  Vector< MultiFab > mflux_cc;

  // mflux on nodes in 2d
  // mflux on edges in 3d
  Vector<std::array< MultiFab, NUM_EDGE > >  mflux_ed;

public:

  StochMFlux(BoxArray,DistributionMapping,Geometry,int,int);
  void fillMStochastic();

};

StochMFlux::StochMFlux(BoxArray ba_in, DistributionMapping dmap_in, Geometry geom_in, int n_rngs_in, int stoch_stress_form_in=1) {
  n_rngs = n_rngs_in;
  stoch_stress_form = stoch_stress_form_in;

  mflux_cc.resize(n_rngs);
  mflux_ed.resize(n_rngs);
  
  geom = geom_in;
  ba = ba_in;
  dmap = dmap_in;

  // Temporary assumption: filtering_width = 0
  for (int i=0; i<n_rngs; ++i) {
    mflux_cc[i].define(ba, dmap, ncomp_cc, 1);
  }

#if (AMREX_SPACEDIM == 2)
  for (int i=0; i<n_rngs; ++i) {
    mflux_ed[i][0].define(convert(ba,nodal_flag), dmap, ncomp_ed, 0);
  }
#elif (AMREX_SPACEDIM == 3)
  for (int i=0; i<n_rngs; ++i) {
    mflux_ed[i][0].define(convert(ba,nodal_flag_xy), dmap, ncomp_ed, 0);
    mflux_ed[i][1].define(convert(ba,nodal_flag_xz), dmap, ncomp_ed, 0);
    mflux_ed[i][2].define(convert(ba,nodal_flag_yz), dmap, ncomp_ed, 0);
  }
#endif
} 

void StochMFlux::fillMStochastic() {  
  for (int i=0; i<n_rngs; ++i){ 
    switch(stoch_stress_form) {
    case 0: // Non-symmetric
      for (int n=0; n<ncomp_cc; ++n) {
	MultiFABFillRandom(mflux_cc[i],n,1.0,geom);
      }

      for (int d=0; d<NUM_EDGE; ++d) {
	for (int n=0; n<ncomp_ed; ++n) {
	  MultiFABFillRandom(mflux_ed[i][d],n,1.0,geom);
	}
      }
      break;
    default: // Symmetric
      for (int n=0; n<ncomp_cc; ++n) {
	MultiFABFillRandom(mflux_cc[i],n,2.0,geom);
      }

      for (int d=0; d<NUM_EDGE; ++d) {
	MultiFABFillRandom(mflux_ed[i][d],0,1.0,geom);
        MultiFab::Copy(mflux_ed[i][d], mflux_ed[i][d], 0, 1, ncomp_ed-1, 0);
      }
      break;
    }
  }
} 

#endif
