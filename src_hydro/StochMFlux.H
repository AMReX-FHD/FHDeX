#ifndef _StochMFlux_H
#define _StochMFlux_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Vector.H>

#include <AMReX_MultiFabUtil.H>
#include <AMReX_VisMF.H>

#include <common_namespace.H>

using namespace amrex;
using namespace common;

class StochMFlux {

  const int ncomp_ed = 2;
  int n_rngs;
  
  Geometry geom;

  // Vectors containing mflux Gaussian white noise field at each stage:
  // mflux cell-centered
  Vector< MultiFab > mflux_cc;
  // mflux on nodes in 2d, on edges in 3d
  Vector<std::array< MultiFab, NUM_EDGE > >  mflux_ed;

  // Weighted mflux cell-centered & nodal ...
  // also, multiplied by variance*sqrt(eta*temperature):
  MultiFab mflux_cc_weighted;
  std::array< MultiFab, NUM_EDGE >  mflux_ed_weighted;

public:

  StochMFlux(BoxArray,DistributionMapping,Geometry,int);
  void fillMStochastic();
  void weightMflux(Vector< amrex::Real >);
  void multbyVarSqrtEtaTemp(const MultiFab& eta_cc,  
			    const std::array< MultiFab, NUM_EDGE >& eta_ed, 
			    const MultiFab& temp_cc, 
			    const std::array< MultiFab, NUM_EDGE >& temp_ed,
			    const amrex::Real& dt);
  void stochMforce(std::array< MultiFab, AMREX_SPACEDIM >& mfluxdiv, 
		   const MultiFab& eta_cc,  
		   const std::array< MultiFab, NUM_EDGE >& eta_ed, 
		   const MultiFab& temp_cc, 
		   const std::array< MultiFab, NUM_EDGE >& temp_ed, 
		   const Vector< amrex::Real >& weights,
		   const amrex::Real& dt);
  void writeMFs(std::array< MultiFab, AMREX_SPACEDIM >& mfluxdiv);

};

#endif

