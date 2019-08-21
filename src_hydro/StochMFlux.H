#ifndef _StochMFlux_H
#define _StochMFlux_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Vector.H>

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

  StochMFlux(                 BoxArray,DistributionMapping,Geometry,int);
  void fillMStochastic(       );
  void weightMflux(           Vector< amrex::Real >);
  void multbyVarSqrtEtaTemp(  const amrex::MultiFab&, const std::array< MultiFab, NUM_EDGE >&,
			      const amrex::MultiFab&, const std::array< MultiFab, NUM_EDGE >&,
			      const amrex::Real&);
  void StochMFluxDiv(         std::array< amrex::MultiFab, AMREX_SPACEDIM >&,
                              const int&, const amrex::MultiFab&,
			      const std::array< amrex::MultiFab, NUM_EDGE >&, const amrex::MultiFab&,
			      const std::array< amrex::MultiFab, NUM_EDGE >&,
			      const amrex::Vector< amrex::Real >&, const amrex::Real&);
  void addMfluctuations(      std::array< amrex::MultiFab, AMREX_SPACEDIM >&,
			      const amrex::MultiFab&, const amrex::MultiFab&,
			      const amrex::Real&);
  void addMfluctuations_stag( std::array< amrex::MultiFab, AMREX_SPACEDIM >&,
			      const std::array< amrex::MultiFab, AMREX_SPACEDIM >&,
			      const std::array< amrex::MultiFab, AMREX_SPACEDIM >&,
			      const amrex::Real&);
  void writeMFs(              std::array< MultiFab, AMREX_SPACEDIM >&);

};

#endif
