
#include "hydro_functions.H"
#include "hydro_functions_F.H"

#include "common_functions.H"
#include "common_functions_F.H"

#include "common_namespace.H"

#include "gmres_functions.H"
#include "gmres_functions_F.H"

#include "gmres_namespace.H"

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_MultiFabUtil.H>

using namespace amrex;
using namespace common;
using namespace gmres;

// argv contains the name of the inputs file entered at the command line
void advance(  std::array< MultiFab, AMREX_SPACEDIM >& umac,
	       MultiFab& pres,
	       const std::array< MultiFab, AMREX_SPACEDIM >& stochMfluxdiv,
	       std::array< MultiFab, AMREX_SPACEDIM >& alpha_fc,
	       MultiFab& beta,
               MultiFab& gamma,
	       std::array< MultiFab, NUM_EDGE >& beta_ed,
	       const Geometry geom, const Real& dt)
{

  BL_PROFILE_VAR("advance()",advance);

  const Real* dx = geom.CellSize();
  const Real dtinv = 1.0/dt;
  Real theta_alpha = 0.;
  Real norm_pre_rhs;

  const BoxArray& ba = beta.boxArray();
  const DistributionMapping& dmap = beta.DistributionMap();

   // rhs_p GMRES solve
   MultiFab gmres_rhs_p(ba, dmap, 1, 0);
   gmres_rhs_p.setVal(0.);

  // rhs_u GMRES solve
  std::array< MultiFab, AMREX_SPACEDIM > gmres_rhs_u;
  AMREX_D_TERM(gmres_rhs_u[0].define(convert(ba,nodal_flag_x), dmap, 1, 0);,
	       gmres_rhs_u[1].define(convert(ba,nodal_flag_y), dmap, 1, 0);,
	       gmres_rhs_u[2].define(convert(ba,nodal_flag_z), dmap, 1, 0););
  AMREX_D_TERM(gmres_rhs_u[0].setVal(0.);,
	       gmres_rhs_u[1].setVal(0.);,
	       gmres_rhs_u[2].setVal(0.););

  //////////////////////////////////////////////////
  // ADVANCE velocity field
  //////////////////////////////////////////////////

  // add stochastic forcing to gmres_rhs_u
  AMREX_D_TERM(MultiFab::Add(gmres_rhs_u[0], stochMfluxdiv[0], 0, 0, 1, 0);,
	       MultiFab::Add(gmres_rhs_u[1], stochMfluxdiv[1], 0, 0, 1, 0);,
	       MultiFab::Add(gmres_rhs_u[2], stochMfluxdiv[2], 0, 0, 1, 0););

  // HERE is where you would add the particle forcing to gmres_rhs_u
  //
  //
  //
  
  AMREX_D_TERM(umac[0].setVal(0.);,
	       umac[1].setVal(0.);,
	       umac[2].setVal(0.););
  pres.setVal(0.);

  
  // call GMRES
  GMRES(gmres_rhs_u,gmres_rhs_p,umac,pres,alpha_fc,beta,beta_ed,gamma,theta_alpha,geom,norm_pre_rhs);

  // fill periodic ghost cells
  AMREX_D_TERM(umac[0].FillBoundary(geom.periodicity());,
	       umac[1].FillBoundary(geom.periodicity());,
	       umac[2].FillBoundary(geom.periodicity()););

  //////////////////////////////////////////////////

}
