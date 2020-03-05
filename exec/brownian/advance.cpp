
#include "hydro_functions.H"
#include "hydro_functions_F.H"

#include "common_functions.H"

#include "common_namespace.H"

#include "gmres_functions.H"


#include <AMReX_ParallelDescriptor.H>
#include <AMReX_MultiFabUtil.H>

using namespace amrex;
using namespace common;

// argv contains the name of the inputs file entered at the command line
void advance(  std::array< MultiFab, AMREX_SPACEDIM >& umac,
	       MultiFab& pres,
	       const std::array< MultiFab, AMREX_SPACEDIM >& stochMfluxdiv,
	       const std::array< MultiFab, AMREX_SPACEDIM >& sourceTerms,
	       std::array< MultiFab, AMREX_SPACEDIM >& alpha_fc,
	       MultiFab& beta,
               MultiFab& gamma,
	       std::array< MultiFab, NUM_EDGE >& beta_ed,
	       const Geometry geom, const Real& dt)
{

  BL_PROFILE_VAR("advance()",advance);

  Real theta_alpha = 0.;
  Real norm_pre_rhs;

  const BoxArray& ba = beta.boxArray();
  const DistributionMapping& dmap = beta.DistributionMap();

   // rhs_p GMRES solve
   MultiFab gmres_rhs_p(ba, dmap, 1, 0);
   gmres_rhs_p.setVal(0.);

  // rhs_u GMRES solve
  std::array< MultiFab, AMREX_SPACEDIM > gmres_rhs_u;
  for (int d=0; d<AMREX_SPACEDIM; ++d) {
      gmres_rhs_u[d].define(convert(ba,nodal_flag_dir[d]), dmap, 1, 0);
      gmres_rhs_u[d].setVal(0.);
  }

  //////////////////////////////////////////////////
  // ADVANCE velocity field
  //////////////////////////////////////////////////

  // add stochastic forcing to gmres_rhs_u
  for (int d=0; d<AMREX_SPACEDIM; ++d) {
      MultiFab::Add(gmres_rhs_u[d], stochMfluxdiv[d], 0, 0, 1, 0);
      MultiFab::Add(gmres_rhs_u[d], sourceTerms[d], 0, 0, 1, 0);
  }

  // HERE is where you would add the particle forcing to gmres_rhs_u
  //
  //
  //
  
  // call GMRES
  GMRES(gmres_rhs_u,gmres_rhs_p,umac,pres,alpha_fc,beta,beta_ed,gamma,theta_alpha,geom,norm_pre_rhs);

  // fill periodic ghost cells
  for (int d=0; d<AMREX_SPACEDIM; ++d) {
      umac[d].FillBoundary(geom.periodicity());
  }
}
