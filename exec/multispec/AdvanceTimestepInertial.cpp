
#include "hydro_functions.H"

#include "common_functions.H"

#include "gmres_functions.H"

#include "multispec_functions.H"


#include <AMReX_ParallelDescriptor.H>
#include <AMReX_MultiFabUtil.H>


// argv contains the name of the inputs file entered at the command line
void AdvanceTimestepInertial(std::array< MultiFab, AMREX_SPACEDIM >& umac,
                             MultiFab& rho_old, MultiFab& rho_new,
                             MultiFab& rhotot_old, MultiFab& rhotot_new,
                             MultiFab& pi, const MultiFab& eta, 
                             const std::array< MultiFab, NUM_EDGE >&  eta_ed,
                             const MultiFab& kappa, const MultiFab& Temp,
                             const std::array< MultiFab, NUM_EDGE >& Temp_ed,
                             MultiFab& diff_mass_fluxdiv,
                             MultiFab& stoch_mass_fluxdiv_mass_fluxdiv,
                             std::array< MultiFab, AMREX_SPACEDIM >& stoch_mass_flux,
                             const Real& dt, const Real& time, const int& istep,
                             const Geometry& geom)
{
  
  BL_PROFILE_VAR("AdvanceTimestepInertial()",AdvanceTimestepInertial);

  Real theta_alpha = 1./dt;

}
