#include "multispec_functions.H"

/*

this routine is only called for all inertial simulations (both restart and non-restart)
it does the following:
  1. fill mass random numbers
  2. computes mass fluxes and flux divergences
     if restarting, the subroutine ends; otherwise
  3. perform an initial projection
  
overdamped schemes need to do 1. and 2. within the advance_timestep routine
in principle, performing an initial projection for overdamped will change
the reference state for the GMRES solver
For overdamped the first ever solve cannot have a good reference state
so in general there is the danger it will be less accurate than subsequent solves
but I do not see how one can avoid that
From this perspective it may be useful to keep initial_projection even in overdamped
because different gmres tolerances may be needed in the first step than in the rest

*/
void InitialProjection(std::array< MultiFab, AMREX_SPACEDIM >& umac,
                       MultiFab& rho, MultiFab& rhotot,
                       MultiFab& diff_mass_fluxdiv,
                       MultiFab& stoch_mass_fluxdiv,
                       std::array< MultiFab, AMREX_SPACEDIM >& stoch_mass_flux,
                       StochMassFlux& sMassFlux,
                       const MultiFab& Temp, const MultiFab& eta,
                       const std::array< MultiFab, NUM_EDGE >& eta_ed,
                       const Real& dt, const Real& time, const Geometry& geom)
{

    if (algorithm_type == 2) {
        Abort("InitialProjection.cpp: should not call initial_projection for overdamped scheme");
    }

    BoxArray ba = rho.boxArray();
    DistributionMapping dmap = rho.DistributionMap();
    
    Real dt_eff;
    
    Vector<Real> weights;
    if (algorithm_type == 5 || algorithm_type == 6) {
        weights = {1., 0.};
        // for midpoint scheme where predictor goes to t^{n+1/2}
        dt_eff = 0.5*dt;
    }
    else {
        weights = {1.};
        dt_eff = dt;
    }

    MultiFab mac_rhs(ba,dmap,1,0);
    MultiFab divu   (ba,dmap,1,0);
    MultiFab phi    (ba,dmap,1,1);

    std::array< MultiFab, AMREX_SPACEDIM > rhotot_fc;
    std::array< MultiFab, AMREX_SPACEDIM > rhototinv_fc;
    std::array< MultiFab, AMREX_SPACEDIM > diff_mass_flux;
    std::array< MultiFab, AMREX_SPACEDIM > total_mass_flux;
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        rhotot_fc[d]      .define(convert(ba,nodal_flag_dir[d]), dmap,        1, 0);
        rhototinv_fc[d]   .define(convert(ba,nodal_flag_dir[d]), dmap,        1, 0);
        diff_mass_flux[d] .define(convert(ba,nodal_flag_dir[d]), dmap, nspecies, 0);
        total_mass_flux[d].define(convert(ba,nodal_flag_dir[d]), dmap, nspecies, 0);
    }
    
    if (variance_coef_mass != 0.) {
        sMassFlux.fillMassStochastic();
    }
        
    ComputeMassFluxdiv(rho,rhotot,diff_mass_fluxdiv,diff_mass_flux,dt,time,geom);







    // set inhomogeneous velocity bc's to values supplied in inhomogeneous_bc_val
    //
    //
    //




}
