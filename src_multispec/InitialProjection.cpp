#include "multispec_functions.H"
#include "gmres_functions.H"

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
                       const Real& dt, const Real& time, const Geometry& geom,
                       MultiFab& charge_old,
                       std::array<MultiFab,AMREX_SPACEDIM>& grad_Epot_old,
                       MultiFab& Epot,
                       MultiFab& permittivity)
{
    BL_PROFILE_VAR("InitialProjection()",InitialProjection);

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
        // predictor integrates over full time step
        dt_eff = dt;
    }

    /*
    if (nreactions > 0) {
        if (algorithm_type != 5) {
            Abort("Error: only algorithm_type=5 allowed for nreactions>0");
        } else if (use_Poisson_rng == 2) {
            Abort("Error: currently use_Poisson_rng=2 not allowed for algorithm_type=5 and nreactions>0");
        }
    }
    */

    MultiFab mac_rhs(ba,dmap,1,0);
    MultiFab divu   (ba,dmap,1,0);
    MultiFab phi    (ba,dmap,1,1);

    std::array< MultiFab, AMREX_SPACEDIM > rhotot_fc;
    std::array< MultiFab, AMREX_SPACEDIM > rhototinv_fc;
    std::array< MultiFab, AMREX_SPACEDIM > diff_mass_flux;
    std::array< MultiFab, AMREX_SPACEDIM > total_mass_flux;
    std::array< MultiFab, AMREX_SPACEDIM > gradp;
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        rhotot_fc[d]      .define(convert(ba,nodal_flag_dir[d]), dmap,        1, 0);
        rhototinv_fc[d]   .define(convert(ba,nodal_flag_dir[d]), dmap,        1, 0);
        diff_mass_flux[d] .define(convert(ba,nodal_flag_dir[d]), dmap, nspecies, 0);
        total_mass_flux[d].define(convert(ba,nodal_flag_dir[d]), dmap, nspecies, 0);
        gradp[d]          .define(convert(ba,nodal_flag_dir[d]), dmap,        1, 0);
    }

    // set inhomogeneous velocity bc's to values supplied in inhomogeneous_bc_val
    //
    //
    //

    // don't fill random numbers for Boussinesq (algorithm_type=6)
    // since rhobars are equal the RHS of constraint is zero regardless
    // by not filling we preserve the random number sequences for regression purposes
    if (variance_coef_mass != 0. && algorithm_type != 6) {
        sMassFlux.fillMassStochastic();
    }

    ComputeMassFluxdiv(rho,rhotot,Temp,diff_mass_fluxdiv,stoch_mass_fluxdiv,
                       diff_mass_flux,stoch_mass_flux,sMassFlux,dt_eff,time,geom,weights,
                       charge_old,grad_Epot_old,Epot,permittivity);

    // assumble total fluxes to be used in reservoirs
    //
    //
    //

    // compute chemical rates m_i*R_i
    //
    //
    //

    // set the Dirichlet velocity value on reservoir faces
    // call reservoir_bc_fill
    //
    //

    if (restart < 0) {

        // project the velocities
        // only for non-restarting runs
        mac_rhs.setVal(0.);

        // set mac_rhs to -S = sum_i div(F_i)/rhobar_i
        for (int i=0; i<nspecies; ++i) {
            MultiFab::Saxpy(mac_rhs,-1./rhobar[i],diff_mass_fluxdiv,i,0,1,0);
            if (variance_coef_mass != 0.) {
                MultiFab::Saxpy(mac_rhs,-1./rhobar[i],stoch_mass_fluxdiv,i,0,1,0);
            }
            // if nreactions>0, also add sum_i -(m_i*R_i)/rhobar_i
        }

        // build rhs = div(v^init) - S^0

        for (int i=0; i<AMREX_SPACEDIM; ++i) {
            // to deal with reservoirs
            // set normal velocity on physical domain boundaries
            MultiFabPhysBCDomainVel(umac[i],geom,i);
            // fill periodic and interior ghost cells
            umac[i].FillBoundary(geom.periodicity());
        }

        // set divu = div(v^init)
        ComputeDiv(divu,umac,0,0,1,geom,0);

        // add div(v^init) to mac_rhs
        // now mac_rhs = div(v^init - S)
        MultiFab::Add(mac_rhs,divu,0,0,1,0);

        // average rhotot on faces
        AverageCCToFace(rhotot,rhotot_fc,0,1,RHO_BC_COMP,geom);

        // compute (1/rhotot on faces)
        for (int i=0; i<AMREX_SPACEDIM; ++i) {
            rhototinv_fc[i].setVal(1.);
            MultiFab::Divide(rhototinv_fc[i],rhotot_fc[i],0,0,1,0);
        }

        // solve div (1/rhotot) grad phi = div(v^init) - S^0
        // solve to completion, i.e., use the 'full' solver
        phi.setVal(0.);

        MacProj macproj;
        macproj.Define(ba,dmap,geom);
        macproj.Solve(rhototinv_fc,mac_rhs,phi,geom,true);

        // v^0 = v^init - (1/rho^0) grad phi
        SubtractWeightedGradP(umac,rhototinv_fc,phi,gradp,geom);

        // fill ghost cells
        for (int i=0; i<AMREX_SPACEDIM; ++i) {
            // set normal velocity of physical domain boundaries
            MultiFabPhysBCDomainVel(umac[i],geom,i);
            // set transverse velocity behind physical boundaries
            int is_inhomogeneous = 1;
            MultiFabPhysBCMacVel(umac[i],geom,i,is_inhomogeneous);
            // fill periodic and interior ghost cells
            umac[i].FillBoundary(geom.periodicity());
        }
    }
}
