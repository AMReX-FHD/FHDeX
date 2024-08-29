#include "reactDiff_functions.H"
#include "chemistry_functions.H"

void AdvanceReactionDiffusion(MultiFab& n_old,
                              MultiFab& n_new,
                              const MultiFab& ext_src,
                              const Real& dt,
                              const Real& time,
                              const Geometry& geom) {

    BoxArray ba = n_old.boxArray();
    DistributionMapping dmap = n_old.DistributionMap();

    // store D_Fick on faces
    std::array< MultiFab, AMREX_SPACEDIM > diff_coef_face;
    AMREX_D_TERM(diff_coef_face[0].define(convert(ba,nodal_flag_x), dmap, nspecies, 0);,
                 diff_coef_face[1].define(convert(ba,nodal_flag_y), dmap, nspecies, 0);,
                 diff_coef_face[2].define(convert(ba,nodal_flag_z), dmap, nspecies, 0););

    for (int i=0; i<nspecies; ++i) {
        // load D_fick for species i
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            diff_coef_face[d].setVal(D_Fick[i],i,1,0);
        }
    }

    MultiFab rate1(ba,dmap,nspecies,0);

    if (temporal_integrator == -3) { // multinomial diffusion
        Abort("AdvanceReactionDiffusion() - temporal_integrator=-3 not supported yet");
    }

    Vector<Real> mattingly_lin_comb_coef(2);
    mattingly_lin_comb_coef[0] = 1.;
    mattingly_lin_comb_coef[1] = 0.;

    MultiFab diff_fluxdiv (ba,dmap,nspecies,0);
    MultiFab stoch_fluxdiv(ba,dmap,nspecies,0);

    DiffusiveNFluxdiv(n_old,diff_fluxdiv,diff_coef_face,geom,time);

    if (variance_coef_mass > 0.) {
        Abort("AdvanceReactionDiffusion() - write stochastic case");
    } else {
        stoch_fluxdiv.setVal(0.);
    }

    //!!!!!!!!!!!!!!!
    // time advance !
    //!!!!!!!!!!!!!!!

    if (temporal_integrator == -1) { // forward Euler

        // calculate rates
        // rates could be deterministic or stochastic depending on use_Poisson_rng
        ChemicalRates(n_old,rate1,geom,dt,n_old,mattingly_lin_comb_coef,volume_factor);

        // n_k^{n+1} = n_k^n + dt div (D_k grad n_k)^n
        //                   + dt div (sqrt(2 D_k n_k^n dt) Z) ! Gaussian noise
        //                   + 1/dV * P( f(n_k)*dt*dV )        ! Poisson noise
        //                   + dt ext_src
        MultiFab::LinComb(n_new,1,n_old,0,dt,diff_fluxdiv,0,0,nspecies,0);
        MultiFab::Saxpy(n_new,dt,stoch_fluxdiv,0,0,nspecies,0);
        MultiFab::Saxpy(n_new,dt,rate1,0,0,nspecies,0);
        MultiFab::Saxpy(n_new,dt,ext_src,0,0,nspecies,0);
        n_new.FillBoundary(geom.periodicity());
        MultiFabPhysBC(n_new, geom, 0, nspecies, SPEC_BC_COMP, time);
        
    } else if (temporal_integrator == -2) { // explicit midpoint

        // temporary storage for second rate
        MultiFab rate2(ba,dmap,nspecies,0);

        if (reaction_type == 2) { // explicit midpoint with SSA

            Abort("AdvanceReactionDiffusion() - temporal_integrator=-2 (SSA) not written yet");
            
        } else {

            Abort("AdvanceReactionDiffusion() - temporal_integrator=-2 (non-SSA) not written yet");

        } // explicit midpoint for det/tau/CLE
        
    } else if (temporal_integrator == -4) { // implicit midpoint

        if (reaction_type == 2) { // implicit midpoint with SSA

            Abort("AdvanceReactionDiffusion() - temporal_integrator=-4 (SSA) not written yet");

        } else {

            Abort("AdvanceReactionDiffusion() - temporal_integrator=-4 (non-SSA) not written yet");

        }
    } else {

        Abort("AdvanceReactionDiffusion() - invalid temporal_integrator");

    }
}
