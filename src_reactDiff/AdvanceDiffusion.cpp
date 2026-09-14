#include "reactDiff_functions.H"

// Solves n_t = div ( D grad (n)) + div (sqrt(2*variance*D*n)*W) + g
// where g is a constant in time external source
void AdvanceDiffusion(MultiFab& n_old,
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

    // do not do diffusion if only one cell (well-mixed system)
    // there is no restriction on the number of cells
    // but we can shortcut the single cell case anyway for simplicity
    if (n_cells[0] == 0 && n_cells[1] == 0) {
        Abort("AdvanceDiffusion() - fix one cell case");
    }

    if (reactDiff_diffusion_type == 3) {
        MultinomialDiffusion(n_old,n_new,diff_coef_face,geom,dt,time);
        return;
    }

    MultiFab diff_fluxdiv (ba,dmap,nspecies,0);
    MultiFab stoch_fluxdiv(ba,dmap,nspecies,0);

    DiffusiveNFluxdiv(n_old,diff_fluxdiv,diff_coef_face,geom,time);

    if (variance_coef_mass > 0.) {
        StochasticNFluxdiv(n_old,stoch_fluxdiv,diff_coef_face,geom,dt,time,0);
    } else {
        stoch_fluxdiv.setVal(0.);
    }

    if (reactDiff_diffusion_type == 0 || reactDiff_diffusion_type == 4) {
        // explicit trapezoidal predictor-corrector OR forward Euler

        // forward Euler predictor
        // n_k^{n+1,*} = n_k^n + dt div (D_k grad n_k)^n
        //                     + dt div (sqrt(2 D_k n_k / dt) Z)^n
        //                     + dt ext_src
        MultiFab::Copy(n_new,n_old,0,0,nspecies,0);
        MultiFab::Saxpy(n_new,dt,diff_fluxdiv ,0,0,nspecies,0);
        MultiFab::Saxpy(n_new,dt,stoch_fluxdiv,0,0,nspecies,0);
        MultiFab::Saxpy(n_new,dt,ext_src      ,0,0,nspecies,0);
        n_new.FillBoundary(geom.periodicity());
        MultiFabPhysBC(n_new, geom, 0, nspecies, SPEC_BC_COMP, time);

        if (reactDiff_diffusion_type == 0) {
            /*
              ! Trapezoidal corrector:
              ! n_k^{n+1} = n_k^n + (dt/2) div (D_k grad n_k)^n
              !                   + (dt/2) div (D_k grad n_k)^{n+1,*}
              !                   +  dt    div (sqrt(2 D_k n_k / dt) Z)^n
              !                   +  dt    ext_src
              ! This is the same as stepping to time t+2*dt and then averaging with the state at time t:
              !  n_new = 1/2 * (n_old + n_new + dt*div (D grad n_new) + div (sqrt(2 D_k n_k dt) Z)^n)
              !  which is what we use below
            */

            // compute diffusive flux divergence
            DiffusiveNFluxdiv(n_new,diff_fluxdiv,diff_coef_face,geom,time);

            MultiFab::Saxpy(n_new,1.,n_old,0,0,nspecies,0);
            MultiFab::Saxpy(n_new,dt,diff_fluxdiv ,0,0,nspecies,0);
            MultiFab::Saxpy(n_new,dt,stoch_fluxdiv,0,0,nspecies,0);
            MultiFab::Saxpy(n_new,dt,ext_src      ,0,0,nspecies,0);
            n_new.mult(0.5);
            n_new.FillBoundary(geom.periodicity());
            MultiFabPhysBC(n_new, geom, 0, nspecies, SPEC_BC_COMP, time);
        }

    } else if (reactDiff_diffusion_type == 1) {
        /*
       ! Crank-Nicolson
       ! n_k^{n+1} = n_k^n + (dt/2)(div D_k grad n_k)^n
       !                   + (dt/2)(div D_k grad n_k)^n+1
       !                   +  dt    div (sqrt(2 D_k n_k / dt) Z)^n
       !                   +  dt    ext_src
       !
       ! ( I- (dt/2) div D_k grad) n_k^n+1 = n_k^n
       !                                     + (dt/2)(div D_k grad n_k)^n
       !                                     +  dt    div (sqrt(2 D_k n_k / dt) Z)^n
       !                                     +  dt    ext_src
       ! we combine the entire rhs into stoch_fluxdiv
        */

        MultiFab::Saxpy(stoch_fluxdiv,1.,ext_src,0,0,nspecies,0);
        MultiFab::Saxpy(stoch_fluxdiv,0.5,diff_fluxdiv,0,0,nspecies,0);
        stoch_fluxdiv.mult(dt);
        MultiFab::Saxpy(stoch_fluxdiv,1.,n_old,0,0,nspecies,0);

        ImplicitDiffusion(n_old, n_new, stoch_fluxdiv, diff_coef_face, geom, 0.5*dt, time);

    } else if (reactDiff_diffusion_type == 2) {

        /*
       ! explicit midpoint scheme

       ! n_k^{n+1/2} = n_k^n + (dt/2) div (D_k grad n_k)^n
       !                     + (dt/2) div (sqrt(2 D_k n_k / (dt/2) ) Z_1)^n
       !                     + (dt/2) ext_src
        */

        MultiFab::Copy(n_new,n_old,0,0,nspecies,0);
        MultiFab::Saxpy(n_new,0.5*dt,diff_fluxdiv,0,0,nspecies,0);
        MultiFab::Saxpy(n_new,dt/std::sqrt(2.),stoch_fluxdiv,0,0,nspecies,0);
        MultiFab::Saxpy(n_new,0.5*dt,ext_src,0,0,nspecies,0);
        n_new.FillBoundary(geom.periodicity());
        MultiFabPhysBC(n_new, geom, 0, nspecies, SPEC_BC_COMP, time);

        // compute diffusive flux divergence at t^{n+1/2}
        DiffusiveNFluxdiv(n_new,diff_fluxdiv,diff_coef_face,geom,time);

        if (variance_coef_mass > 0.) {
            GenerateStochasticFluxdivCorrector(n_old,n_new,stoch_fluxdiv,diff_coef_face,dt,time,geom);
        }

       /*
       ! n_k^{n+1} = n_k^n + dt div (D_k grad n_k)^{n+1/2}
       !                   + dt div (sqrt(2 D_k n_k^n dt) Z_1 / sqrt(2) )
       !                   + dt div (sqrt(2 D_k n_k^? dt) Z_2 / sqrt(2) )
       !                   + dt ext_src
       ! where
       ! n_k^? = n_k^n               (midpoint_stoch_flux_type=1)
       !       = n_k^pred            (midpoint_stoch_flux_type=2)
       !       = 2*n_k^pred - n_k^n  (midpoint_stoch_flux_type=3)
       */

        MultiFab::Copy(n_new,n_old,0,0,nspecies,0);
        MultiFab::Saxpy(n_new,dt,diff_fluxdiv,0,0,nspecies,0);
        MultiFab::Saxpy(n_new,dt/std::sqrt(2.),stoch_fluxdiv,0,0,nspecies,0);
        MultiFab::Saxpy(n_new,dt,ext_src,0,0,nspecies,0);
        n_new.FillBoundary(geom.periodicity());
        MultiFabPhysBC(n_new, geom, 0, nspecies, SPEC_BC_COMP, time);

    } else {
        Abort("AdvanceDiffusion() - invalid reactDiff_diffusion_type");
    }

}

void GenerateStochasticFluxdivCorrector(MultiFab& n_old,
                                        MultiFab& n_new,
                                        MultiFab& stoch_fluxdiv,
                                        const std::array< MultiFab, AMREX_SPACEDIM >& diff_coef_face,
                                        const Real& dt,
                                        const Real& time,
                                        const Geometry& geom) {

    // fill random flux multifabs with new random numbers and
    // compute second-stage stochastic flux divergence and
    // add to first-stage stochastic flux divergence
    if (midpoint_stoch_flux_type == 1) {

        // use n_old
        StochasticNFluxdiv(n_old,stoch_fluxdiv,diff_coef_face,geom,dt,time,1);

    } else if (midpoint_stoch_flux_type == 2) {

        // use n_pred
        StochasticNFluxdiv(n_new,stoch_fluxdiv,diff_coef_face,geom,dt,time,1);

    } else if (midpoint_stoch_flux_type == 3) {

        // We use n_new=2*n_pred-n_old here as temporary storage since we will overwrite it shortly
        n_new.mult(2.);
        MultiFab::Subtract(n_new,n_old,0,0,nspecies,1);

        // use n_new=2*n_pred-n_old
        StochasticNFluxdiv(n_new,stoch_fluxdiv,diff_coef_face,geom,dt,time,1);

    } else {
        Abort("GenerateStochasticFluxdivCorrector() - invalid midpoint_stoch_flux_type");
    }
}
