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
        StochasticNFluxdiv(n_old,stoch_fluxdiv,diff_coef_face,geom,dt,time,0);
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
            
        } else { // explicit midpoint for det/tau/CLE

            //!!!!!!!!!!!!!!
            // predictor   !
            //!!!!!!!!!!!!!!

            // calculate rates from a(n_old)
            ChemicalRates(n_old,rate1,geom,0.5*dt,n_old,mattingly_lin_comb_coef,volume_factor);

            /*
         ! n_k^{n+1/2} = n_k^n + (dt/2)       div (D_k grad n_k)^n
         !                     + (dt/sqrt(2)) div sqrt(2 D_k n_k^n / (dt*dV)) Z_1 ! Gaussian noise
         !                     + 1/dV * P_1( f(n_k)*(dt/2)*dV )                   ! Poisson noise
         !                     + (dt/2)        ext_src
            */
            MultiFab::Copy(n_new,n_old,0,0,nspecies,0);
            MultiFab::Saxpy(n_new,0.5*dt,diff_fluxdiv,0,0,nspecies,0);
            MultiFab::Saxpy(n_new,dt/std::sqrt(2.),stoch_fluxdiv,0,0,nspecies,0);
            MultiFab::Saxpy(n_new,0.5*dt,rate1,0,0,nspecies,0);
            MultiFab::Saxpy(n_new,0.5*dt,ext_src,0,0,nspecies,0);
            n_new.FillBoundary(geom.periodicity());
            MultiFabPhysBC(n_new, geom, 0, nspecies, SPEC_BC_COMP, time);

            //!!!!!!!!!!!!!!
            // corrector   !
            //!!!!!!!!!!!!!!

            // Here we do not write this in the form that Mattingly et al do
            // where we just continue the second half of the time step from where we left
            // Rather, we compute terms at the midpoint and then add contributions from both 
            // halves of the time step to n_old
            // This works simpler with diffusion but we have to store both rates1 and rates2
         
            // compute diffusive flux divergence
            DiffusiveNFluxdiv(n_new,diff_fluxdiv,diff_coef_face,geom,time);

            // calculate rates from 2*a(n_pred)-a(n_old)
            mattingly_lin_comb_coef[0] = -1.;
            mattingly_lin_comb_coef[1] = 2.;
            ChemicalRates(n_old,rate2,geom,0.5*dt,n_new,mattingly_lin_comb_coef,volume_factor);

            //compute stochastic flux divergence and add to the ones from the predictor stage
            if (variance_coef_mass > 0.) {
                GenerateStochasticFluxdivCorrector(n_old,n_new,stoch_fluxdiv,diff_coef_face,dt,time,geom);
            }

            /*
         ! n_k^{n+1} = n_k^n + dt div (D_k grad n_k)^{n+1/2}
         !                   + dt div (sqrt(2 D_k n_k^n / (dt*dV)) Z_1 / sqrt(2) ) ! Gaussian noise
         !                   + dt div (sqrt(2 D_k n_k^? / (dt*dV)) Z_2 / sqrt(2) ) ! Gaussian noise
         !                   + 1/dV * P_1( f(n_k)*(dt/2)*dV )                        ! Poisson noise
         !                   + 1/dV * P_2( (2*f(n_k^pred)-f(n_k))*(dt/2)*dV )        ! Poisson noise
         !                   + dt ext_src
         ! where
         ! n_k^? = n_k^n               (midpoint_stoch_flux_type=1)
         !       = n_k^pred            (midpoint_stoch_flux_type=2)
         !       = 2*n_k^pred - n_k^n  (midpoint_stoch_flux_type=3)
            */

            MultiFab::Copy(n_new,n_old,0,0,nspecies,0);
            MultiFab::Saxpy(n_new,dt,diff_fluxdiv,0,0,nspecies,0);
            MultiFab::Saxpy(n_new,dt/std::sqrt(2.),stoch_fluxdiv,0,0,nspecies,0);
            MultiFab::Saxpy(n_new,0.5*dt,rate1,0,0,nspecies,0);
            MultiFab::Saxpy(n_new,0.5*dt,rate2,0,0,nspecies,0);
            MultiFab::Saxpy(n_new,dt,ext_src,0,0,nspecies,0);
            n_new.FillBoundary(geom.periodicity());
            MultiFabPhysBC(n_new, geom, 0, nspecies, SPEC_BC_COMP, time);
                                     
        } // explicit midpoint for det/tau/CLE
        
    } else if (temporal_integrator == -4) { // implicit midpoint

        if (reaction_type == 2) { // implicit midpoint with SSA

            Abort("AdvanceReactionDiffusion() - temporal_integrator=-4 (SSA) not written yet");

        } else { // implicit midpoint for det/tau/CLE

/*
         ! backward Euler predictor to half-time
         ! n_k^{n+1/2} = n_k^n + (dt/2)       div (D_k grad n_k)^{n+1/2}
         !                     + (dt/sqrt(2)) div sqrt(2 D_k n_k^n / (dt*dV)) Z_1 ! Gaussian noise
         !                     + 1/dV * P_1( f(n_k)*(dt/2)*dV )                   ! Poisson noise
         !                     + (dt/2)       ext_src
         !
         ! in delta form
         !
         ! (I - div (dt/2) D_k grad) delta n_k =   (dt/2)       div (D_k grad n_k^n)
         !                                       + (dt/sqrt(2)) div (sqrt(2 D_k n_k^n / (dt*dV)) Z_1
         !                                       + 1/dV * P_1( f(n_k)*(dt/2)*dV )
         !                                       + (dt/2) ext_src
         
         do n=1,nlevs
            call multifab_build(rhs(n),mla%la(n),nspecies,0)
            call multifab_build(rate2(n),mla%la(n),nspecies,0)
         end do

         ! calculate rates
         ! rates could be deterministic or stochastic depending on use_Poisson_rng
         call chemical_rates(mla,n_old,rate1,dx,dt/2.d0,vol_fac_in=volume_factor)

         do n=1,nlevs
            call multifab_setval(rhs(n),0.d0)
            call multifab_saxpy_3(rhs(n),dt/2.d0,diff_fluxdiv(n))
            call multifab_saxpy_3(rhs(n),dt/sqrt(2.d0),stoch_fluxdiv(n))
            call multifab_saxpy_3(rhs(n),dt/2.d0,rate1(n))
            if(present(ext_src)) call multifab_saxpy_3(rhs(n),dt/2.d0,ext_src(n))
         end do

         call implicit_diffusion(mla,n_old,n_new,rhs,diff_coef_face,dx,dt,the_bc_tower)

         ! corrector

         ! calculate rates from 2*a(n_pred)-a(n_old)
         call chemical_rates(mla,n_old,rate2,dx,dt/2.d0,n_new,mattingly_lin_comb_coef,vol_fac_in=volume_factor)

         ! compute stochastic flux divergence and add to the ones from the predictor stage
         if (variance_coef_mass .gt. 0.d0) then
            
            ! first, fill random flux multifabs with new random numbers
            call fill_mass_stochastic(mla,the_bc_tower%bc_tower_array)

            ! compute n on faces to use in the stochastic flux in the corrector
            ! three possibilities
            call generate_stochastic_fluxdiv_corrector()

         end if

         ! Crank-Nicolson
         ! n_k^{n+1} = n_k^n + (dt/2) div (D_k grad n_k)^n
         !                   + (dt/2) div (D_k grad n_k)^{n+1}
         !                   + dt div (sqrt(2 D_k n_k^n / (dt*dV)) Z_1 / sqrt(2) ) ! Gaussian noise
         !                   + dt div (sqrt(2 D_k n_k^? / (dt*dV)) Z_2 / sqrt(2) ) ! Gaussian noise
         !                   + 1/dV * P_1( f(n_k)*(dt/2)*dV )                        ! Poisson noise
         !                   + 1/dV * P_2( (2*f(n_k^pred)-f(n_k))*(dt/2)*dV )        ! Poisson noise
         !                   + dt ext_src
         !
         ! in delta form
         !
         ! (I - div (dt/2) D_k grad) delta n_k =   dt div (D_k grad n_k^n)
         !                   + dt div (sqrt(2 D_k n_k^n / (dt*dV)) Z_1 / sqrt(2) ) ! Gaussian noise
         !                   + dt div (sqrt(2 D_k n_k^? / (dt*dV)) Z_2 / sqrt(2) ) ! Gaussian noise
         !                   + 1/dV * P_1( f(n_k)*(dt/2)*dV )                        ! Poisson noise
         !                   + 1/dV * P_2( (2*f(n_k^pred)-f(n_k))*(dt/2)*dV )        ! Poisson noise
         !                   + dt ext_src

         do n=1,nlevs
            call multifab_setval(rhs(n),0.d0)
            call multifab_saxpy_3(rhs(n),dt,diff_fluxdiv(n))
            call multifab_saxpy_3(rhs(n),dt/sqrt(2.d0),stoch_fluxdiv(n))
            call multifab_saxpy_3(rhs(n),dt/2.d0,rate1(n))
            call multifab_saxpy_3(rhs(n),dt/2.d0,rate2(n))
            if(present(ext_src)) call multifab_saxpy_3(rhs(n),dt,ext_src(n))
         end do
      
         call implicit_diffusion(mla,n_old,n_new,rhs,diff_coef_face,dx,dt,the_bc_tower)
*/
             Abort("AdvanceReactionDiffusion() - temporal_integrator=-4 (non-SSA) not written yet");

        }
    } else {

        Abort("AdvanceReactionDiffusion() - invalid temporal_integrator");

    }
}
