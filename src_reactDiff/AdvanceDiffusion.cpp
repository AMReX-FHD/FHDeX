#include "reactDiff_functions.H"

#include "AMReX_MLMG.H"
#include <AMReX_MLABecLaplacian.H>

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

    // store for one component of D_Fick
    std::array< MultiFab, AMREX_SPACEDIM > diff_coef_face;
    AMREX_D_TERM(diff_coef_face[0].define(convert(ba,nodal_flag_x), dmap, 1, 0);,
                 diff_coef_face[1].define(convert(ba,nodal_flag_y), dmap, 1, 0);,
                 diff_coef_face[2].define(convert(ba,nodal_flag_z), dmap, 1, 0););

    // do not do diffusion if only one cell (well-mixed system)
    // there is no restriction on the number of cells
    // but we can shortcut the single cell case anyway for simplicity
    if (n_cells[0] == 0 && n_cells[1] == 0) {
        Abort("AdvanceDiffusion() - fix one cell case");
    }

    if (reactDiff_diffusion_type == 3) {
        Abort("AdvanceDiffusion() - write multinomial case");
        return;
    }

    MultiFab diff_fluxdiv (ba,dmap,nspecies,0);
    MultiFab stoch_fluxdiv(ba,dmap,nspecies,0);

    Abort("Write DiffusiveNFluxdiv()");
    // DiffusiveNFluxdiv();

    if (variance_coef_mass > 0.) {
        Abort("AdvanceDiffusion() - write stochastic case");
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

        if (reactDiff_diffusion_type == 4) {
            Abort("AdvanceDiffusion() - write trapezoidal corrector");

            /*

          ! Trapezoidal corrector:
          ! n_k^{n+1} = n_k^n + (dt/2) div (D_k grad n_k)^n
          !                   + (dt/2) div (D_k grad n_k)^{n+1,*}
          !                   +  dt    div (sqrt(2 D_k n_k / dt) Z)^n
          !                   +  dt    ext_src
          ! This is the same as stepping to time t+2*dt and then averaging with the state at time t:
          !  n_new = 1/2 * (n_old + n_new + dt*div (D grad n_new) + div (sqrt(2 D_k n_k dt) Z)^n)
          !  which is what we use below

          ! compute diffusive flux divergence
          call diffusive_n_fluxdiv(mla,n_new,diff_coef_face,diff_fluxdiv,dx,the_bc_tower)

          do n=1,nlevs
             call multifab_plus_plus_c(n_new(n),1,n_old(n),1,nspecies,0)
             call multifab_saxpy_3(n_new(n),dt,diff_fluxdiv(n))
             call multifab_saxpy_3(n_new(n),dt,stoch_fluxdiv(n))
             call multifab_saxpy_3(n_new(n),dt,ext_src(n))
             call multifab_mult_mult_s_c(n_new(n),1,0.5d0,nspecies,0)
             call multifab_fill_boundary(n_new(n))
             call multifab_physbc(n_new(n),1,scal_bc_comp,nspecies, &
                                  the_bc_tower%bc_tower_array(n),dx_in=dx(n,:))
          end do

            */
        }
        
    } else if (reactDiff_diffusion_type == 1) {
        Abort("AdvanceDiffusion() - write Crank-Nicolson");

        /*
       ! Crank-Nicolson
       ! n_k^{n+1} = n_k^n + (dt/2)(div D_k grad n_k)^n
       !                   + (dt/2)(div D_k grad n_k)^n+1
       !                   +  dt    div (sqrt(2 D_k n_k / dt) Z)^n
       !                   +  dt    ext_src
       !
       ! in delta formulation:
       !
       ! (I - div (dt/2) D_k grad) delta n_k =   dt div (D_k grad n_k^n)
       !                                       + dt div (sqrt(2 D_k n_k / dt) Z)^n
       !                                       + dt ext_src
       !
       ! we combine the entire rhs into stoch_fluxdiv
       do n=1,nlevs
          call multifab_plus_plus(stoch_fluxdiv(n),ext_src(n),0)
          call multifab_plus_plus(stoch_fluxdiv(n),diff_fluxdiv(n),0)
          call multifab_mult_mult_s(stoch_fluxdiv(n),dt)
       end do
       call implicit_diffusion(mla,n_old,n_new,stoch_fluxdiv,diff_coef_face,dx,dt,the_bc_tower)
        */
    } else if (reactDiff_diffusion_type == 2) {
        Abort("AdvanceDiffusion() - write explicit midpoint scheme");

        /*
! explicit midpoint scheme

       ! n_k^{n+1/2} = n_k^n + (dt/2) div (D_k grad n_k)^n
       !                     + (dt/2) div (sqrt(2 D_k n_k / (dt/2) ) Z_1)^n
       !                     + (dt/2) ext_src
       do n=1,nlevs
          call multifab_copy_c(n_new(n),1,n_old(n),1,nspecies,0)
          call multifab_saxpy_3(n_new(n),dt/2.d0      ,diff_fluxdiv(n))
          call multifab_saxpy_3(n_new(n),dt/sqrt(2.d0),stoch_fluxdiv(n))
          call multifab_saxpy_3(n_new(n),dt/2.d0      ,ext_src(n))
          call multifab_fill_boundary(n_new(n))
          call multifab_physbc(n_new(n),1,scal_bc_comp,nspecies, &
                               the_bc_tower%bc_tower_array(n),dx_in=dx(n,:))
       end do

       ! compute diffusive flux divergence at t^{n+1/2}
       call diffusive_n_fluxdiv(mla,n_new,diff_coef_face,diff_fluxdiv,dx,the_bc_tower)

       if (variance_coef_mass .gt. 0.d0) then
          ! fill random flux multifabs with new random numbers
          call fill_mass_stochastic(mla,the_bc_tower%bc_tower_array)

          ! compute second-stage stochastic flux divergence and
          ! add to first-stage stochastic flux divergence
          select case (midpoint_stoch_flux_type)
          case (1)
             ! use n_old
             call stochastic_n_fluxdiv(mla,n_old,diff_coef_face,stoch_fluxdiv,dx,dt, &
                                       the_bc_tower,increment_in=.true.)
          case (2)
             ! use n_pred 
             call stochastic_n_fluxdiv(mla,n_new,diff_coef_face,stoch_fluxdiv,dx,dt, &
                                       the_bc_tower,increment_in=.true.)
          case (3)
             ! We use n_new=2*n_pred-n_old here as temporary storage since we will overwrite it shortly
             do n=1,nlevs
                call multifab_mult_mult_s_c(n_new(n),1,2.d0,nspecies,n_new(n)%ng)
                call multifab_sub_sub_c(n_new(n),1,n_old(n),1,nspecies,n_new(n)%ng)
             end do
             ! use n_new=2*n_pred-n_old
             call stochastic_n_fluxdiv(mla,n_new,diff_coef_face,stoch_fluxdiv,dx,dt, &
                                       the_bc_tower,increment_in=.true.)
          case default
             call bl_error("advance_diffusion: invalid midpoint_stoch_flux_type")
          end select
       end if
       
       ! n_k^{n+1} = n_k^n + dt div (D_k grad n_k)^{n+1/2}
       !                   + dt div (sqrt(2 D_k n_k^n dt) Z_1 / sqrt(2) )
       !                   + dt div (sqrt(2 D_k n_k^? dt) Z_2 / sqrt(2) )
       !                   + dt ext_src
       ! where
       ! n_k^? = n_k^n               (midpoint_stoch_flux_type=1)
       !       = n_k^pred            (midpoint_stoch_flux_type=2)
       !       = 2*n_k^pred - n_k^n  (midpoint_stoch_flux_type=3)
       do n=1,nlevs
          call multifab_copy_c(n_new(n),1,n_old(n),1,nspecies,0)
          call multifab_saxpy_3(n_new(n),dt           ,diff_fluxdiv(n))
          call multifab_saxpy_3(n_new(n),dt/sqrt(2.d0),stoch_fluxdiv(n))
          call multifab_saxpy_3(n_new(n),dt           ,ext_src(n))
          call multifab_fill_boundary(n_new(n))
          call multifab_physbc(n_new(n),1,scal_bc_comp,nspecies, &
                               the_bc_tower%bc_tower_array(n),dx_in=dx(n,:))
       end do
         */

    } else {
        Abort("AdvanceDiffusion() - invalid reactDiff_diffusion_type");
    }
    
}


void DiffusiveNFluxdiv(MultiFab& n_in,
                       MultiFab& diff_fluxdiv,
                       const Geometry& geom,
                       const Real& time) {

    // fill n ghost cells
    n_in.FillBoundary(geom.periodicity());
    MultiFabPhysBC(n_in, geom, 0, nspecies, SPEC_BC_COMP, time);

    BoxArray ba = n_in.boxArray();
    DistributionMapping dmap = n_in.DistributionMap();
    
    // don't need to set much here for explicit evaluations
    LPInfo info;

    // operator of the form (ascalar * acoef - bscalar div bcoef grad) phi
    MLABecLaplacian mlabec({geom}, {ba}, {dmap}, info);
    mlabec.setMaxOrder(2);

    // store one component at a time and take L(phi) one component at a time
    MultiFab phi (ba,dmap,1,1);
    MultiFab Lphi(ba,dmap,1,0);

    MultiFab acoef(ba,dmap,1,0);
    std::array< MultiFab, AMREX_SPACEDIM > bcoef;
    AMREX_D_TERM(bcoef[0].define(convert(ba,nodal_flag_x), dmap, 1, 0);,
                 bcoef[1].define(convert(ba,nodal_flag_y), dmap, 1, 0);,
                 bcoef[2].define(convert(ba,nodal_flag_z), dmap, 1, 0););
    
    // build array of boundary conditions needed by MLABecLaplacian
    std::array<LinOpBCType, AMREX_SPACEDIM> lo_mlmg_bc;
    std::array<LinOpBCType, AMREX_SPACEDIM> hi_mlmg_bc;

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {
        if (bc_mass_lo[idim] == -1 || bc_mass_hi[idim] == -1) {
            if ( !(bc_mass_lo[idim] == -1 && bc_mass_hi[idim] == -1) ) {
                Abort("Both bc_mass_lo and bc_mass_hi must be periodic in a given direction if the other one is");
            }            
            lo_mlmg_bc[idim] = LinOpBCType::Periodic;            
            hi_mlmg_bc[idim] = LinOpBCType::Periodic;
        }

        if (bc_mass_lo[idim] == 0) {
            lo_mlmg_bc[idim] = LinOpBCType::inhomogNeumann;
        } else if (bc_mass_lo[idim] == 1) {
            lo_mlmg_bc[idim] = LinOpBCType::Dirichlet;
        } else if (bc_mass_lo[idim] != -1) {
            Abort("Invalid bc_mass_lo");
        }

        if (bc_mass_hi[idim] == 0) {
            hi_mlmg_bc[idim] = LinOpBCType::inhomogNeumann;
        } else if (bc_mass_hi[idim] == 1) {
            hi_mlmg_bc[idim] = LinOpBCType::Dirichlet;
        } else if (bc_mass_hi[idim] != -1) {
            Abort("Invalid bc_mass_hi");
        }
    }

    mlabec.setDomainBC(lo_mlmg_bc,hi_mlmg_bc);

    // set acoeff to 0and bcoeff to -1
    mlabec.setScalars(0., -1.);

    acoef.setVal(0.);
    mlabec.setACoeffs(0, acoef);

    for (int i=0; i<nspecies; ++i) {

        // copy ith component of n_in into phi, including ghost cells
        MultiFab::Copy(phi,n_in,i,0,1,1);

        // load D_fick for species i into bcoef
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            bcoef[d].setVal(D_Fick[i]);
        }
        mlabec.setBCoeffs(0, amrex::GetArrOfConstPtrs(bcoef));

        MLMG mlmg(mlabec);

        mlmg.apply({&Lphi},{&phi});

        MultiFab::Copy(diff_fluxdiv,Lphi,0,i,1,0);
        
    }
    
}
