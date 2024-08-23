#include "reactDiff_functions.H"
#include "chemistry_functions.H"

void AdvanceDiffusion(const MultiFab& n_old,
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
