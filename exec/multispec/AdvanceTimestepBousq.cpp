
#include "hydro_functions.H"
#include "common_functions.H"
#include "gmres_functions.H"
#include "multispec_functions.H"

#include "StochMomFlux.H"


#include <AMReX_ParallelDescriptor.H>
#include <AMReX_MultiFabUtil.H>


// argv contains the name of the inputs file entered at the command line
void AdvanceTimestepBousq(std::array< MultiFab, AMREX_SPACEDIM >& umac,
                          MultiFab& rho_old,
                          MultiFab& rho_new,
                          MultiFab& rhotot_old,
                          MultiFab& rhotot_new,
                          MultiFab& pi,
                          MultiFab& eta, 
                          std::array< MultiFab, NUM_EDGE >&  eta_ed,
                          MultiFab& kappa, MultiFab& Temp,
                          std::array< MultiFab, NUM_EDGE >& Temp_ed,
                          MultiFab& diff_mass_fluxdiv,
                          MultiFab& stoch_mass_fluxdiv,
                          std::array< MultiFab, AMREX_SPACEDIM >& stoch_mass_flux,
                          std::array< MultiFab, AMREX_SPACEDIM >& grad_Epot_old,
                          std::array< MultiFab, AMREX_SPACEDIM >& grad_Epot_new,
                          MultiFab& charge_old,
                          MultiFab& charge_new,
                          MultiFab& Epot,
                          MultiFab& permittivity,
                          StochMassFlux& sMassFlux,
                          StochMomFlux& sMomFlux,
                          const Real& dt,
                          const Real& time,
                          const int& istep,
                          const Geometry& geom)
{
  
    BL_PROFILE_VAR("AdvanceTimestepBousq()",AdvanceTimestepBousq);

    BoxArray ba = rho_old.boxArray();
    DistributionMapping dmap = rho_old.DistributionMap();

    const Real* dx = geom.CellSize();
    
    Vector<Real> weights_mom(1);
    Vector<Real> weights_mass(2);

    if (barodiffusion_type != 0) {
        Abort("AdvanceTimestepBousq: barodiffusion not supported yet");
    }

    // For BDS we need to project edge states onto constraints
    int proj_type = (use_charged_fluid && electroneutral) ? 4 : 3;
    
    Real theta_alpha = 1./dt;

    Real relxn_param_charge_in;
    Real norm_pre_rhs;

    MultiFab adv_mass_fluxdiv(ba,dmap,nspecies,0);
    MultiFab gmres_rhs_p     (ba,dmap,       1,0);
    MultiFab dpi             (ba,dmap,       1,1);
    
    std::array< MultiFab, AMREX_SPACEDIM > umac_old;
    std::array< MultiFab, AMREX_SPACEDIM > mtemp;
    std::array< MultiFab, AMREX_SPACEDIM > adv_mom_fluxdiv_old;
    std::array< MultiFab, AMREX_SPACEDIM > adv_mom_fluxdiv_new;
    std::array< MultiFab, AMREX_SPACEDIM > diff_mom_fluxdiv_old;
    std::array< MultiFab, AMREX_SPACEDIM > diff_mom_fluxdiv_new;
    std::array< MultiFab, AMREX_SPACEDIM > stoch_mom_fluxdiv;
    std::array< MultiFab, AMREX_SPACEDIM > gmres_rhs_v;
    std::array< MultiFab, AMREX_SPACEDIM > dumac;
    std::array< MultiFab, AMREX_SPACEDIM > gradpi;
    std::array< MultiFab, AMREX_SPACEDIM > rho_fc;
    std::array< MultiFab, AMREX_SPACEDIM > rhotot_fc_old;
    std::array< MultiFab, AMREX_SPACEDIM > rhotot_fc_new;
    std::array< MultiFab, AMREX_SPACEDIM > diff_mass_flux;

    // only used when variance_coef_mass>0 and midpoint_stoch_mass_flux_type=2
    MultiFab stoch_mass_fluxdiv_old;
    std::array< MultiFab, AMREX_SPACEDIM > stoch_mass_flux_old;

    // only used when use_charged_fluid=1
    std::array< MultiFab, AMREX_SPACEDIM > Lorentz_force;

    // only used when use_multiphase=1
    std::array< MultiFab, AMREX_SPACEDIM > div_reversible_stress;

    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        umac_old[d]            .define(convert(ba,nodal_flag_dir[d]), dmap,        1, 1);
        mtemp[d]               .define(convert(ba,nodal_flag_dir[d]), dmap,        1, 1);
        adv_mom_fluxdiv_old[d] .define(convert(ba,nodal_flag_dir[d]), dmap,        1, 0);
        adv_mom_fluxdiv_new[d] .define(convert(ba,nodal_flag_dir[d]), dmap,        1, 0);
        diff_mom_fluxdiv_old[d].define(convert(ba,nodal_flag_dir[d]), dmap,        1, 0);
        diff_mom_fluxdiv_new[d].define(convert(ba,nodal_flag_dir[d]), dmap,        1, 0);
        stoch_mom_fluxdiv[d]   .define(convert(ba,nodal_flag_dir[d]), dmap,        1, 0);
        gmres_rhs_v[d]         .define(convert(ba,nodal_flag_dir[d]), dmap,        1, 0);
        dumac[d]               .define(convert(ba,nodal_flag_dir[d]), dmap,        1, 1);
        gradpi[d]              .define(convert(ba,nodal_flag_dir[d]), dmap,        1, 0);
        rho_fc[d]              .define(convert(ba,nodal_flag_dir[d]), dmap, nspecies, 0);
        rhotot_fc_old[d]       .define(convert(ba,nodal_flag_dir[d]), dmap,        1, 1);
        rhotot_fc_new[d]       .define(convert(ba,nodal_flag_dir[d]), dmap,        1, 1);
    }

    // for ito interpretation we need to save stoch_mass_fluxdiv_old 
    if (variance_coef_mass != 0. && midpoint_stoch_mass_flux_type == 2) {
        stoch_mass_fluxdiv_old.define(ba,dmap,nspecies,0);
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            stoch_mass_flux_old[d].define(convert(ba,nodal_flag_dir[d]), dmap, nspecies, 0);
        }
    }

    if (use_charged_fluid) {
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            Lorentz_force[d].define(convert(ba,nodal_flag_dir[d]), dmap, 1, 0);
        }
    }
    
    if (use_multiphase) {
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            div_reversible_stress[d].define(convert(ba,nodal_flag_dir[d]), dmap, 1, 0);
        }
    }

    // make a copy of umac at t^n
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        MultiFab::Copy(umac_old[d],umac[d],0,0,1,1);
    }

    //////////////////////////////////////////////
    // Step 1: solve for v^{n+1,*} and pi^{n+1/2,*} using GMRES
    //////////////////////////////////////////////

    //  average rho_i^n and rho^n to faces
    AverageCCToFace(   rho_old,   rho_fc    ,0,nspecies,SPEC_BC_COMP,geom);
    AverageCCToFace(rhotot_old,rhotot_fc_old,0,1       , RHO_BC_COMP,geom);

    // compute mtemp = (rho*v)^n
    ConvertMToUmac(rhotot_fc_old,umac,mtemp,0);

    // set gmres_rhs_v = mtemp / dt
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        MultiFab::Copy(gmres_rhs_v[d],mtemp[d],0,0,1,0);
        gmres_rhs_v[d].mult(1./dt,0);
    }

    // compute adv_mom_fluxdiv_old = -rho*v^n*v^n
    // save this for use in the corrector GMRES solve
    MkAdvMFluxdiv(umac,mtemp,adv_mom_fluxdiv_old,dx,0);

    // add -rho*v^n,v^n to gmres_rhs_v
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        MultiFab::Add(gmres_rhs_v[d],adv_mom_fluxdiv_old[d],0,0,1,0);
    }
    
    // compute diff_mom_fluxdiv_old = L_0^n v^n
    // save this for use in the corrector GMRES solve
    MkDiffusiveMFluxdiv(diff_mom_fluxdiv_old,umac,eta,eta_ed,kappa,geom,dx,0);


    // add (1/2)*diff_mom_fluxdiv_old to gmres_rhs_v
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        MultiFab::Saxpy(gmres_rhs_v[d],0.5,diff_mom_fluxdiv_old[d],0,0,1,0);
    }

    if (variance_coef_mom != 0.) {

        // fill the stochastic momentum multifabs with new sets of random numbers
        sMomFlux.fillMomStochastic();

        // compute stoch_mom_fluxdiv = div (sqrt(eta^n...) Wbar^n)
        // save this for use in the corrector GMRES solve
        weights_mom[0] = 1.;
        sMomFlux.StochMomFluxDiv(stoch_mom_fluxdiv,0,eta,eta_ed,Temp,Temp_ed,weights_mom,dt);

        // add stochastic momentum fluxes to gmres_rhs_v
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            MultiFab::Saxpy(gmres_rhs_v[d],1.,stoch_mom_fluxdiv[d],0,0,1,0);
        }

    }

    // gravity
    bool any_grav = false;
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        if (grav[d] != 0.) any_grav = true;
    }
    if (any_grav) {
        Abort("AdvanceTimestepBousq.cpp gravity not implemented");
    }

    if (variance_coef_mass != 0.) {
        sMassFlux.fillMassStochastic();
        if (chk_int > 0) {
            Abort("AdvanceTimestepBousq.cpp checkpointing not implemented");
        }
    }

    // compute diffusive, stochastic, potential mass fluxes
    // with barodiffusion and thermodiffusion
    // this computes "-F = rho W chi [Gamma grad x... ]"
    weights_mass[0] = 1.;
    weights_mass[1] = 0.;

    // For electroneutral, we only have to do charge relaxation in the corrector
    relxn_param_charge_in=relxn_param_charge;
    relxn_param_charge=0.0; // Don't correct in predictor

    ComputeMassFluxdiv(rho_old,rhotot_old,Temp,diff_mass_fluxdiv,stoch_mass_fluxdiv,
                       diff_mass_flux,stoch_mass_flux,sMassFlux,0.5*dt,time,geom,weights_mass,
                       charge_old,grad_Epot_old,Epot,permittivity);
    
    // here is a reasonable place to call something to compute in reversible stress term
    // in this case want to get divergence so it looks like a add to rhs for stokes solver
    if (use_multiphase) {

        // compute reversible stress tensor ---added term
        // call compute_div_reversible_stress(mla,div_reversible_stress,rhotot_old,rho_old,dx,the_bc_tower)
        Abort("AdvanceTimestepBousq.cpp write ComputeDivReversibleStress()");

        // add divergence of reversible stress to gmres_rhs_v
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            MultiFab::Saxpy(gmres_rhs_v[d],1.,div_reversible_stress[d],0,0,1,0);
        }

    }

    if (use_charged_fluid) {

        // compute old Lorentz force
        ComputeLorentzForce(Lorentz_force,grad_Epot_old,permittivity,charge_old,geom);

        // add Lorentz force to gmres_rhs_v
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            MultiFab::Saxpy(gmres_rhs_v[d],1.0,Lorentz_force[d],0,0,1,0);
        }

    }

    // compute grad pi^{n-1/2}
    ComputeGrad(pi,gradpi,0,0,1,PRES_BC_COMP,geom);
    
    // subtract grad pi^{n-1/2} from gmres_rhs_v
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        MultiFab::Subtract(gmres_rhs_v[d],gradpi[d],0,0,1,0);
    }

    // set inhomogeneous velocity bc's to values supplied in inhomogeneous_bc_val
    //
    //
    //

    // modify umac to respect the boundary conditions we want after the next gmres solve
    // thus when we add L_0^n vbar^n to gmres_rhs_v and add div vbar^n to gmres_rhs_p we
    // are automatically putting the system in delta form WITH homogeneous boundary conditions
    for (int i=0; i<AMREX_SPACEDIM; ++i) {
        // set normal velocity of physical domain boundaries
        MultiFabPhysBCDomainVel(umac[i],geom,i);
        // set transverse velocity behind physical boundaries
        MultiFabPhysBCMacVel(umac[i],geom,i);
        // fill periodic and interior ghost cells
        umac[i].FillBoundary(geom.periodicity());
    }

    // compute mtemp = (rho*vbar)^n
    ConvertMToUmac(rhotot_fc_old,umac,mtemp,0);

    // add -(mtemp/dt) to gmres_rhs_v
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        MultiFab::Saxpy(gmres_rhs_v[d],-1./dt,mtemp[d],0,0,1,0);
    }

    // compute diff_mom_fluxdiv_new = L_0^n vbar^n
    MkDiffusiveMFluxdiv(diff_mom_fluxdiv_new,umac,eta,eta_ed,kappa,geom,dx,0);

    // add (1/2)*diff_mom_fluxdiv_new to gmres_rhs_v
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        MultiFab::Saxpy(gmres_rhs_v[d],0.5,diff_mom_fluxdiv_new[d],0,0,1,0);
    }

    // compute div(vbar^n) and store in gmres_rhs_p
    // the sign convention is correct since we solve -div(delta v) = div(vbar^n)
    ComputeDiv(gmres_rhs_p,umac,0,0,1,geom,0);

    // multiply eta and kappa by 1/2 to put in proper form for gmres solve
    eta.mult  (0.5,0,1,1);
    kappa.mult(0.5,0,1,1);
    for (int d=0; d<NUM_EDGE; ++d) {
        eta_ed[d].mult(0.5,0,1,0);
    }

    // set the initial guess to zero
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        dumac[d].setVal(0.);
    }
    dpi.setVal(0.);

    // zero gmres_rhs_v on physical boundaries
    ZeroEdgevalPhysical(gmres_rhs_v, geom, 0, 1);

    Real gmres_abs_tol_in = gmres_abs_tol; // save this

    // This relies entirely on relative tolerance and can fail if the rhs is roundoff error only:
    // gmres_abs_tol = 0.d0 ! It is better to set gmres_abs_tol in namelist to a sensible value

    // call gmres to compute delta v and delta pi
    GMRES gmres(ba,dmap,geom);
    gmres.Solve(gmres_rhs_v, gmres_rhs_p, dumac, dpi, rhotot_fc_old, eta, eta_ed,
                kappa, theta_alpha, geom, norm_pre_rhs);

    // for the corrector gmres solve we want the stopping criteria based on the
    // norm of the preconditioned rhs from the predictor gmres solve.  otherwise
    // for cases where du in the corrector should be small the gmres stalls
    gmres_abs_tol = amrex::max(gmres_abs_tol_in, norm_pre_rhs*gmres_rel_tol);

    // compute v^{n+1,*} = vbar^n + dumac
    // compute pi^{n+1/2,*} = pi^{n-1/2} + dpi
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        MultiFab::Add(umac[d],dumac[d],0,0,1,0);
    }
    MultiFab::Add(pi,dpi,0,0,1,0);

    // pressure ghost cells
    pi.FillBoundary(geom.periodicity());
    MultiFabPhysBC(pi,geom,0,1,PRES_BC_COMP);
    
    for (int i=0; i<AMREX_SPACEDIM; ++i) {
        // set normal velocity of physical domain boundaries
        MultiFabPhysBCDomainVel(umac[i],geom,i);
        // set transverse velocity behind physical boundaries
        MultiFabPhysBCMacVel(umac[i],geom,i);
        // fill periodic and interior ghost cells
        umac[i].FillBoundary(geom.periodicity());
    }

    // restore eta and kappa
    eta.mult  (2.,0,1,1);
    kappa.mult(2.,0,1,1);
    for (int d=0; d<NUM_EDGE; ++d) {
        eta_ed[d].mult(2.,0,1,0);
    }

    //////////////////////////////////////////////
    // Step 2: compute reactions and mass fluxes at t^n
    //////////////////////////////////////////////

    // compute chemical rates m_i*R^n_i
    /*
    if (nreactions > 0) then

    end if
    */

    if (advection_type >= 1) {

        Abort("AdvanceTimestepBousq.cpp advection_type >= 1 not supported");

    } else {

        // compute adv_mass_fluxdiv = -rho_i^n * v^n and then
        // increment adv_mass_fluxdiv by -rho_i^n * v^{n+1,*}
        MkAdvSFluxdiv(umac_old,rho_fc,adv_mass_fluxdiv,geom,0,nspecies,false);
        MkAdvSFluxdiv(umac    ,rho_fc,adv_mass_fluxdiv,geom,0,nspecies,true);
    }

    //////////////////////////////////////////////
    /// Step 3: density prediction to t^{n+1/2}
    //////////////////////////////////////////////

    if (advection_type >= 1) {
        /*
        do n=1,nlevs
          call multifab_saxpy_4(rho_new(n),rho_old(n),0.5d0*dt,rho_update(n))
        end do
        */
    } else {

        // compute rho_i^{n+1/2} (store in rho_new)
        // multiply adv_mass_fluxdiv by (1/4) since it contains -rho_i^n * (v^n + v^{n+1,*})
        MultiFab::Copy(rho_new,rho_old,0,0,nspecies,0);
        MultiFab::Saxpy(rho_new,0.25*dt, adv_mass_fluxdiv,0,0,nspecies,0);
        MultiFab::Saxpy(rho_new,0.50*dt,diff_mass_fluxdiv,0,0,nspecies,0);
        if (variance_coef_mass != 0.) {
            MultiFab::Saxpy(rho_new,0.50*dt,stoch_mass_fluxdiv,0,0,nspecies,0);
        }
        /*
        if (nreactions > 0) then
           call multifab_saxpy_3_cc(rho_new(n),1,0.5d0*dt,chem_rate(n),1,nspecies)
        end if
        */
    }

    // compute rhotot^{n+1/2} from rho^{n+1/2} in VALID REGION
    ComputeRhotot(rho_new,rhotot_new);

    // fill rho and rhotot ghost cells at t^{n+1/2}
    FillRhoRhototGhost(rho_new,rhotot_new,geom);
    
    // average rho_i^{n+1/2} to faces
    AverageCCToFace(rho_new,rho_fc,0,nspecies,SPEC_BC_COMP,geom);

    if (use_charged_fluid) {
        // compute total charge at t^{n+1/2}
        DotWithZ(rho_new,charge_new);

        // compute permittivity at t^{n+1/2}
        if (dielectric_type != 0) {
            Abort("AdvanceTimestepBousq dielectric_type != 0");
        }
    }

    // compute (eta,kappa)^{n+1/2}
    if (mixture_type != 0) {
        Abort("AdvanceTimestepBousq mixture_type != 0");
    }

    //////////////////////////////////////////////
    // Step 4: compute mass fluxes and reactions at t^{n+1/2}
    //////////////////////////////////////////////

    // compute mass fluxes and reactions at t^{n+1/2}
    // For electroneutral, enable charge correction in the corrector
    relxn_param_charge=relxn_param_charge_in; // Default value is 1

#if 0
    if (midpoint_stoch_mass_flux_type .eq. 1) then
       ! strato

       ! compute diffusive, stochastic, potential mass fluxes
       ! with barodiffusion and thermodiffusion
       ! this computes "-F = rho W chi [Gamma grad x... ]"
       weights(:) = 1.d0/sqrt(2.d0)
       call compute_mass_fluxdiv(mla,rho_new,rhotot_new,gradp_baro,Temp, &
                                 diff_mass_fluxdiv,stoch_mass_fluxdiv, &
                                 diff_mass_flux,stoch_mass_flux, &
                                 dt,time,dx,weights,the_bc_tower, &
                                 charge_new,grad_Epot_new,Epot, &
                                 permittivity,zero_initial_Epot_in=.false.)

       ! begin to assemble total fluxes (diffusive plus stochastic);
       ! we still need to add add advective flux later
       if(present(total_mass_flux)) then
          do n=1,nlevs
             do i=1,dm
                call multifab_copy_c(total_mass_flux(n,i),1,diff_mass_flux(n,i),1,nspecies,0)
                if (variance_coef_mass .ne. 0.d0) then
                   call multifab_plus_plus_c(total_mass_flux(n,i),1,stoch_mass_flux(n,i),1,nspecies,0)
                end if
             end do
          end do
       end if

    else if (midpoint_stoch_mass_flux_type .eq. 2) then
       ! ito

       if(present(total_mass_flux)) then
          call bl_error('Code cannot presently return total_mass_flux if midpoint_stoch_mass_flux_type!=1')
       end if

       if (variance_coef_mass .ne. 0.d0) then
          ! for ito interpretation we need to save stoch_mass_fluxdiv_old here
          ! then later add it to stoch_mass_fluxdiv and multiply by 1/2
          do n=1,nlevs
             call multifab_copy_c(stoch_mass_fluxdiv_old(n),1,stoch_mass_fluxdiv(n),1,nspecies,0)
             do i=1,dm
                call multifab_copy_c(stoch_mass_flux_old(n,i),1,stoch_mass_flux(n,i),1,nspecies,0)
             end do
          end do
       end if

       ! compute diffusive, stochastic, potential mass fluxes
       ! with barodiffusion and thermodiffusion
       ! this computes "-F = rho W chi [Gamma grad x... ]"
       weights(1) = 0.d0
       weights(2) = 1.d0
       call compute_mass_fluxdiv(mla,rho_new,rhotot_new,gradp_baro,Temp, &
                                 diff_mass_fluxdiv,stoch_mass_fluxdiv, &
                                 diff_mass_flux,stoch_mass_flux, &
                                 0.5d0*dt,time,dx,weights,the_bc_tower, &
                                 charge_new,grad_Epot_new,Epot, &
                                 permittivity,zero_initial_Epot_in=.false.)

       if (variance_coef_mass .ne. 0.d0) then
          ! add stoch_mass_fluxdiv_old to stoch_mass_fluxdiv and multiply by 1/2
          do n=1,nlevs
             call multifab_plus_plus_c(stoch_mass_fluxdiv(n),1,stoch_mass_fluxdiv_old(n),1,nspecies,0)
             call multifab_mult_mult_s_c(stoch_mass_fluxdiv(n),1,0.5d0,nspecies,0)
          end do
       end if

    end if
    
    if(use_multiphase) then

       !compute reversible stress tensor ---added term (will add to gmres_rhs_v later)
       call compute_div_reversible_stress(mla,div_reversible_stress,rhotot_new,rho_new,dx,the_bc_tower)

    end if

    if (use_charged_fluid) then

       ! compute Lorentz force (using midpoint value not trapezoidal, will add to gmres_rhs_v later)
       call compute_Lorentz_force(mla,Lorentz_force,grad_Epot_new,permittivity, &
                                  charge_new,dx,the_bc_tower)

    end if

    ! compute chemical rates m_i*R^{n+1/2}_i
    if (nreactions > 0) then
       ! convert rho_old (at n) and rho_new (at n+1/2) (mass densities rho_i)
       ! into n_old and n_new (number densities n_i=rho_i/m_i)
       do n=1,nlevs
          call multifab_copy_c(n_old(n),1,rho_old(n),1,nspecies,0)
          call multifab_copy_c(n_new(n),1,rho_new(n),1,nspecies,0)
          do i=1,nspecies
             call multifab_div_div_s_c(n_old(n),i,molmass(i),1,0)
             call multifab_div_div_s_c(n_new(n),i,molmass(i),1,0)
          end do
       end do

       ! compute chemical rates R_i (units=[number density]/[time]) for the second half step
       call chemical_rates(mla,n_old,chem_rate_temp,dx,0.5d0*dt,n_new,mattingly_lin_comb_coef)

       ! convert chemical rates R_i into m_i*R_i (units=[mass density]/[time])
       do n=1,nlevs
          do i=1,nspecies
             call multifab_mult_mult_s_c(chem_rate_temp(n),i,molmass(i),1,0)
          end do
       end do

       ! compute chemical rates m_i*R^{n+1/2}_i for the full step
       do n=1,nlevs
          call multifab_plus_plus_c(chem_rate(n),1,chem_rate_temp(n),1,nspecies,0)
          call multifab_mult_mult_s_c(chem_rate(n),1,0.5d0,nspecies,0)
       end do

       if (use_bl_rng) then
          ! save random state for restart
          call bl_rng_copy_engine(rng_eng_reaction_chk,rng_eng_reaction)
       end if

    end if

    if (advection_type .ge. 1) then

       if(present(total_mass_flux)) then
          call bl_error('Code cannot presently return total_mass_flux if using BDS advection')
       end if

       ! add the diff/stoch/react terms to rho_update
       do n=1,nlevs
          call multifab_copy_c(rho_update(n),1,diff_mass_fluxdiv(n),1,nspecies)
          if (variance_coef_mass .ne. 0.d0) then
             call multifab_plus_plus_c(rho_update(n),1,stoch_mass_fluxdiv(n),1,nspecies)
          end if
          if (nreactions > 0) then
             call multifab_plus_plus_c(rho_update(n),1,chem_rate(n),1,nspecies,0)
          end if
       end do

       do n=1,nlevs
          ! set to zero to make sure ghost cells behind physical boundaries don't have NaNs
          call setval(bds_force(n),0.d0,all=.true.)
          call multifab_copy_c(bds_force(n),1,rho_update(n),1,nspecies,0)
          call multifab_fill_boundary(bds_force(n))
       end do

       if (advection_type .eq. 1 .or. advection_type .eq. 2) then
          ! bds increments rho_update with the advection term
          call bds(mla,umac_tmp,rho_tmp,rho_update,bds_force,rho_fc,rho_nd,dx,dt,1, &
                   nspecies,c_bc_comp,the_bc_tower,proj_type_in=proj_type)

       else if (advection_type .eq. 3 .or. advection_type .eq. 4) then
          call bds_quad(mla,umac_tmp,rho_old,rho_update,bds_force,rho_fc,dx,dt,1,nspecies, &
                        c_bc_comp,the_bc_tower,proj_type_in=proj_type)
       end if

       do n=1,nlevs
          call multifab_destroy(rho_tmp(n))
          call multifab_destroy(bds_force(n))
          call multifab_destroy(rho_nd(n))
          do i=1,dm
             call multifab_destroy(umac_tmp(n,i))
          end do
       end do

    else

       ! compute adv_mass_fluxdiv = -rho_i^{n+1/2} * v^n and
       ! increment adv_mass_fluxdiv by -rho_i^{n+1/2} * v^{n+1,*}
       call mk_advective_s_fluxdiv(mla,umac_old,rho_fc,adv_mass_fluxdiv,.false.,dx,1,nspecies)
       call mk_advective_s_fluxdiv(mla,umac    ,rho_fc,adv_mass_fluxdiv,.true. ,dx,1,nspecies)

       ! Add advective fluxes to total fluxes and fix sign for fluxes to be actual fluxes not negative of flux ;-)
       if(present(total_mass_flux)) then

          do n=1,nlevs
             do i=1,dm

                ! compute advective mass flux contribution, rho_i^{n+1/2} * v^n
                do comp=1,nspecies
                   call multifab_copy_c(adv_mass_flux(n,i),comp,umac_old(n,i),1,1,0)
                end do
                call multifab_mult_mult_c(adv_mass_flux(n,i),1,rho_fc(n,i),1,nspecies,0)

                ! now *subtract* off half of these advective fluxes from total_mass_flux
                call multifab_saxpy_3_cc(total_mass_flux(n,i),1,-0.5d0,adv_mass_flux(n,i),1,nspecies)

                ! compute advective mass flux contribution, rho_i^{n+1/2} * v^{n+1,*}
                do comp=1,nspecies
                   call multifab_copy_c(adv_mass_flux(n,i),comp,umac(n,i),1,1,0)
                end do
                call multifab_mult_mult_c(adv_mass_flux(n,i),1,rho_fc(n,i),1,nspecies,0)

                ! now *subtract* off half of these advective fluxes from total_mass_flux
                call multifab_saxpy_3_cc(total_mass_flux(n,i),1,-0.5d0,adv_mass_flux(n,i),1,nspecies)
                
                ! multiply the total mass flux w adv by -1, so that we take the 
                ! standard convention that the mass flux has the same sign as the fluid flow 
                call multifab_mult_mult_s_c(total_mass_flux(n,i),1,-1.d0,nspecies,0)

             enddo
          end do

       do n=1,nlevs
          do i=1,dm
                
          enddo
       enddo        

       end if

    end if

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Step 5: density integration to t^{n+1}
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (advection_type .ge. 1) then

       do n=1,nlevs
          call multifab_saxpy_4(rho_new(n),rho_old(n),dt,rho_update(n))
       end do

       do n=1,nlevs
          call multifab_destroy(rho_update(n))
       end do

    else

       ! compute rho_i^{n+1}
       ! multiply adv_mass_fluxdiv by (1/2) since it contains -rho_i^{n+1/2} * (v^n + v^{n+1,*})
       do n=1,nlevs
          call multifab_copy_c(rho_new(n),1,rho_old(n),1,nspecies,0)
          call multifab_saxpy_3_cc(rho_new(n),1,0.5d0*dt, adv_mass_fluxdiv(n),1,nspecies)
          call multifab_saxpy_3_cc(rho_new(n),1,      dt,diff_mass_fluxdiv(n),1,nspecies)
          if (variance_coef_mass .ne. 0.d0) then
             call multifab_saxpy_3_cc(rho_new(n),1,dt,stoch_mass_fluxdiv(n),1,nspecies)
          end if
          if (nreactions > 0) then
             call multifab_saxpy_3_cc(rho_new(n),1,dt,chem_rate(n),1,nspecies)
          end if
       end do

    end if

    if ( project_eos_int .gt. 0 .and. mod(istep,project_eos_int) .eq. 0) then
       call project_onto_eos(mla,rho_new)
    end if

    ! compute rhotot^{n+1} from rho^{n+1} in VALID REGION
    call compute_rhotot(mla,rho_new,rhotot_new)

    ! fill rho and rhotot ghost cells at t^{n+1}
    call fill_rho_rhotot_ghost(mla,rho_new,rhotot_new,dx,the_bc_tower)

    ! average rho^{n+1} to faces
    call average_cc_to_face(nlevs,rhotot_new,rhotot_fc_new,1,scal_bc_comp,1,the_bc_tower%bc_tower_array)

    if (use_charged_fluid) then
       ! compute total charge at t^{n+1}
       call dot_with_z(mla,rho_new,charge_new)
       ! compute permittivity at t^{n+1}
       call compute_permittivity(mla,permittivity,rho_new,rhotot_new,the_bc_tower)
    end if

    ! compute (eta,kappa)^{n+1}
    call compute_eta_kappa(mla,eta,eta_ed,kappa,rho_new,rhotot_new,Temp,dx, &
                           the_bc_tower%bc_tower_array)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Step 6: solve for v^{n+1} and pi^{n+1/2} using GMRES
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! compute mtemp = (rho*v)^n
    call convert_m_to_umac(mla,rhotot_fc_old,mtemp,umac_old,.false.)

    ! set gmres_rhs_v = mtemp / dt
    do n=1,nlevs
       do i=1,dm
          call multifab_copy_c(gmres_rhs_v(n,i),1,mtemp(n,i),1,1,0)
          call multifab_mult_mult_s_c(gmres_rhs_v(n,i),1,1.d0/dt,1,0)
       end do
    end do

    ! compute mtemp = rho^{n+1}*v^{n+1,*} for adv_mom_fluxdiv computation
    call convert_m_to_umac(mla,rhotot_fc_new,mtemp,umac,.false.)

    ! compute adv_mom_fluxdiv_new = -rho^{n+1}*v^{n+1,*}*v^{n+1,*}
    call mk_advective_m_fluxdiv(mla,umac,mtemp,adv_mom_fluxdiv_new,.false., &
                                dx,the_bc_tower%bc_tower_array)

    ! add (1/2) (-rho^n*v^n*v^n -rho^{n+1}*v^{n+1,*}*v^{n+1,*}) to gmres_rhs_v
    do n=1,nlevs
       do i=1,dm
          call multifab_saxpy_3_cc(gmres_rhs_v(n,i),1,0.5d0,adv_mom_fluxdiv_old(n,i),1,1)
          call multifab_saxpy_3_cc(gmres_rhs_v(n,i),1,0.5d0,adv_mom_fluxdiv_new(n,i),1,1)
       end do
    end do

    ! add (1/2)*diff_mom_fluxdiv_old to gmres_rhs_v
    do n=1,nlevs
       do i=1,dm
          call multifab_saxpy_3_cc(gmres_rhs_v(n,i),1,0.5d0,diff_mom_fluxdiv_old(n,i),1,1)
       end do
    end do

    if (variance_coef_mom .ne. 0.d0) then

       ! increment stoch_mom_fluxdiv by = div (sqrt(eta^{n+1}...) Wbar^n)
       weights(1) = 1.d0
       call stochastic_m_fluxdiv(mla,the_bc_tower%bc_tower_array,stoch_mom_fluxdiv,.true., &
                                 eta,eta_ed,Temp,Temp_ed,dx,dt,weights)

       ! add stochastic momentum fluxes to gmres_rhs_v
       ! note the factor of (1/2) since stoch_mom_fluxdiv contains
       ! div (sqrt(eta^{n+1}...) Wbar^n) + div (sqrt(eta^n...) Wbar^n)
       do n=1,nlevs
          do i=1,dm
             call multifab_saxpy_3_cc(gmres_rhs_v(n,i),1,0.5d0,stoch_mom_fluxdiv(n,i),1,1)
          end do
       end do

    end if

    ! gravity
    if (any(grav(1:dm) .ne. 0.d0)) then
       call mk_grav_force_bousq(mla,gmres_rhs_v,.true.,rho_fc,the_bc_tower)
    end if

    if(use_multiphase) then

      ! add divergence of reversible stress to gmres_rhs_v
      do n=1,nlevs
         do i=1,dm
            call multifab_saxpy_3(gmres_rhs_v(n,i),1.d0,div_reversible_stress(n,i))
         end do
      end do

    end if
   

    if (use_charged_fluid) then

       ! add Lorentz force to gmres_rhs_v
       do n=1,nlevs
          do i=1,dm
             call multifab_saxpy_3(gmres_rhs_v(n,i),1.d0,Lorentz_force(n,i))
          end do
       end do

    end if

    ! compute grad pi^{n+1/2,*}
    call compute_grad(mla,pi,gradpi,dx,1,pres_bc_comp,1,1,the_bc_tower%bc_tower_array)

    ! subtract grad pi^{n+1/2,*} from gmres_rhs_v
    do n=1,nlevs
       do i=1,dm
          call multifab_sub_sub_c(gmres_rhs_v(n,i),1,gradpi(n,i),1,1,0)
       end do
    end do

    ! modify umac to respect the boundary conditions we want after the next gmres solve
    ! thus when we add L_0^{n+1} vbar^{n+1,*} to gmres_rhs_v and add div vbar^{n+1,*} to gmres_rhs_p we
    ! are automatically putting the system in delta form WITH homogeneous boundary conditions
    do n=1,nlevs
       do i=1,dm
          ! set normal velocity on physical domain boundaries
          call multifab_physbc_domainvel(umac(n,i),vel_bc_comp+i-1, &
                                         the_bc_tower%bc_tower_array(n), &
                                         dx(n,:),vel_bc_n(n,:))
          ! set transverse velocity behind physical boundaries
          call multifab_physbc_macvel(umac(n,i),vel_bc_comp+i-1, &
                                      the_bc_tower%bc_tower_array(n), &
                                      dx(n,:),vel_bc_t(n,:))
          ! fill periodic and interior ghost cells
          call multifab_fill_boundary(umac(n,i))
       end do
    end do

    ! compute mtemp = rho^{n+1}*vbar^{n+1,*}
    call convert_m_to_umac(mla,rhotot_fc_new,mtemp,umac,.false.)

    ! add -(mtemp/dt) to gmres_rhs_v
    do n=1,nlevs
       do i=1,dm
          call multifab_saxpy_3_cc(gmres_rhs_v(n,i),1,-1.d0/dt,mtemp(n,i),1,1)
       end do
    end do

    ! compute diff_mom_fluxdiv_new = L_0^{n+1} vbar^{n+1,*}
    call diffusive_m_fluxdiv(mla,diff_mom_fluxdiv_new,.false.,umac,eta,eta_ed,kappa,dx, &
                             the_bc_tower%bc_tower_array)

    ! add (1/2)*diff_mom_fluxdiv_new to gmres_rhs_v
    do n=1,nlevs
       do i=1,dm
          call multifab_saxpy_3_cc(gmres_rhs_v(n,i),1,0.5d0,diff_mom_fluxdiv_new(n,i),1,1)
       end do
    end do

    ! compute div(vbar^{n+1,*}) and store in gmres_rhs_p
    ! the sign convention is correct since we solve -div(delta v) = div(vbar^{n+1,*})
    call compute_div(mla,umac,gmres_rhs_p,dx,1,1,1)

    ! multiply eta and kappa by 1/2 to put in proper form for gmres solve
    do n=1,nlevs
       call multifab_mult_mult_s_c(eta(n)  ,1,0.5d0,1,eta(n)%ng)
       call multifab_mult_mult_s_c(kappa(n),1,0.5d0,1,kappa(n)%ng)
       do i=1,size(eta_ed,dim=2)
          call multifab_mult_mult_s_c(eta_ed(n,i),1,0.5d0,1,eta_ed(n,i)%ng)
       end do
    end do

    ! set the initial guess to zero
    do n=1,nlevs
       do i=1,dm
          call multifab_setval(dumac(n,i),0.d0,all=.true.)
       end do
       call multifab_setval(dpi(n),0.d0,all=.true.)
    end do

    ! zero gmres_rhs_v on physical boundaries
    do n=1,nlevs
       call zero_edgeval_physical(gmres_rhs_v(n,:),1,1,the_bc_tower%bc_tower_array(n))
    end do

    ! call gmres to compute delta v and delta pi
    call gmres(mla,the_bc_tower,dx,gmres_rhs_v,gmres_rhs_p,dumac,dpi,rhotot_fc_new, &
               eta,eta_ed,kappa,theta_alpha,norm_pre_rhs)

    gmres_abs_tol = gmres_abs_tol_in ! Restore the desired tolerance   
       
    ! compute v^{n+1} = vbar^{n+1,*} + dumac
    ! compute pi^{n+1/2}= pi^{n+1/2,*} + dpi
    do n=1,nlevs
       do i=1,dm
          call multifab_plus_plus_c(umac(n,i),1,dumac(n,i),1,1,0)
       end do
       call multifab_plus_plus_c(pi(n),1,dpi(n),1,1,0)
    end do
       
    do n=1,nlevs
       ! presure ghost cells
       call multifab_fill_boundary(pi(n))
       call multifab_physbc(pi(n),1,pres_bc_comp,1,the_bc_tower%bc_tower_array(n), &
                            dx_in=dx(n,:))
       do i=1,dm
          ! set normal velocity on physical domain boundaries
          call multifab_physbc_domainvel(umac(n,i),vel_bc_comp+i-1, &
                                         the_bc_tower%bc_tower_array(n), &
                                         dx(n,:),vel_bc_n(n,:))
          ! set transverse velocity behind physical boundaries
          call multifab_physbc_macvel(umac(n,i),vel_bc_comp+i-1, &
                                      the_bc_tower%bc_tower_array(n), &
                                      dx(n,:),vel_bc_t(n,:))
          ! fill periodic and interior ghost cells
          call multifab_fill_boundary(umac(n,i))
       end do
    end do

    ! restore eta and kappa
    do n=1,nlevs
       call multifab_mult_mult_s_c(eta(n)  ,1,2.d0,1,eta(n)%ng)
       call multifab_mult_mult_s_c(kappa(n),1,2.d0,1,kappa(n)%ng)
       do i=1,size(eta_ed,dim=2)
          call multifab_mult_mult_s_c(eta_ed(n,i),1,2.d0,1,eta_ed(n,i)%ng)
       end do
    end do

    ! if writing a plotfile for electro-explicit option, compute Epot using new densities
    ! use existing machinery to compute all mass fluxes - yes this is overkill
    ! but better than writing a separate routine for now
    ! diff/stoch_mass_flux and fluxdiv are not persistent so changing values doesn't
    ! matter.  grad_Epot is persistent so we overwrite the old values,
    ! which are thrown away
    if (use_charged_fluid .and. (.not. electroneutral) ) then
    if (plot_int.gt.0 .and. ( (mod(istep,plot_int).eq.0) .or. (istep.eq.max_step)) ) then
       ! compute the new charge and store it in charge_old (it could get modified to
       ! subtract off the average in compute_mass_fluxdiv)
       call dot_with_z(mla,rho_new,charge_old)
       call compute_mass_fluxdiv(mla,rho_new,rhotot_new,gradp_baro,Temp, &
                                 diff_mass_fluxdiv,stoch_mass_fluxdiv, &
                                 diff_mass_flux,stoch_mass_flux, &
                                 dt,time,dx,weights,the_bc_tower, &
                                 charge_old,grad_Epot_old,Epot, &
                                 permittivity)
    end if
    end if
#endif

}
