#include "hydro_functions.H"
#include "common_functions.H"
#include "gmres_functions.H"
#include "multispec_functions.H"

#include "StochMomFlux.H"


#include <AMReX_ParallelDescriptor.H>
#include <AMReX_MultiFabUtil.H>


// argv contains the name of the inputs file entered at the command line
void AdvanceTimestepInertial(std::array< MultiFab, AMREX_SPACEDIM >& umac,
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

    BL_PROFILE_VAR("AdvanceTimestepInertial()",AdvanceTimestepInertial);

    BoxArray ba = rho_old.boxArray();
    DistributionMapping dmap = rho_old.DistributionMap();

    const Real* dx = geom.CellSize();

    Vector<Real> weights;
    weights = {1.0};

    Real theta_alpha = 1./dt;
    Real norm_pre_rhs;

    MultiFab rho_update (ba,dmap,nspecies,0);
    MultiFab bds_force;
    if (advection_type > 0) {
        bds_force.define(ba,dmap,nspecies,1);
    }
    MultiFab gmres_rhs_p(ba,dmap,       1,0);
    MultiFab dpi        (ba,dmap,       1,1);

    std::array< MultiFab, AMREX_SPACEDIM > mold;
    std::array< MultiFab, AMREX_SPACEDIM > mtemp;
    std::array< MultiFab, AMREX_SPACEDIM > adv_mom_fluxdiv;
    std::array< MultiFab, AMREX_SPACEDIM > diff_mom_fluxdiv;
    std::array< MultiFab, AMREX_SPACEDIM > stoch_mom_fluxdiv;
    std::array< MultiFab, AMREX_SPACEDIM > gmres_rhs_v;
    std::array< MultiFab, AMREX_SPACEDIM > dumac;
    std::array< MultiFab, AMREX_SPACEDIM > gradpi;
    std::array< MultiFab, AMREX_SPACEDIM > rhotot_fc_old;
    std::array< MultiFab, AMREX_SPACEDIM > rhotot_fc_new;
    std::array< MultiFab, AMREX_SPACEDIM > rho_fc;
    std::array< MultiFab, AMREX_SPACEDIM > diff_mass_flux;
    std::array< MultiFab, AMREX_SPACEDIM > total_mass_flux;
    std::array< MultiFab, AMREX_SPACEDIM > umac_tmp;


    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        mold[d]             .define(convert(ba,nodal_flag_dir[d]), dmap,        1, 1);
        mtemp[d]            .define(convert(ba,nodal_flag_dir[d]), dmap,        1, 1);
        adv_mom_fluxdiv[d]  .define(convert(ba,nodal_flag_dir[d]), dmap,        1, 0);
        diff_mom_fluxdiv[d] .define(convert(ba,nodal_flag_dir[d]), dmap,        1, 0);
        stoch_mom_fluxdiv[d].define(convert(ba,nodal_flag_dir[d]), dmap,        1, 0);
        gmres_rhs_v[d]      .define(convert(ba,nodal_flag_dir[d]), dmap,        1, 0);
        dumac[d]            .define(convert(ba,nodal_flag_dir[d]), dmap,        1, 1);
        gradpi[d]           .define(convert(ba,nodal_flag_dir[d]), dmap,        1, 0);
        rhotot_fc_old[d]    .define(convert(ba,nodal_flag_dir[d]), dmap,        1, 1);
        rhotot_fc_new[d]    .define(convert(ba,nodal_flag_dir[d]), dmap,        1, 1);
        rho_fc[d]           .define(convert(ba,nodal_flag_dir[d]), dmap, nspecies, 0);
        diff_mass_flux[d]   .define(convert(ba,nodal_flag_dir[d]), dmap, nspecies, 0);
        total_mass_flux[d]  .define(convert(ba,nodal_flag_dir[d]), dmap, nspecies, 0);

        if (advection_type > 0) {
            umac_tmp[d].define(convert(ba,nodal_flag_dir[d]), dmap,        1, 1);
        }
    }

    // only used when use_charged_fluid=T
    std::array< MultiFab, AMREX_SPACEDIM > Lorentz_force_old;
    std::array< MultiFab, AMREX_SPACEDIM > Lorentz_force_new;

    if (use_charged_fluid) {
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            Lorentz_force_old[d].define(convert(ba,nodal_flag_dir[d]), dmap, 1, 0);
            Lorentz_force_new[d].define(convert(ba,nodal_flag_dir[d]), dmap, 1, 0);
        }
    }

    // make copies of old quantities
    if (advection_type > 0) {
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            MultiFab::Copy(umac_tmp[d],umac[d],0,0,1,1);
        }
    }

    //////////////////////////////////////////////
    /// Step 1 - Calculate Predictor Diffusive and Stochastic Fluxes
    //////////////////////////////////////////////

    // diff/stoch_mass_fluxdiv already contain F_i
    // this was already done in Step 0 (initialization) or Step 6 from the previous time step

    //////////////////////////////////////////////
    // Step 2 - Predictor Euler Step
    //////////////////////////////////////////////

    if (use_charged_fluid) {
        // compute old Lorentz force
        ComputeLorentzForce(Lorentz_force_old,grad_Epot_old,permittivity,charge_old,geom);
    }

    // average rho_old and rhotot_old to faces
    AverageCCToFace(rho_old,rho_fc,0,nspecies,SPEC_BC_COMP,geom);
    AverageCCToFace(rhotot_old,rhotot_fc_old,0,1,RHO_BC_COMP,geom);

    // add D^n and St^n to rho_update
    MultiFab::Copy(rho_update,diff_mass_fluxdiv,0,0,nspecies,0);
    if (variance_coef_mass != 0.) {
        MultiFab::Add(rho_update,stoch_mass_fluxdiv,0,0,nspecies,0);
    }

    // add A^n to rho_update
    if (advection_type == 0) {
        MkAdvSFluxdiv(umac,rho_fc,rho_update,geom,0,nspecies,true);
    }
    else if (advection_type == 1 || advection_type == 2) {
        bds_force.setVal(0.);
        MultiFab::Copy(bds_force,rho_update,0,0,nspecies,0);
        bds_force.FillBoundary(geom.periodicity());

        BDS(rho_update, nspecies, SPEC_BC_COMP, rho_old, umac, bds_force, geom, dt, 2);

    }
    else {
        Abort("Invalid advection_type");
    }

    // set rho_new = rho_old + dt * (A^n + D^n + St^n)
    MultiFab::LinComb(rho_new,1.,rho_old,0,dt,rho_update,0,0,nspecies,0);

    // compute rhotot from rho in VALID REGION
    ComputeRhotot(rho_new,rhotot_new);

    // fill rho and rhotot ghost cells
    FillRhoRhototGhost(rho_new,rhotot_new,geom);

    // average rho_new and rhotot_new to faces
    AverageCCToFace(rho_new,rho_fc,0,nspecies,SPEC_BC_COMP,geom);
    AverageCCToFace(rhotot_new,rhotot_fc_new,0,1,RHO_BC_COMP,geom);

    if (use_charged_fluid) {
        // compute total charge
        DotWithZ(rho_new,charge_new);

        if (dielectric_type != 0) {
            // compute permittivity
            Abort("AdvanceTimestepInertial dielectric_type != 0");
        }
    }

    //////////////////////////////////////////////
    // Step 3 - Calculate Corrector Diffusive and Stochastic Fluxes
    // Step 4 - Predictor Crank-Nicolson Step
    //////////////////////////////////////////////

    // compute mold
    ConvertMToUmac(rhotot_fc_old,umac,mold,0);

    // build up rhs_v for gmres solve: first set gmres_rhs_v to mold/dt
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        MultiFab::Copy(gmres_rhs_v[d],mold[d],0,0,1,0);
        gmres_rhs_v[d].mult(1./dt,0);
    }

    // compute grad pi^n
    ComputeGrad(pi,gradpi,0,0,1,PRES_BC_COMP,geom);

    /*
    if (barodiffusion_type == 2) {
       // barodiffusion uses lagged grad(pi)
    }
    else if (barodiffusion_type == 3) {
       // compute p0 from rho0*g
    }
    */

    // subtract grad pi^n from gmres_rhs_v
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        MultiFab::Subtract(gmres_rhs_v[d],gradpi[d],0,0,1,0);
    }

    // compute adv_mom_fluxdiv = A^n for momentum
    MkAdvMFluxdiv(umac,mold,adv_mom_fluxdiv,dx,0);

    // add A^n for momentum into gmres_rhs_v
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        MultiFab::Add(gmres_rhs_v[d],adv_mom_fluxdiv[d],0,0,1,0);
    }

    // compute diff_mom_fluxdiv = A_0^n v^n
    MkDiffusiveMFluxdiv(diff_mom_fluxdiv,umac,eta,eta_ed,kappa,geom,dx,0);

    // add (1/2) A_0^n v^n to gmres_rhs_v
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        MultiFab::Saxpy(gmres_rhs_v[d],0.5,diff_mom_fluxdiv[d],0,0,1,0);
    }

    if (variance_coef_mom != 0.) {

        // fill the stochastic multifabs with a new set of random numbers
        sMomFlux.fillMomStochastic();

        // compute and save stoch_mom_fluxdiv = div(Sigma^n) (save for later)
        sMomFlux.StochMomFluxDiv(stoch_mom_fluxdiv,0,eta,eta_ed,Temp,Temp_ed,weights,dt);

        // add div(Sigma^n) to gmres_rhs_v
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            MultiFab::Add(gmres_rhs_v[d],stoch_mom_fluxdiv[d],0,0,1,0);
        }
    }

    // add rho^n*g to gmres_rhs_v
    bool any_grav = false;
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        if (grav[d] != 0.) any_grav = true;
    }
    if (any_grav) {
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            MultiFab::Saxpy(gmres_rhs_v[d],grav[d],rhotot_fc_old[d],0,0,1,0);
        }
    }

    // compute (eta,kappa)^{*,n+1}
    ComputeEta(rho_new, rhotot_new, eta);
    if (AMREX_SPACEDIM == 2) {
        AverageCCToNode(eta,eta_ed[0],0,1,SPEC_BC_COMP,geom);
    }
    else {
        AverageCCToEdge(eta,eta_ed,0,1,SPEC_BC_COMP,geom);
    }

    // set inhomogeneous velocity bc's to values supplied in inhomogeneous_bc_val
    //
    //

    // compute diffusive, stochastic, potential mass fluxes
    // with barodiffusion and thermodiffusion
    // this computes "-F = rho W chi [Gamma grad x... ]"
    ComputeMassFluxdiv(rho_new,rhotot_new,Temp,diff_mass_fluxdiv,stoch_mass_fluxdiv,
                       diff_mass_flux,stoch_mass_flux,sMassFlux,dt,time,geom,weights,
                       charge_new,grad_Epot_new,Epot,permittivity);

    // assemble total fluxes to be used in reservoirs
    //
    //


    // set the Dirichlet velocity value on reservoir faces
    //
    //

    if (use_charged_fluid == 1) {

        // compute new Lorentz force
        ComputeLorentzForce(Lorentz_force_new,grad_Epot_new,permittivity,charge_new,geom);

        // add (1/2) old and (1/2) new to gmres_rhs_v
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            MultiFab::Saxpy(gmres_rhs_v[d],0.5,Lorentz_force_old[d],0,0,1,0);
            MultiFab::Saxpy(gmres_rhs_v[d],0.5,Lorentz_force_new[d],0,0,1,0);
        }
    }

    // compute gmres_rhs_p
    // put "-S = div(F_i/rho_i)" into gmres_rhs_p (we will later add divu)
    gmres_rhs_p.setVal(0.);
    for (int i=0; i<nspecies; ++i) {
        MultiFab::Saxpy(gmres_rhs_p,-1/rhobar[i],diff_mass_fluxdiv,i,0,1,0);
        if (variance_coef_mass != 0.) {
            MultiFab::Saxpy(gmres_rhs_p,-1/rhobar[i],stoch_mass_fluxdiv,i,0,1,0);
        }
    }

    // modify umac to respect the boundary conditions we want after the next gmres solve
    // thus when we add A_0^n vbar^n to gmres_rhs_v and add div vbar^n to gmres_rhs_p we
    // are automatically putting the system in delta form WITH homogeneous boundary conditions
    for (int i=0; i<AMREX_SPACEDIM; ++i) {
        // set normal velocity of physical domain boundaries
        MultiFabPhysBCDomainVel(umac[i],geom,i);
        // set transverse velocity behind physical boundaries
        int is_inhomogeneous = 0;
        MultiFabPhysBCMacVel(umac[i],geom,i,is_inhomogeneous);
        // fill periodic and interior ghost cells
        umac[i].FillBoundary(geom.periodicity());
    }

    // compute mtemp = rho^{*,n+1} * vbar^n
    ConvertMToUmac(rhotot_fc_new,umac,mtemp,0);

    // subtract rho^{*,n+1} * vbar^n / dt from gmres_rhs_v
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        MultiFab::Saxpy(gmres_rhs_v[d],-1./dt,mtemp[d],0,0,1,0);
    }

    // compute mtemp = A_0^n vbar^n
    MkDiffusiveMFluxdiv(mtemp,umac,eta,eta_ed,kappa,geom,dx,0);

    // add (1/2) A_0^n vbar^n to gmres_rhs_v
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        MultiFab::Saxpy(gmres_rhs_v[d],0.5,mtemp[d],0,0,1,0);
    }

    // set physical boundary values to zero
    ZeroEdgevalPhysical(gmres_rhs_v, geom, 0, 1);

    // set rho_update to F^{*,n+1} = div(rho*chi grad c)^{*,n+1} + div(Psi^n)
    // it is used in Step 5 below
    MultiFab::Copy(rho_update,diff_mass_fluxdiv,0,0,nspecies,0);
    if (variance_coef_mass != 0.) {
        MultiFab::Add(rho_update,stoch_mass_fluxdiv,0,0,nspecies,0);
    }

    // compute div vbar^n and add to gmres_rhs_p
    // now gmres_rhs_p = div vbar^n - S^{*,n+1}
    // the sign convention is correct since we solve -div(delta v) = gmres_rhs_p
    ComputeDiv(gmres_rhs_p,umac,0,0,1,geom,1);

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

    Real gmres_abs_tol_in = gmres_abs_tol; // save this

    // This relies entirely on relative tolerance and can fail if the rhs is roundoff error only:
    // gmres_abs_tol = 0.d0 ! It is better to set gmres_abs_tol in namelist to a sensible value

    // call gmres to compute delta v and delta pi
    GMRES gmres(ba,dmap,geom);
    gmres.Solve(gmres_rhs_v, gmres_rhs_p, dumac, dpi, rhotot_fc_new, eta, eta_ed,
                kappa, theta_alpha, geom, norm_pre_rhs);

    // for the corrector gmres solve we want the stopping criteria based on the
    // norm of the preconditioned rhs from the predictor gmres solve.  otherwise
    // for cases where du in the corrector should be small the gmres stalls
    gmres_abs_tol = amrex::max(gmres_abs_tol_in, norm_pre_rhs*gmres_rel_tol);

    // restore eta and kappa
    eta.mult  (2.,0,1,1);
    kappa.mult(2.,0,1,1);
    for (int d=0; d<NUM_EDGE; ++d) {
        eta_ed[d].mult(2.,0,1,0);
    }

    // compute v^{*,n+1} = v^n + dumac
    // compute pi^{*,n+1}= pi^n + dpi
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
        int is_inhomogeneous = 1;
        MultiFabPhysBCMacVel(umac[i],geom,i,is_inhomogeneous);
        // fill periodic and interior ghost cells
        umac[i].FillBoundary(geom.periodicity());
    }

    // convert v^{*,n+1} to rho^{*,n+1}v^{*,n+1} in valid and ghost region
    // now mnew has properly filled ghost cells
    ConvertMToUmac(rhotot_fc_new,umac,mtemp,0);

    //////////////////////////////////////////////
    // Step 5 - Trapezoidal Scalar Corrector
    //////////////////////////////////////////////

    // rho_update already contains D^{*,n+1} + St^{*,n+1} for rho from above
    // add A^{*,n+1} for rho to rho_update
    if (advection_type == 0) {
        MkAdvSFluxdiv(umac,rho_fc,rho_update,geom,0,nspecies,true);

        // snew = s^{n+1}
        //      = (1/2)*(s^n + s^{*,n+1} + dt*(A^{*,n+1} + D^{*,n+1} + St^{*,n+1}))
        MultiFab::Add(rho_new,rho_old,0,0,nspecies,0);
        MultiFab::Saxpy(rho_new,dt,rho_update,0,0,nspecies,0);
        rho_new.mult(0.5,0,nspecies,0);
    }
    else if (advection_type == 1 || advection_type == 2) {

        // bds force currently contains D^n + St^n
        // add D^{*,n+1} + St^{*,n+1} and then multiply by 1/2
        MultiFab::Add(bds_force,rho_update,0,0,nspecies,0);
        bds_force.mult(0.5,0,nspecies,0);
        bds_force.FillBoundary(geom.periodicity());
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            MultiFab::Add(umac_tmp[d],umac[d],0,0,1,1);
            umac_tmp[d].mult(0.5,0,1,1);
        }
        rho_update.setVal(0.);

        BDS(rho_update, nspecies, SPEC_BC_COMP, rho_old, umac_tmp, bds_force, geom, dt, 2);
    }
    else {
        Abort("Invalid advection_type");
    }

    // need to project rho onto eos here and use this rho to compute S
    // if you do this in main_driver, the fluxes don't match the state
    // they were derived from and the Poisson solver has tolerance
    // convergence issues
    if (project_eos_int > 0 && istep%project_eos_int == 0) {
        ProjectOntoEOS(rho_new);
    }

    // compute rhotot from rho in VALID REGION
    ComputeRhotot(rho_new,rhotot_new);

    // fill rho and rhotot ghost cells
    FillRhoRhototGhost(rho_new,rhotot_new,geom);

    // average rho_new and rhotot_new to faces
    AverageCCToFace(rho_new,rho_fc,0,nspecies,SPEC_BC_COMP,geom);
    AverageCCToFace(rhotot_new,rhotot_fc_new,0,1,RHO_BC_COMP,geom);

    if (use_charged_fluid) {
        // compute total charge
        DotWithZ(rho_new,charge_new);

        if (dielectric_type != 0) {
            // compute permittivity
            Abort("AdvanceTimestepInertial dielectric_type != 0");
        }
    }

    // compute (eta,kappa)^{n+1}
    ComputeEta(rho_new, rhotot_new, eta);
    if (AMREX_SPACEDIM == 2) {
        AverageCCToNode(eta,eta_ed[0],0,1,SPEC_BC_COMP,geom);
    }
    else {
        AverageCCToEdge(eta,eta_ed,0,1,SPEC_BC_COMP,geom);
    }

    //////////////////////////////////////////////
    // Step 6 - Calculate Diffusive and Stochastic Fluxes
    // Step 7 - Corrector Crank-Nicolson Step
    //////////////////////////////////////////////

    // build up rhs_v for gmres solve: first set gmres_rhs_v to mold/dt
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        MultiFab::Copy(gmres_rhs_v[d],mold[d],0,0,1,0);
        gmres_rhs_v[d].mult(1./dt,0);
    }

    // compute grad pi^{*,n+1}
    ComputeGrad(pi,gradpi,0,0,1,PRES_BC_COMP,geom);

    /*
    if (barodiffusion_type == 2) {
       // barodiffusion uses lagged grad(pi)
    }
    else if (barodiffusion_type == 3) {
       // compute p0 from rho0*g
    }
    */

    // subtract grad pi^{*,n+1} from gmres_rhs_v
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        MultiFab::Subtract(gmres_rhs_v[d],gradpi[d],0,0,1,0);
    }


    // adv_mom_fluxdiv already contains A^n for momentum
    // add A^{*,n+1} = -rho^{*,n+1} v^{*,n+1} v^{*,n+1} for momentum to adv_mom_fluxdiv
    MkAdvMFluxdiv(umac,mtemp,adv_mom_fluxdiv,dx,1);

    // add (1/2) adv_mom_fluxdiv to gmres_rhs_v
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        MultiFab::Saxpy(gmres_rhs_v[d],0.5,adv_mom_fluxdiv[d],0,0,1,0);
    }

    // add (1/2) A_0^n v^n to gmres_rhs_v
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        MultiFab::Saxpy(gmres_rhs_v[d],0.5,diff_mom_fluxdiv[d],0,0,1,0);
    }

    if (variance_coef_mom != 0.) {

        // compute div(Sigma^n') by incrementing existing stochastic flux and
        // dividing by 2 before adding to gmres_rhs_v
        sMomFlux.StochMomFluxDiv(stoch_mom_fluxdiv,1,eta,eta_ed,Temp,Temp_ed,weights,dt);

        // divide by 2 and add the resulting div(Sigma^n') to gmres_rhs_v
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            MultiFab::Saxpy(gmres_rhs_v[d],0.5,stoch_mom_fluxdiv[d],0,0,1,0);
        }
    }

    // add gravity term
    if (any_grav) {
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            MultiFab::Saxpy(gmres_rhs_v[d],0.5*grav[d],rhotot_fc_old[d],0,0,1,0);
            MultiFab::Saxpy(gmres_rhs_v[d],0.5*grav[d],rhotot_fc_new[d],0,0,1,0);
        }
    }

    // set inhomogeneous velocity bc's to values supplied in inhomogeneous_bc_val
    //
    //
    //

    // fill the stochastic multifabs with a new set of random numbers
    if (variance_coef_mass != 0.) {
        // keep this random number engine state for checkpointing
        if (chk_int > 0) {
            Abort("AdvanceTimestepInertial.cpp checkpointing not implemented");
        }
        sMassFlux.fillMassStochastic();
    }

    // compute diffusive, stochastic, potential mass fluxes
    // with barodiffusion and thermodiffusion
    // this computes "-F = rho W chi [Gamma grad x... ]"
    ComputeMassFluxdiv(rho_new,rhotot_new,Temp,diff_mass_fluxdiv,stoch_mass_fluxdiv,
                       diff_mass_flux,stoch_mass_flux,sMassFlux,dt,time,geom,weights,
                       charge_new,grad_Epot_new,Epot,permittivity);

    // assemble total fluxes to be used in reservoirs
    //
    //


    // set the Dirichlet velocity value on reservoir faces
    //
    //

    if (use_charged_fluid == 1) {
        // compute new Lorentz force
        ComputeLorentzForce(Lorentz_force_new,grad_Epot_new,permittivity,charge_new,geom);

        // add (1/2) old and (1/2) new to gmres_rhs_v
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            MultiFab::Saxpy(gmres_rhs_v[d],0.5,Lorentz_force_old[d],0,0,1,0);
            MultiFab::Saxpy(gmres_rhs_v[d],0.5,Lorentz_force_new[d],0,0,1,0);
        }
    }

    // compute gmres_rhs_p
    // put "-S = div(F_i/rho_i)" into gmres_rhs_p (we will later add divu)
    gmres_rhs_p.setVal(0.);
    for (int i=0; i<nspecies; ++i) {
        MultiFab::Saxpy(gmres_rhs_p,-1/rhobar[i],diff_mass_fluxdiv,i,0,1,0);
        if (variance_coef_mass != 0.) {
            MultiFab::Saxpy(gmres_rhs_p,-1/rhobar[i],stoch_mass_fluxdiv,i,0,1,0);
        }
    }

    // modify umac to respect the boundary conditions we want after the next gmres solve
    // thus when we add A_0^n vbar^n to gmres_rhs_v and add div vbar^n to gmres_rhs_p we
    // are automatically putting the system in delta form WITH homogeneous boundary conditions
    for (int i=0; i<AMREX_SPACEDIM; ++i) {
        // set normal velocity of physical domain boundaries
        MultiFabPhysBCDomainVel(umac[i],geom,i);
        // set transverse velocity behind physical boundaries
        int is_inhomogeneous = 0;
        MultiFabPhysBCMacVel(umac[i],geom,i,is_inhomogeneous);
        // fill periodic and interior ghost cells
        umac[i].FillBoundary(geom.periodicity());
    }

    // compute mtemp = rho^{n+1} * vbar^{*,n+1}
    ConvertMToUmac(rhotot_fc_new,umac,mtemp,0);

    // subtract rho^{n+1} * vbar^{*,n+1} / dt from gmres_rhs_v
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        MultiFab::Saxpy(gmres_rhs_v[d],-1./dt,mtemp[d],0,0,1,0);
    }

    // set diff_mom_fluxdiv = A_0^{n+1} vbar^{n+1,*}
    MkDiffusiveMFluxdiv(diff_mom_fluxdiv,umac,eta,eta_ed,kappa,geom,dx,0);

    // add (1/2) A_0^{n+1} vbar^{n+1,*} to gmres_rhs_v
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        MultiFab::Saxpy(gmres_rhs_v[d],0.5,diff_mom_fluxdiv[d],0,0,1,0);
    }

    // set physical boundary values to zero
    ZeroEdgevalPhysical(gmres_rhs_v, geom, 0, 1);

    // compute div(vbar^{n+1,*}) and add to gmres_rhs_p
    // now gmres_rhs_p = div(vbar^{n+1,*}) - S^{n+1}
    // the sign convention is correct since we solve -div(delta v) = gmres_rhs_p
    ComputeDiv(gmres_rhs_p,umac,0,0,1,geom,1);

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

    // call gmres to compute delta v and delta pi
    gmres.Solve(gmres_rhs_v, gmres_rhs_p, dumac, dpi, rhotot_fc_new, eta, eta_ed,
                kappa, theta_alpha, geom, norm_pre_rhs);

    gmres_abs_tol = gmres_abs_tol_in; // Restore the desired tolerance

    // restore eta and kappa
    eta.mult  (2.,0,1,1);
    kappa.mult(2.,0,1,1);
    for (int d=0; d<NUM_EDGE; ++d) {
        eta_ed[d].mult(2.,0,1,0);
    }

    // compute v^{n+1} = v^{n+1,*} + dumac
    // compute pi^{n+1}= pi^{n+1,*} + dpi
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
        int is_inhomogeneous = 1;
        MultiFabPhysBCMacVel(umac[i],geom,i,is_inhomogeneous);
        // fill periodic and interior ghost cells
        umac[i].FillBoundary(geom.periodicity());
    }

    //////////////////////////////////////////////
    // End Time-Advancement
    //////////////////////////////////////////////
}