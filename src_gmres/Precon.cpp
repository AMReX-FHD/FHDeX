#include "gmres_functions.H"

Precon::Precon() {}


void Precon::Define(const Vector<BoxArray>& ba_in,
              const Vector<DistributionMapping>& dmap_in,
              const Vector<Geometry>& geom_in, int nlev) {

    nlevels = nlev;
    phi.resize(2);
    mac_rhs.resize(2);
    gradp.resize(2);
    
    for(int lev=0;lev<nlevels;++lev)
    {
        phi[lev].define    (ba_in[lev],dmap_in[lev],1,1);
        mac_rhs[lev].define(ba_in[lev],dmap_in[lev],1,0);

        for (int d=0; d<AMREX_SPACEDIM; d++) {
            gradp[lev][d].define(convert(ba_in[lev], nodal_flag_dir[d]), dmap_in[lev], 1, 0);
        }
    }
    
    macproj.Define(ba_in[0],dmap_in[0],geom_in[0]);    
   
}

void Precon::Define(const BoxArray& ba_in,
                    const DistributionMapping& dmap_in,
                    const Geometry& geom_in) {

    BL_PROFILE_VAR("Precon::Define()",Precon);
    
    nlevels = 1;
    phi.resize(1);
    mac_rhs.resize(1);
    gradp.resize(1);
    
    phi[0].define    (ba_in,dmap_in,1,1);
    mac_rhs[0].define(ba_in,dmap_in,1,0);

    for (int d=0; d<AMREX_SPACEDIM; d++) {
        gradp[0][d].define(convert(ba_in, nodal_flag_dir[d]), dmap_in, 1, 0);
    }

    macproj.Define(ba_in,dmap_in,geom_in);


}

void Precon::Apply(std::array<MultiFab, AMREX_SPACEDIM> & b_u,
                   MultiFab & b_p,
                   std::array<MultiFab, AMREX_SPACEDIM> & x_u,
                   MultiFab & x_p,
                   std::array<MultiFab, AMREX_SPACEDIM> & alpha_fc,
                   std::array<MultiFab, AMREX_SPACEDIM> & alphainv_fc,
                   MultiFab & beta, std::array<MultiFab, NUM_EDGE> & beta_ed,
                   MultiFab & gamma,
                   Real & theta_alpha,
                   Geometry & geom,
                   StagMGSolver& StagSolver)
{

    std::array<MultiFab, AMREX_SPACEDIM>* b_up = &b_u;
    std::array<MultiFab, AMREX_SPACEDIM>* x_up = &x_u;
    std::array<MultiFab, NUM_EDGE>* beta_edp = &beta_ed;    
    std::array<MultiFab, AMREX_SPACEDIM>* alpha_fcp = &alpha_fc;
    std::array<MultiFab, AMREX_SPACEDIM>* alphainv_fcp = &alphainv_fc;
         
    MultiFab* b_pp = &b_p;
    MultiFab* x_pp = &x_p;
    MultiFab* gammap = &gamma;   
    MultiFab* betap = &beta;             
    Geometry* geomp = &geom;     
        
    Apply(b_up,b_pp,x_up,x_pp,alpha_fcp,alphainv_fcp,betap,beta_edp,gammap,theta_alpha,geomp,1, StagSolver);


}


void Precon::Apply(std::array<MultiFab, AMREX_SPACEDIM>* & b_u,
                   MultiFab*                             & b_p,
                   std::array<MultiFab, AMREX_SPACEDIM>* & x_u,
                   MultiFab*                             & x_p,
                   std::array<MultiFab, AMREX_SPACEDIM>* & alpha_fc,
                   std::array<MultiFab, AMREX_SPACEDIM>* & alphainv_fc,
                   MultiFab*                             & beta, 
                   std::array<MultiFab, NUM_EDGE>*      & beta_ed,
                   MultiFab*                            & gamma,
                   const Real                           & theta_alpha,
                   Geometry*                            & geom,
                   int                                    nlevels,
                   StagMGSolver& StagSolver)
{

    BL_PROFILE_VAR("Precon::Apply()",Precon_Apply);

    Vector<BoxArray> ba(nlevels);
    Vector<DistributionMapping> dmap(nlevels);
    
    for(int lev=0;lev<nlevels;lev++)
    {
        ba[lev] = b_p[lev].boxArray();
        dmap[lev] = b_p[lev].DistributionMap();
        
        phi[lev].setVal(0.);
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            x_u[lev][d].setVal(0.);
        }
    }
     
    Real         mean_val_pres;
    Vector<Real> mean_val_umac(AMREX_SPACEDIM);

    // set the initial guess for Phi in the Poisson solve to 0
    // set x_u = 0 as initial guess


    // 1 = projection preconditioner
    // 2 = lower triangular preconditioner
    // 3 = upper triangular preconditioner
    // 4 = block diagonal preconditioner
    // 5 = Uzawa-type approximation (see paper)
    // 6 = upper triangular + viscosity-based BFBt Schur complement (from Georg Stadler)

    // projection preconditioner
    if (amrex::Math::abs(precon_type) == 1) {

        ////////////////////
        // STEP 1: Solve for an intermediate state, x_u^star, using an implicit viscous solve
        //         x_u^star = A^{-1} b_u
        ////////////////////

        // x_u^star = A^{-1} b_u
        StagSolver.Solve(alpha_fc[0],beta[0],beta_ed[0],gamma[0],x_u[0],b_u[0],theta_alpha);

        ////////////////////
        // STEP 2: Construct RHS for pressure Poisson problem
        ////////////////////

        // set mac_rhs = D(x_u^star)
        ComputeDiv(mac_rhs[0],x_u[0],0,0,1,geom[0],0);

        // add b_p to mac_rhs
        MultiFab::Add(mac_rhs[0],b_p[0],0,0,1,0);

        ////////////////////
        // STEP 3: Compute x_u
        ////////////////////

        // use multigrid to solve for Phi
        // x_u^star is only passed in to get a norm for absolute residual criteria
        macproj.Solve(alphainv_fc[0],mac_rhs[0],phi[0],geom[0]);

        // x_u = x_u^star - (alpha I)^-1 grad Phi
        SubtractWeightedGradP(x_u[0],alphainv_fc[0],phi[0],gradp[0],geom[0]);

        ////////////////////
        // STEP 4: Compute x_p by applying the Schur complement approximation
        ////////////////////

        if (visc_schur_approx == 0) {
            // if precon_type = +1, or theta_alpha=0 then x_p = theta_alpha*Phi - c*beta*(mac_rhs)
            // if precon_type = -1                   then x_p = theta_alpha*Phi - c*beta*L_alpha Phi

            if (precon_type == 1 || theta_alpha == 0) {
                // first set x_p = -mac_rhs
                MultiFab::Copy(x_p[0],mac_rhs[0],0,0,1,0);
                x_p[0].mult(-1.,0,1,0);
            }
            else {
                // first set x_p = -L_alpha Phi
                CCApplyNegLap(phi[0],x_p[0],alphainv_fc[0],geom[0]);
            }

            if ( amrex::Math::abs(visc_type) == 1 || amrex::Math::abs(visc_type) == 2) {
                // multiply x_p by beta; x_p = -beta L_alpha Phi
                MultiFab::Multiply(x_p[0],beta[0],0,0,1,0);

                if (amrex::Math::abs(visc_type) == 2) {
                    // multiply by c=2; x_p = -2*beta L_alpha Phi
                    x_p[0].mult(2.,0,1,0);
                }
            }
            else if (amrex::Math::abs(visc_type) == 3) {

                // multiply x_p by gamma, use mac_rhs a temparary to save x_p
                MultiFab::Copy(mac_rhs[0],x_p[0],0,0,1,0);
                MultiFab::Multiply(mac_rhs[0],gamma[0],0,0,1,0);
                // multiply x_p by beta; x_p = -beta L_alpha Phi
                MultiFab::Multiply(x_p[0],beta[0],0,0,1,0);
                // multiply by c=4/3; x_p = -(4/3) beta L_alpha Phi
                x_p[0].mult(4./3.,0,1,0);
                // x_p = -(4/3) beta L_alpha Phi - gamma L_alpha Phi
                MultiFab::Add(x_p[0],mac_rhs[0],0,0,1,0);
            }

            // multiply Phi by theta_alpha
            phi[0].mult(theta_alpha,0,1,0);

            // add theta_alpha*Phi to x_p
            MultiFab::Add(x_p[0],phi[0],0,0,1,0);
        }
        else {
            Abort("StagApplyOp: visc_schur_approx != 0 not supported");
        }
    }
    else {
        Abort("StagApplyOp: unsupposed precon_type");
    }

    ////////////////////
    // STEP 5: Handle null-space issues in MG solvers
    ////////////////////

    // subtract off mean value: Single level only! No need for ghost cells
    SumStag(x_u[0],mean_val_umac,true);
    SumCC(x_p[0],0,mean_val_pres,true);

    // The pressure Poisson problem is always singular:
    x_p[0].plus(-mean_val_pres,0,1,0);

    // The velocity problem is also singular under these cases
    if (theta_alpha == 0.) {

        bool no_wall_is_no_slip = true;

        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            if (bc_vel_lo[d] == 2 || bc_vel_hi[d] == 2) {
                no_wall_is_no_slip = false;
            }
        }


        if (no_wall_is_no_slip) {
            for (int d=0; d<AMREX_SPACEDIM; ++d) {
                if (geom[0].isPeriodic(d)) {
                    x_u[0][d].plus(-mean_val_umac[d],0,1,0);
                }
            }
        }
    }
}
