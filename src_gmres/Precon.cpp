#include "gmres_functions.H"

Precon::Precon() {}


void Precon::Define(const Vector<BoxArray>& ba_in,
              const Vector<DistributionMapping>& dmap_in,
              const Vector<Geometry>& geom_in) {

    nlevels = ba_in.size();
    phi.resize(nlevels);
    mac_rhs.resize(nlevels);
    gradp.resize(nlevels);
    
    for(int lev=0;lev<nlevels;++lev)
    {
        phi[lev].define    (ba_in[lev],dmap_in[lev],1,1);
        mac_rhs[lev].define(ba_in[lev],dmap_in[lev],1,0);

        for (int d=0; d<AMREX_SPACEDIM; d++) {
            gradp[lev][d].define(convert(ba_in[lev], nodal_flag_dir[d]), dmap_in[lev], 1, 0);
        }
    }
    
//    macproj.Define(ba_in,dmap_in,geom_in, nlevels);    
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

    MultiFab cc_mask;
    std::array< MultiFab, AMREX_SPACEDIM > fc_mask;
    
    cc_mask.define(gamma.boxArray(),gamma.DistributionMap(), 1, 0);
    cc_mask.setVal(1);
    for(int d=0; d<AMREX_SPACEDIM; ++d) {

        fc_mask[d].define(convert(x_u[d].boxArray(),nodal_flag_dir[d]), x_u[d].DistributionMap(), 1, 0);
        fc_mask[d].setVal(1);       
    }
    
    MultiFab* cc_maskp = &cc_mask;
    std::array<MultiFab, AMREX_SPACEDIM>* fc_maskp = &fc_mask;     
        
    Apply(b_up,b_pp,x_up,x_pp,alpha_fcp,alphainv_fcp,betap,beta_edp,gammap,cc_maskp,fc_maskp,theta_alpha,geomp, StagSolver);


}

void Precon::Apply(std::array<MultiFab, AMREX_SPACEDIM>* & b_u,
                   MultiFab*                             & b_p,
                   Vector<std::array<MultiFab, AMREX_SPACEDIM>> & x_u,
                   Vector<MultiFab>                             & x_p,
                   std::array<MultiFab, AMREX_SPACEDIM>* & alpha_fc,
                   Vector<std::array<MultiFab, AMREX_SPACEDIM>> & alphainv_fc,
                   MultiFab*                             & beta, 
                   std::array<MultiFab, NUM_EDGE>*      & beta_ed,
                   MultiFab*                            & gamma,
                   MultiFab*                            & cc_mask,
                   std::array<MultiFab, AMREX_SPACEDIM>* & fc_mask,                   
                   const Real                           & theta_alpha,
                   Geometry*                            & geom,
                   StagMGSolver& StagSolver)
{

    std::array<MultiFab, AMREX_SPACEDIM>* x_up = &x_u[0];
    std::array<MultiFab, AMREX_SPACEDIM>* alphainv_fcp = &alphainv_fc[0];
    MultiFab* x_pp = &x_p[0];    
        
    Apply(b_u,b_p,x_up,x_pp,alpha_fc,alphainv_fcp,beta,beta_ed,gamma,cc_mask, fc_mask, theta_alpha,geom,StagSolver);


}

void Precon::Apply(Vector<std::array<MultiFab, AMREX_SPACEDIM>> & b_u,
                   Vector<MultiFab>                             & b_p,
                   Vector<std::array<MultiFab, AMREX_SPACEDIM>> & x_u,
                   Vector<MultiFab>                             & x_p,
                   std::array<MultiFab, AMREX_SPACEDIM>*        & alpha_fc,
                   Vector<std::array<MultiFab, AMREX_SPACEDIM>> & alphainv_fc,
                   MultiFab*                                    & beta, 
                   std::array<MultiFab, NUM_EDGE>*              & beta_ed,
                   MultiFab*                                    & gamma,
                   MultiFab*                                    & cc_mask,
                   std::array<MultiFab, AMREX_SPACEDIM>*        & fc_mask,                   
                   const Real                                   & theta_alpha,
                   Geometry*                                    & geom,
                   StagMGSolver& StagSolver)
{

    std::array<MultiFab, AMREX_SPACEDIM>* b_up = &b_u[0];
    std::array<MultiFab, AMREX_SPACEDIM>* x_up = &x_u[0];
    std::array<MultiFab, AMREX_SPACEDIM>* alphainv_fcp = &alphainv_fc[0];
    MultiFab* b_pp = &b_p[0];
    MultiFab* x_pp = &x_p[0];        
        
    Apply(b_up,b_pp,x_up,x_pp,alpha_fc,alphainv_fcp,beta,beta_ed,gamma,cc_mask, fc_mask, theta_alpha,geom,StagSolver);


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
                   MultiFab*                            & cc_mask,
                   std::array<MultiFab, AMREX_SPACEDIM>* & fc_mask,                   
                   const Real                           & theta_alpha,
                   Geometry*                            & geom,
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
        

        if(nlevels>1)
        {

           FaceFillCoarse(b_u,1);
           FaceFillGhost(b_u,geom,1);    
           StagSolver.TopSolve(alpha_fc,beta,beta_ed,gamma,x_u,b_u,geom,theta_alpha,0);
        }
        StagSolver.Solve(alpha_fc[0],beta[0],beta_ed[0],gamma[0],x_u[0],b_u[0],theta_alpha);
        if(nlevels>1)
        {
           StagSolver.TopSolve(alpha_fc,beta,beta_ed,gamma,x_u,b_u,geom,theta_alpha,1);
           FaceFillCoarse(x_u,0);
           FaceFillGhost(x_u,geom,0);
        }

        ////////////////////
        // STEP 2: Construct RHS for pressure Poisson problem
        ////////////////////
      

        // set mac_rhs = D(x_u^star)
        for(int lev=0;lev<nlevels;++lev)
        {
            ComputeDiv(mac_rhs[lev],x_u[lev],0,0,1,geom[lev],0);

        // add b_p to mac_rhs
            MultiFab::Add(mac_rhs[lev],b_p[lev],0,0,1,0);
        }
       

        ////////////////////
        // STEP 3: Compute x_u
        ////////////////////

        // use multigrid to solve for Phi
        // x_u^star is only passed in to get a norm for absolute residual criteria
        macproj.Solve(alphainv_fc[0],mac_rhs[0],phi[0],geom[0]);
        if(nlevels>1)
        {
            CellFillFine(phi, geom);
            CellFillGhost(phi,geom);    
        }
        // x_u = x_u^star - (alpha I)^-1 grad Phi
        for(int lev=0;lev<nlevels;++lev)
        {
            SubtractWeightedGradP(x_u[lev],alphainv_fc[lev],phi[lev],gradp[lev],geom[lev]);
        }

        ////////////////////
        // STEP 4: Compute x_p by applying the Schur complement approximation
        ////////////////////

        if (visc_schur_approx == 0) {
            // if precon_type = +1, or theta_alpha=0 then x_p = theta_alpha*Phi - c*beta*(mac_rhs)
            // if precon_type = -1                   then x_p = theta_alpha*Phi - c*beta*L_alpha Phi

            if (precon_type == 1 || theta_alpha == 0) {
                // first set x_p = -mac_rhs
                for(int lev=0;lev<nlevels;++lev)
                {
                    MultiFab::Copy(x_p[lev],mac_rhs[lev],0,0,1,0);
                    x_p[lev].mult(-1.,0,1,0);
                }
            }
            else {
                // first set x_p = -L_alpha Phi
                for(int lev=0;lev<nlevels;++lev)
                {
                    CCApplyNegLap(phi[lev],x_p[lev],alphainv_fc[lev],geom[lev]);
                }
            }

            if ( amrex::Math::abs(visc_type) == 1 || amrex::Math::abs(visc_type) == 2) {
                // multiply x_p by beta; x_p = -beta L_alpha Phi
                for(int lev=0;lev<nlevels;++lev)
                {
                    MultiFab::Multiply(x_p[lev],beta[lev],0,0,1,0);
                }

                if (amrex::Math::abs(visc_type) == 2) {
                    // multiply by c=2; x_p = -2*beta L_alpha Phi
                    for(int lev=0;lev<nlevels;++lev)
                    {
                        x_p[lev].mult(2.,0,1,0);
                    }
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


            for(int lev=0;lev<nlevels;++lev)
            {   
                // multiply Phi by theta_alpha            
                phi[lev].mult(theta_alpha,0,1,0);
                // add theta_alpha*Phi to x_p
                MultiFab::Add(x_p[lev],phi[lev],0,0,1,0);                
            }

        }
        else {
            Abort("StagApplyOp: visc_schur_approx != 0 not supported");
        }
    }
    else {
        Abort("StagApplyOp: unsupposed precon_type");
    }
    
    if(nlevels>1)
    {
        //FaceFillCoarse(x_u,0);
        CellFillGhost(x_p,geom);    
    }

    ////////////////////
    // STEP 5: Handle null-space issues in MG solvers
    ////////////////////

    // subtract off mean value: Single level only! No need for ghost cells
    if(nlevels==1)
    {
        SumStag(x_u[0],mean_val_umac,true);
        SumCC(x_p[0],0,mean_val_pres,true);
    }    
    else
    {
        for(int d=0;d<AMREX_SPACEDIM;++d)
        {
            MultiFab::Multiply(x_u[0][d],fc_mask[0][d],0,0,1,0);
        }
        
        MultiFab::Multiply(x_p[0],cc_mask[0],0,0,1,0);
        
        SumStag(x_u[0],mean_val_umac,false);
        SumCC(x_p[0],0,mean_val_pres,false);
        
        Vector<Real> mean_val_umac_t(AMREX_SPACEDIM);
        Real         mean_val_pres_t;
        
        SumStag(x_u[1],mean_val_umac_t,false);
        SumCC(x_p[1],0,mean_val_pres_t,false);
        
        BoxArray ba_xuc = x_u[0][0].boxArray();
        ba_xuc.enclosedCells();
        long numpts_xuc = ba_xuc.numPts();
        BoxArray ba_xuf = x_u[1][0].boxArray();
        ba_xuf.enclosedCells();
        long numpts_xuf = ba_xuf.numPts();  
        
        for(int d=0;d<AMREX_SPACEDIM;++d)
        {
            mean_val_umac[d] = (numpts_xuc*mean_val_umac[d]+numpts_xuf*mean_val_umac_t[d])/(numpts_xuc+numpts_xuf);
        }
        
        BoxArray ba_xpc = x_p[0].boxArray();
        ba_xpc.enclosedCells();
        long numpts_xpc = ba_xpc.numPts();
        BoxArray ba_xpf = x_p[1].boxArray();
        ba_xpf.enclosedCells();
        long numpts_xpf = ba_xpf.numPts();
        
        mean_val_pres = (numpts_xpc*mean_val_pres + numpts_xpf*mean_val_pres_t)/(numpts_xpc+numpts_xpf);
        
        FaceFillCoarse(x_u,0);
        CellFillCoarse(x_p,geom);
    
    }

    for(int lev=0;lev<nlevels;++lev)
    {
        // The pressure Poisson problem is always singular:
        x_p[lev].plus(-mean_val_pres,0,1,0);

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
                    if (geom[lev].isPeriodic(d)) {
                        x_u[lev][d].plus(-mean_val_umac[d],0,1,0);
                    }
                }
            }
        }
    }
}
