#include "common_functions.H"

#include "gmres_functions.H"

#include <AMReX_MLMG.H>
#include <AMReX_MLABecLaplacian.H>
#include <AMReX_MLLinOp.H>


void MacProj(const std::array<MultiFab, AMREX_SPACEDIM>& alphainv_fc,
             MultiFab& mac_rhs,
             MultiFab& phi,
             const Geometry& geom,
             bool full_solve)
{
    BL_PROFILE_VAR("MacProj()",MacProj);

    int lev=0;

    BoxArray ba = phi.boxArray();
    DistributionMapping dmap = phi.DistributionMap();
    LPInfo info;

    MLABecLaplacian mlabec({geom}, {ba}, {dmap}, info);

    // order of stencil
    int stencil_order = 2;
    mlabec.setMaxOrder(stencil_order);

    mlabec.setScalars(0.0, -1.0);

    // build array of boundary conditions needed by MLABecLaplacian
    std::array<LinOpBCType,AMREX_SPACEDIM> lo_mlmg_bc;
    std::array<LinOpBCType,AMREX_SPACEDIM> hi_mlmg_bc;

    for (int i=0; i<AMREX_SPACEDIM; ++i) {
        if (geom.isPeriodic(i)) {
            lo_mlmg_bc[i] = hi_mlmg_bc[i] = LinOpBCType::Periodic;
        } else {
            lo_mlmg_bc[i] = hi_mlmg_bc[i] = LinOpBCType::Neumann;
        }
    }
    mlabec.setDomainBC(lo_mlmg_bc,hi_mlmg_bc);
    mlabec.setLevelBC(lev, &phi);

    // coefficients for solver (alpha already set to zero via setScalars)
    mlabec.setBCoeffs(lev,amrex::GetArrOfConstPtrs(alphainv_fc));

    MLMG mlmg(mlabec);

    mlmg.setVerbose(mg_verbose);

    // for the preconditioner, we do 1 v-cycle and the bottom solver is smooths
    if (!full_solve) {
        mlmg.setBottomSolver(amrex::MLMG::BottomSolver::smoother);
        mlmg.setFixedIter(1);
        mlmg.setBottomSmooth(8);
    }

    mlmg.solve({&phi}, {&mac_rhs}, mg_rel_tol, mg_abs_tol);

    phi.FillBoundary(geom.periodicity());
}

void SubtractWeightedGradP(std::array<MultiFab, AMREX_SPACEDIM>& x_u,
                           const std::array<MultiFab, AMREX_SPACEDIM>& alphainv_fc,
                           MultiFab& phi,
                           std::array<MultiFab, AMREX_SPACEDIM>& gradp,
                           const Geometry& geom)
{
    BL_PROFILE_VAR("SubtractWeightedGradP()",SubtractWeightedGradP);
    
    BoxArray ba = phi.boxArray();
    DistributionMapping dmap = phi.DistributionMap();

    ComputeGrad(phi,gradp,0,0,1,0,geom);

    for (int i=0; i<AMREX_SPACEDIM; ++i) {
        MultiFab::Multiply(gradp[i],alphainv_fc[i],0,0,1,0);
        MultiFab::Saxpy(x_u[i],-1.,gradp[i],0,0,1,0);
    }
}

void CCApplyNegLap(MultiFab& phi,
                   MultiFab& Lphi,
                   const std::array<MultiFab, AMREX_SPACEDIM>& beta_fc,
                   const Geometry& geom)
{
    BL_PROFILE_VAR("CCApplyOp()",CCApplyOp);
    
    int lev=0;

    BoxArray ba = phi.boxArray();
    DistributionMapping dmap = phi.DistributionMap();
    LPInfo info;
    info.setMaxCoarseningLevel(0); // turn off mg coarsening since no actual solve is performed

    MLABecLaplacian mlabec({geom}, {ba}, {dmap}, info);

    // order of stencil
    int stencil_order = 2;
    mlabec.setMaxOrder(stencil_order);

    mlabec.setScalars(0.0, 1.0);

    // build array of boundary conditions needed by MLABecLaplacian
    std::array<LinOpBCType,AMREX_SPACEDIM> lo_mlmg_bc;
    std::array<LinOpBCType,AMREX_SPACEDIM> hi_mlmg_bc;

    for (int i=0; i<AMREX_SPACEDIM; ++i) {
        if (geom.isPeriodic(i)) {
            lo_mlmg_bc[i] = hi_mlmg_bc[i] = LinOpBCType::Periodic;
        }
        else {
            Abort("CCApplyNegLap only works for periodic");
        }
    }
    mlabec.setDomainBC(lo_mlmg_bc,hi_mlmg_bc);
    mlabec.setLevelBC(lev, &phi);

    // coefficients for solver
    mlabec.setBCoeffs(lev,amrex::GetArrOfConstPtrs(beta_fc));

    MLMG mlmg(mlabec);

    mlmg.apply({&Lphi},{&phi});
}
