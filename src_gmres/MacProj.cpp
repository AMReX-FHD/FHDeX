#include "common_functions.H"
#include "gmres_functions.H"
#include "MacProj.H"

MacProj::MacProj() {}

void MacProj::Define(const BoxArray& ba,
                     const DistributionMapping& dmap,
                     const Geometry& geom) {

    LPInfo info;
    mlabec.define({geom}, {ba}, {dmap}, info);

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
}


void MacProj::Solve(const std::array<MultiFab, AMREX_SPACEDIM>& alphainv_fc,
                    MultiFab& mac_rhs,
                    MultiFab& phi,
                    const Geometry& geom,
                    bool full_solve)
{
    BL_PROFILE_VAR("MacProj()",MacProj);

    int lev=0;
    mlabec.setLevelBC(lev, &phi);

    // coefficients for solver (alpha already set to zero via setScalars)
    mlabec.setBCoeffs(lev,amrex::GetArrOfConstPtrs(alphainv_fc));

    MLMG mlmg(mlabec);

    mlmg.setVerbose(mg_verbose);
    mlmg.setBottomVerbose(cg_verbose);

    // for the preconditioner, we do 1 v-cycle and the bottom solver is smooths
    if (!full_solve) {
        if (mg_bottom_solver == 0) {
            mlmg.setBottomSolver(amrex::MLMG::BottomSolver::smoother);
        }
        else if (mg_bottom_solver == 1) {
            mlmg.setBottomSolver(amrex::MLMG::BottomSolver::bicgstab);
        }
        else {
            Abort("MacProj.cpp: only mg_bottom_solver=0");
        }
        mlmg.setFixedIter(mg_max_vcycles);
        mlmg.setPreSmooth(mg_nsmooths_down);
        mlmg.setPostSmooth(mg_nsmooths_up);
        mlmg.setFinalSmooth(mg_nsmooths_bottom);
    }

    mlmg.solve({&phi}, {&mac_rhs}, mg_rel_tol, mg_abs_tol);

    phi.FillBoundary(geom.periodicity());
}

