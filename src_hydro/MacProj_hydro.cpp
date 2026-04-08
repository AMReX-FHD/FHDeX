#include "hydro_functions.H"

#include "common_functions.H"

#include <AMReX_BoxArray.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_Vector.H>

#include <AMReX_MLABecLaplacian.H>
#include <AMReX_MLMG.H>


// umac enters with face-centered, time-centered Utilde^* and should leave with Utilde
// macphi is the solution to the elliptic solve

void
MacProj_hydro (std::array< MultiFab, AMREX_SPACEDIM >& umac,
               const MultiFab& rho,
               const Geometry& geom,
               const bool& full_solve)
{
    // timer for profiling
    BL_PROFILE_VAR("MacProj_hydro()",MacProj_hydro);

    BoxArray grids = rho.boxArray();
    DistributionMapping dmap = rho.DistributionMap();

    MultiFab solverrhs; // this will hold solver RHS = macrhs - div(umac)
    solverrhs.define(grids, dmap, 1, 0);

    MultiFab macphi(grids,dmap,1,1);
    MultiFab macrhs(grids,dmap,1,1);
    macrhs.setVal(0.0);
    macphi.setVal(0.0);

    // compute the RHS for the solve, RHS = macrhs - div(umac)
    ComputeMACSolverRHS(solverrhs,macrhs,umac,geom);

    // coefficients for solver
    MultiFab acoef;
    std::array< MultiFab, AMREX_SPACEDIM > face_bcoef;
    acoef.define(grids, dmap, 1, 0);
    AMREX_D_TERM(face_bcoef[0].define(convert(grids,nodal_flag_x), dmap, 1, 0);,
                 face_bcoef[1].define(convert(grids,nodal_flag_y), dmap, 1, 0);,
                 face_bcoef[2].define(convert(grids,nodal_flag_z), dmap, 1, 0););

    // set cell-centered A coefficient to zero
    acoef.setVal(0.);

    // 1) average face-centered B coefficients to rho
    AverageCCToFace(rho, face_bcoef, 0, 1, RHO_BC_COMP, geom);

    // 2) invert B coefficients to 1/rho
    for (int idim=0; idim<AMREX_SPACEDIM; ++idim) {
        face_bcoef[idim].invert(1.0,0,1);
    }

    //
    // Set up implicit solve using MLABecLaplacian class
    //
    LPInfo info;
    MLABecLaplacian mlabec({geom}, {grids}, {dmap}, info);

    // order of stencil
    int linop_maxorder = 2;
    mlabec.setMaxOrder(linop_maxorder);

    // set boundaries for mlabec using velocity bc's
    SetMacSolverBCs(mlabec);

    mlabec.setLevelBC(0, &macphi);
    mlabec.setScalars(0.0, 1.0);
    mlabec.setACoeffs(0, acoef);
    mlabec.setBCoeffs(0, amrex::GetArrOfConstPtrs(face_bcoef));

    // solve -div B grad phi = RHS

    // build an MLMG solver
    MLMG mac_mlmg(mlabec);

    // set solver parameters
    mac_mlmg.setVerbose(mg_verbose);
    mac_mlmg.setBottomVerbose(cg_verbose);

    // for the preconditioner, we do 1 v-cycle and the bottom solver is smooths
    if (!full_solve) {
        mac_mlmg.setBottomSolver(amrex::MLMG::BottomSolver::smoother);
        mac_mlmg.setFixedIter(1);
        mac_mlmg.setBottomSmooth(8);
    }

    // solve for phi
    mac_mlmg.solve({&macphi}, {&solverrhs}, mg_rel_tol, mg_abs_tol);

    // update velocity, Utilde = Utilde^* - B grad phi

    // storage for "-B grad_phi"
    std::array<MultiFab,AMREX_SPACEDIM> mac_fluxes;
    AMREX_D_TERM(mac_fluxes[0].define(convert(grids,nodal_flag_x), dmap, 1, 0);,
                 mac_fluxes[1].define(convert(grids,nodal_flag_y), dmap, 1, 0);,
                 mac_fluxes[2].define(convert(grids,nodal_flag_z), dmap, 1, 0););

    std::array<MultiFab*,AMREX_SPACEDIM> mac_fluxptr;
    // fluxes computed are "-B grad phi"
    mac_fluxptr = GetArrOfPtrs(mac_fluxes);
    mac_mlmg.getFluxes({mac_fluxptr});

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        // add -B grad phi to Utilde
        MultiFab::Add(umac[idim], mac_fluxes[idim], 0, 0, 1, 0);
    }

    // fill periodic ghost cells
    for (int d = 0; d < AMREX_SPACEDIM; d++) {
        // Do apply BCs so that all ghost cells are filled
        MultiFabPhysBCDomainVel(umac[d], geom, d);
        int is_inhomogeneous = 1;
        MultiFabPhysBCMacVel(umac[d], geom, d, is_inhomogeneous);
        umac[d].FillBoundary(geom.periodicity());
    }
}

// compute the RHS for the solve, solverrhs = macrhs - div(umac)
void ComputeMACSolverRHS (MultiFab& solverrhs,
                          const MultiFab& macrhs,
                          const std::array< MultiFab, AMREX_SPACEDIM >& umac,
                          const Geometry& geom)
{
    // timer for profiling
    BL_PROFILE_VAR("ComputeMACSolverRHS()",ComputeMACSolverRHS);

    const GpuArray<Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();

    for ( MFIter mfi(solverrhs); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.validbox();

        AMREX_D_TERM(Array4<Real const> const& umac_fab = umac[0].array(mfi);,
                     Array4<Real const> const& vmac_fab = umac[1].array(mfi);,
                     Array4<Real const> const& wmac_fab = umac[2].array(mfi););


        Array4<Real>       const& solverrhs_fab = solverrhs.array(mfi);
        Array4<Real const> const&    macrhs_fab =    macrhs.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            solverrhs_fab(i,j,k) = macrhs_fab(i,j,k) -
                AMREX_D_TERM(   ( umac_fab(i+1,j,k) - umac_fab(i,j,k) ) / dx[0],
                              - ( vmac_fab(i,j+1,k) - vmac_fab(i,j,k) ) / dx[1],
                              - ( wmac_fab(i,j,k+1) - wmac_fab(i,j,k) ) / dx[2] );;
        });
    }
}

// Set boundaries for MAC velocities
void SetMacSolverBCs(MLABecLaplacian& mlabec)
{
    // timer for profiling
    BL_PROFILE_VAR("SetMacSolverBCs()", SetMacSolverBCs);

    // build array of boundary conditions needed by MLABecLaplacian
    std::array<LinOpBCType, AMREX_SPACEDIM> lo_mlmg_bc;
    std::array<LinOpBCType, AMREX_SPACEDIM> hi_mlmg_bc;

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {
        if (DefaultGeometry().isPeriodic(idim)) {
            lo_mlmg_bc[idim] = hi_mlmg_bc[idim] = LinOpBCType::Periodic;
        } else {
            //amrex::Error("Invalid BC");
            lo_mlmg_bc[idim] = hi_mlmg_bc[idim] = LinOpBCType::Neumann;
            Print() << "Warning! non-periodic boundary conditions in MacProj_hydro." << std::endl
                    << " => Assuming Neumann. But preconditioner might not work properly." << std::endl;
        }
    }

    mlabec.setDomainBC(lo_mlmg_bc,hi_mlmg_bc);
}