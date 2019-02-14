#include "hydro_functions.H"
#include "hydro_functions_F.H"

#include "common_functions.H"
#include "common_functions_F.H"
#include "common_namespace.H"

#include "gmres_namespace.H"

#include <AMReX_BoxArray.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_Vector.H>

#include <AMReX_FluxRegister.H>
#include <AMReX_MLABecLaplacian.H>
#include <AMReX_MLMG.H>

using namespace common;
using namespace gmres;

// ADAPTED FROM MAESTROeX
// umac enters with face-centered, time-centered Utilde^* and should leave with Utilde
// macphi is the solution to the elliptic solve

void
MacProj (std::array< MultiFab, AMREX_SPACEDIM >& umac,
	 const MultiFab& rho,
	 const Geometry& geom,
	 const bool& full_solve)
{
    // timer for profiling
    BL_PROFILE_VAR("MacProj()",MacProj);

    BoxArray grids = umac[0].boxArray();
    // Cell-centered grids based on velocity (staggered) grids
    grids = grids.enclosedCells();
    DistributionMapping dmap = umac[0].DistributionMap();

    MultiFab solverrhs; // this will hold solver RHS = macrhs - div(beta0*umac)
    solverrhs.define(grids, dmap, 1, 0);

    MultiFab macphi(grids,dmap,1,1);
    MultiFab macrhs(grids,dmap,1,1);
    macrhs.setVal(0.0);

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

    // OR 1) average face-centered B coefficients to rho
    AverageCCToFace(rho, 0, face_bcoef, 0, 1);

    // AND 2) invert B coefficients to 1/rho
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
    mac_mlmg.setCGVerbose(cg_verbose);

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
      umac[d].FillBoundary(geom.periodicity());
    }

}

// compute the RHS for the solve, RHS = macrhs - div(beta0*umac)
void ComputeMACSolverRHS (MultiFab& solverrhs,
			  const MultiFab& macrhs,
			  const std::array< MultiFab, AMREX_SPACEDIM >& umac,
			  const Geometry& geom)
{
    // timer for profiling
    BL_PROFILE_VAR("ComputeMACSolverRHS()",ComputeMACSolverRHS);

    // Note that umac = beta0*mac
    // get references to the MultiFabs
    MultiFab& solverrhs_mf = solverrhs;
    const MultiFab& macrhs_mf = macrhs;
    const MultiFab& uedge_mf = umac[0];
#if (AMREX_SPACEDIM >= 2)
    const MultiFab& vedge_mf = umac[1];
#if (AMREX_SPACEDIM == 3)
    const MultiFab& wedge_mf = umac[2];
#endif
#endif

    int lev = 0;
    // loop over boxes
    for ( MFIter mfi(solverrhs_mf); mfi.isValid(); ++mfi) {

      // Get the index space of valid region
      const Box& validBox = mfi.validbox();
      const Real* dx = geom.CellSize();

      // call fortran subroutine
      mac_solver_rhs(&lev,ARLIM_3D(validBox.loVect()),ARLIM_3D(validBox.hiVect()),
		     BL_TO_FORTRAN_3D(solverrhs_mf[mfi]),
		     BL_TO_FORTRAN_3D(macrhs_mf[mfi]),
		     BL_TO_FORTRAN_3D(uedge_mf[mfi]),
#if (AMREX_SPACEDIM >= 2)
		     BL_TO_FORTRAN_3D(vedge_mf[mfi]),
#if (AMREX_SPACEDIM == 3)
		     BL_TO_FORTRAN_3D(wedge_mf[mfi]),
#endif
#endif
		     dx);

    }
}

// Set boundaries for MAC velocities
void SetMacSolverBCs(MLABecLaplacian& mlabec)
{
    // timer for profiling
    BL_PROFILE_VAR("SetMacSolverBCs()", SetMacSolverBCs);

    // build array of boundary conditions needed by MLABecLaplacian
    std::array<LinOpBCType, AMREX_SPACEDIM> mlmg_lobc;
    std::array<LinOpBCType, AMREX_SPACEDIM> mlmg_hibc;

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {
	if (Geometry::isPeriodic(idim)) {
            mlmg_lobc[idim] = mlmg_hibc[idim] = LinOpBCType::Periodic;
        } else {
            //amrex::Error("Invalid BC");
            mlmg_lobc[idim] = mlmg_hibc[idim] = LinOpBCType::Dirichlet;
            Print() << "Warning! non-periodic boundary conditions in MacProj_hydro." << std::endl
                    << " => Assuming Dirichlet. But preconditioner might not work properly." << std::endl;
        }
    }

    mlabec.setDomainBC(mlmg_lobc,mlmg_hibc);
}
