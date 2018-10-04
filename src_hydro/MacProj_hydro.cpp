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
// #include <AMReX_FMultiGrid.H>
#include <AMReX_MLABecLaplacian.H>
#include <AMReX_MLMG.H>

using namespace common;
using namespace gmres;

// ADAPTED FROM MAESTROeX
// umac enters with face-centered, time-centered Utilde^* and should leave with Utilde
// macphi is the solution to the elliptic solve and 
//   enters as either zero, or the solution to the predictor MAC projection

// FIXME: BETA IS CONSTANT
// from MAESTROeX:
// macrhs enters as beta0*(S-Sbar)
// beta0 is a 1d cell-centered array

void
MacProj (std::array< MultiFab, AMREX_SPACEDIM >& umac_in,
	   MultiFab& macphi_in,
	   const MultiFab& macrhs_in,
	   const MultiFab& rho_in,
	   const Geometry& geom_in)
{
    // timer for profiling
    BL_PROFILE_VAR("MacProj()",MacProj);

    // this will hold solver RHS = macrhs - div(beta0*umac)
    Vector<MultiFab> solverrhs(1);

    // Define single component vectors of MultiFabs to "fool" AMReX solvers
    Vector<std::array< MultiFab, AMREX_SPACEDIM > > umac(1);
    Vector<MultiFab> macphi(1);
    Vector<MultiFab> macrhs(1);
    Vector<MultiFab> rho(1);

    BoxArray ba_nodal = umac[0][0].boxArray();
    Vector<BoxArray> grids(1);
    grids[0] = ba_nodal.enclosedCells();
    Vector<DistributionMapping> dmap(1);
    dmap[0] = umac[0][0].DistributionMap();
    Vector<Geometry> geom(1);
    geom[0] = geom_in;
    
    solverrhs[0].define(grids[0], dmap[0], 1, 0);
    
    AMREX_D_TERM(umac[0][0].define(convert(grids[0],nodal_flag_x), dmap[0], 1, 1);,
		 umac[0][1].define(convert(grids[0],nodal_flag_y), dmap[0], 1, 1);,
		 umac[0][2].define(convert(grids[0],nodal_flag_z), dmap[0], 1, 1););
    rho[0].define(grids[0], dmap[0], 1, 0);
    macrhs[0].define(grids[0], dmap[0], 1, 0);
    macphi[0].define(grids[0], dmap[0], 1, 0);

    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
      MultiFab::Copy(umac[0][d],umac_in[d],0,0,1,1);
    }
    MultiFab::Copy(rho[0],rho_in,0,0,1,1);
    MultiFab::Copy(macrhs[0],macrhs_in,0,0,1,1);
    MultiFab::Copy(macphi[0],macphi_in,0,0,1,1);

    Real beta0 = 1.0;
    Real beta0_inv = 1./beta0;
    
    // convert Utilde^* to beta0*Utilde^*
    for (int d=0; d<AMREX_SPACEDIM; d++) {
      umac[0][d].mult(beta0,1);
    }

    // compute the RHS for the solve, RHS = macrhs - div(beta0*umac)
    ComputeMACSolverRHS(solverrhs[0],macrhs[0],umac[0],geom[0]);
    
    // coefficients for solver
    Vector<MultiFab> acoef(1);
    Vector<std::array< MultiFab, AMREX_SPACEDIM > > face_bcoef(1);
    acoef[0].define(grids[0], dmap[0], 1, 0);
    AMREX_D_TERM(face_bcoef[0][0].define(convert(grids[0],nodal_flag_x), dmap[0], 1, 0);,
		 face_bcoef[0][1].define(convert(grids[0],nodal_flag_y), dmap[0], 1, 0);,
		 face_bcoef[0][2].define(convert(grids[0],nodal_flag_z), dmap[0], 1, 0););

    // set cell-centered A coefficient to zero
    acoef[0].setVal(0.);

    // OR 1) average face-centered B coefficients to rho
    AverageCCToFace(rho[0], 0, face_bcoef[0], 0, 1);

    // AND 2) invert B coefficients to 1/rho
    for (int idim=0; idim<AMREX_SPACEDIM; ++idim) {
      face_bcoef[0][idim].invert(1.0,0,1);
    }

    // multiply face-centered B coefficients by beta0 so they contain beta0/rho
    for (int d=0; d<AMREX_SPACEDIM; d++) {
      umac[0][d].mult(beta0,1);
    }

    // 
    // Set up implicit solve using MLABecLaplacian class
    //
    LPInfo info;
    MLABecLaplacian mlabec(geom, grids, dmap, info);

    // order of stencil
    int linop_maxorder = 2;
    mlabec.setMaxOrder(linop_maxorder);

    // set boundaries for mlabec using velocity bc's
    SetMacSolverBCs(mlabec);

    // mlabec.setLevelBC(0, &macphi[0]);  // mlabec.setLevelBC(level, &mf)
    mlabec.setScalars(0.0, 1.0);
    mlabec.setACoeffs(0, acoef[0]);
    mlabec.setBCoeffs(0, amrex::GetArrOfConstPtrs(face_bcoef[0]));

    // solve -div B grad phi = RHS

    // build an MLMG solver
    MLMG mac_mlmg(mlabec);

    // set solver parameters
    mac_mlmg.setVerbose(mg_verbose);
    mac_mlmg.setCGVerbose(cg_verbose);

    // tolerance parameters taken from original MAESTRO fortran code
    const Real mac_tol_abs = -1.e0;
    // FIXME: using GMRES namespace parameter
    const Real mac_tol_rel = gmres_rel_tol;

    // solve for phi
    mac_mlmg.solve(GetVecOfPtrs(macphi), GetVecOfConstPtrs(solverrhs), mac_tol_rel, mac_tol_abs);

    // update velocity, beta0 * Utilde = beta0 * Utilde^* - B grad phi

    // storage for "-B grad_phi"
    std::array<MultiFab,AMREX_SPACEDIM> mac_fluxes;
    AMREX_D_TERM(mac_fluxes[0].define(convert(grids[0],nodal_flag_x), dmap[0], 1, 0);,
		 mac_fluxes[1].define(convert(grids[0],nodal_flag_y), dmap[0], 1, 0);,
		 mac_fluxes[2].define(convert(grids[0],nodal_flag_z), dmap[0], 1, 0););

    Vector< std::array<MultiFab*,AMREX_SPACEDIM> > mac_fluxptr(1);
    // fluxes computed are "-B grad phi"
    mac_fluxptr[0] = GetArrOfPtrs(mac_fluxes);
    mac_mlmg.getFluxes(mac_fluxptr);

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
      // add -B grad phi to beta0*Utilde
      MultiFab::Add(umac[0][idim], mac_fluxes[idim], 0, 0, 1, 0);
    }

    // convert beta0*Utilde to Utilde
    for (int d=0; d<AMREX_SPACEDIM; d++) {
      umac[0][d].mult(beta0_inv,1);
    }

    // fill periodic ghost cells
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
      umac[0][d].FillBoundary(geom[0].periodicity());
    }

    MultiFab::Copy(macphi_in,macphi[0],0,0,1,1);

    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
      MultiFab::Copy(umac_in[d],umac[0][d],0,0,1,1);
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

// Average bcoefs at faces using inverse of rho
void AvgFaceBcoeffsInv(std::array< MultiFab, AMREX_SPACEDIM >& facebcoef,
		       const MultiFab& rhocc)
{
    // timer for profiling
    BL_PROFILE_VAR("AvgFaceBcoeffsInv()",AvgFaceBcoeffsInv);

    // write an MFIter loop 
    // get references to the MultiFabs
    MultiFab& xbcoef_mf = facebcoef[0];
#if (AMREX_SPACEDIM >= 2)
    MultiFab& ybcoef_mf = facebcoef[1];
#if (AMREX_SPACEDIM == 3)
    MultiFab& zbcoef_mf = facebcoef[2];
#endif
#endif

    // Must get cell-centered MultiFab boxes for MIter
    const MultiFab& rhocc_mf = rhocc;

    int lev = 0;
    // loop over boxes
    for ( MFIter mfi(rhocc_mf); mfi.isValid(); ++mfi) {

      // Get the index space of valid region
      const Box& validBox = mfi.validbox();

      // call fortran subroutine
      mac_bcoef_face(&lev,ARLIM_3D(validBox.loVect()),ARLIM_3D(validBox.hiVect()), 
		     BL_TO_FORTRAN_3D(xbcoef_mf[mfi]), 
#if (AMREX_SPACEDIM >= 2)
		     BL_TO_FORTRAN_3D(ybcoef_mf[mfi]),
#if (AMREX_SPACEDIM == 3)
		     BL_TO_FORTRAN_3D(zbcoef_mf[mfi]),
#endif
#endif
		     BL_TO_FORTRAN_3D(rhocc_mf[mfi]));

    }
}

// Set boundaries for MAC velocities
void SetMacSolverBCs(MLABecLaplacian& mlabec) 
{
    // timer for profiling
    BL_PROFILE_VAR("SetMacSolverBCs()",SetMacSolverBCs);

    // build array of boundary conditions needed by MLABecLaplacian
    std::array<LinOpBCType,AMREX_SPACEDIM> mlmg_lobc;
    std::array<LinOpBCType,AMREX_SPACEDIM> mlmg_hibc;

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) 
    {
	if (Geometry::isPeriodic(idim)) {
            mlmg_lobc[idim] = mlmg_hibc[idim] = LinOpBCType::Periodic;
        }
        else {
	  amrex::Error("Invalid BC");
        }	
    }

    mlabec.setDomainBC(mlmg_lobc,mlmg_hibc);
}
