
#include "hydro_test_functions.H"
#include "hydro_test_functions_F.H"

#include "hydro_functions.H"
#include "hydro_functions_F.H"

#include "common_functions.H"
#include "common_functions_F.H"

#include "gmres_functions.H"
#include "gmres_functions_F.H"

#include "common_namespace.H"
#include "common_namespace_declarations.H"

#include "gmres_namespace.H"
#include "gmres_namespace_declarations.H"

#include <AMReX_VisMF.H>
#include <AMReX_PlotFileUtil.H>

using namespace amrex;
using namespace common;
using namespace gmres;

// argv contains the name of the inputs file entered at the command line
void main_driver(const char* argv)
{

    // store the current time so we can later compute total run time.
    Real strt_time = ParallelDescriptor::second();

    std::string inputs_file = argv;

    // read in parameters from inputs file into F90 modules
    // we use "+1" because of amrex_string_c_to_f expects a null char termination
    read_common_namelist(inputs_file.c_str(),inputs_file.size()+1);
    read_gmres_namelist(inputs_file.c_str(),inputs_file.size()+1);

    // copy contents of F90 modules to C++ namespaces
    InitializeCommonNamespace();
    InitializeGmresNamespace();

    // is the problem periodic?
    Vector<int> is_periodic(AMREX_SPACEDIM,0);  // set to 0 (not periodic) by default
    for (int i=0; i<AMREX_SPACEDIM; ++i) {
        if (bc_lo[i] == -1 && bc_hi[i] == -1) {
            is_periodic[i] = 1;
        }
    }

    // make BoxArray and Geometry
    BoxArray ba;
    Geometry geom;
    {
        IntVect dom_lo(AMREX_D_DECL(           0,            0,            0));
        IntVect dom_hi(AMREX_D_DECL(n_cells[0]-1, n_cells[1]-1, n_cells[2]-1));
        Box domain(dom_lo, dom_hi);

        // Initialize the boxarray "ba" from the single box "bx"
        ba.define(domain);

        // Break up boxarray "ba" into chunks no larger than "max_grid_size" along a direction
        // note we are converting "Vector<int> max_grid_size" to an IntVect
        ba.maxSize(IntVect(max_grid_size));

       // This defines the physical box, [-1,1] in each direction.
        RealBox real_box({AMREX_D_DECL(prob_lo[0],prob_lo[1],prob_lo[2])},
                         {AMREX_D_DECL(prob_hi[0],prob_hi[1],prob_hi[2])});

        // This defines a Geometry object
        geom.define(domain,&real_box,CoordSys::cartesian,is_periodic.data());
    }

    Real dt = fixed_dt;
    Real dtinv = 1.0/dt;
    // Print() << "dt = " << dt << ", 1/dt = " << dtinv << "\n";
    const Real* dx = geom.CellSize();
  
    // how boxes are distrubuted among MPI processes
    DistributionMapping dmap(ba);

    // alpha_fc arrays
    Real theta_alpha = 1.;
    std::array< MultiFab, AMREX_SPACEDIM > alpha_fc;
    AMREX_D_TERM(alpha_fc[0].define(convert(ba,nodal_flag_x), dmap, 1, 1);,
                 alpha_fc[1].define(convert(ba,nodal_flag_y), dmap, 1, 1);,
                 alpha_fc[2].define(convert(ba,nodal_flag_z), dmap, 1, 1););
    AMREX_D_TERM(alpha_fc[0].setVal(dtinv);,
                 alpha_fc[1].setVal(dtinv);,
                 alpha_fc[2].setVal(dtinv););

    // beta cell centred
    MultiFab beta(ba, dmap, 1, 1);
    beta.setVal(visc_coef);

    // beta on nodes in 2d
    // beta on edges in 3d
    std::array< MultiFab, NUM_EDGE > beta_ed;
#if (AMREX_SPACEDIM == 2)
    beta_ed[0].define(convert(ba,nodal_flag), dmap, 1, 1);
    beta_ed[0].setVal(visc_coef);
#elif (AMREX_SPACEDIM == 3)
    beta_ed[0].define(convert(ba,nodal_flag_xy), dmap, 1, 1);
    beta_ed[1].define(convert(ba,nodal_flag_xz), dmap, 1, 1);
    beta_ed[2].define(convert(ba,nodal_flag_yz), dmap, 1, 1);
    beta_ed[0].setVal(visc_coef);  
    beta_ed[1].setVal(visc_coef);
    beta_ed[2].setVal(visc_coef);
#endif

    // cell-centered gamma
    MultiFab gamma(ba, dmap, 1, 1);
    gamma.setVal(0.);

    //////// Scaled alpha, beta, gamma: ////////
    // alpha_fc_0 arrays
    std::array< MultiFab, AMREX_SPACEDIM > alpha_fc_0;
    AMREX_D_TERM(alpha_fc_0[0].define(convert(ba,nodal_flag_x), dmap, 1, 1);,
                 alpha_fc_0[1].define(convert(ba,nodal_flag_y), dmap, 1, 1);,
                 alpha_fc_0[2].define(convert(ba,nodal_flag_z), dmap, 1, 1););
    AMREX_D_TERM(alpha_fc_0[0].setVal(0.);,
                 alpha_fc_0[1].setVal(0.);,
                 alpha_fc_0[2].setVal(0.););

    // Scaled by 1/2:
    // beta_hlf cell centered
    MultiFab beta_hlf(ba, dmap, 1, 1);
    MultiFab::Copy(beta_hlf, beta, 0, 0, 1, 1);
    beta_hlf.mult(0.5, 1);

    // beta_hlf on nodes in 2d
    // beta_hlf on edges in 3d
    std::array< MultiFab, NUM_EDGE > beta_ed_hlf;
#if (AMREX_SPACEDIM == 2)
    beta_ed_hlf[0].define(convert(ba,nodal_flag), dmap, 1, 1);
    MultiFab::Copy(beta_ed_hlf[0], beta_ed[0], 0, 0, 1, 1);
    beta_ed_hlf[0].mult(0.5, 1);
#elif (AMREX_SPACEDIM == 3)
    beta_ed_hlf[0].define(convert(ba,nodal_flag_xy), dmap, 1, 1);
    beta_ed_hlf[1].define(convert(ba,nodal_flag_xz), dmap, 1, 1);
    beta_ed_hlf[2].define(convert(ba,nodal_flag_yz), dmap, 1, 1);
    MultiFab::Copy(beta_ed_hlf[0], beta_ed[0], 0, 0, 1, 1);
    MultiFab::Copy(beta_ed_hlf[1], beta_ed[1], 0, 0, 1, 1);
    MultiFab::Copy(beta_ed_hlf[2], beta_ed[2], 0, 0, 1, 1);
    beta_ed_hlf[0].mult(0.5, 1);
    beta_ed_hlf[1].mult(0.5, 1);
    beta_ed_hlf[2].mult(0.5, 1);
#endif

    // cell-centered gamma_hlf
    MultiFab gamma_hlf(ba, dmap, 1, 1);
    MultiFab::Copy(gamma_hlf, gamma, 0, 0, 1, 1);
    gamma_hlf.mult(-0.5, 1);

    // Scaled by -1/2:
    // beta_neghlf cell centered
    MultiFab beta_neghlf(ba, dmap, 1, 1);
    MultiFab::Copy(beta_neghlf, beta, 0, 0, 1, 1);
    beta_neghlf.mult(-0.5, 1);

    // beta_neghlf on nodes in 2d
    // beta_neghlf on edges in 3d
    std::array< MultiFab, NUM_EDGE > beta_ed_neghlf;
#if (AMREX_SPACEDIM == 2)
    beta_ed_neghlf[0].define(convert(ba,nodal_flag), dmap, 1, 1);
    MultiFab::Copy(beta_ed_neghlf[0], beta_ed[0], 0, 0, 1, 1);
    beta_ed_neghlf[0].mult(-0.5, 1);
#elif (AMREX_SPACEDIM == 3)
    beta_ed_neghlf[0].define(convert(ba,nodal_flag_xy), dmap, 1, 1);
    beta_ed_neghlf[1].define(convert(ba,nodal_flag_xz), dmap, 1, 1);
    beta_ed_neghlf[2].define(convert(ba,nodal_flag_yz), dmap, 1, 1);
    MultiFab::Copy(beta_ed_neghlf[0], beta_ed[0], 0, 0, 1, 1);
    MultiFab::Copy(beta_ed_neghlf[1], beta_ed[1], 0, 0, 1, 1);
    MultiFab::Copy(beta_ed_neghlf[2], beta_ed[2], 0, 0, 1, 1);
    beta_ed_neghlf[0].mult(-0.5, 1);
    beta_ed_neghlf[1].mult(-0.5, 1);
    beta_ed_neghlf[2].mult(-0.5, 1);
#endif

    // cell-centered gamma
    MultiFab gamma_neghlf(ba, dmap, 1, 1);
    MultiFab::Copy(gamma_neghlf, gamma, 0, 0, 1, 1);
    gamma_neghlf.mult(-0.5, 1);
    ///////////////////////////////////////////

    // tracer
    MultiFab tracer(ba,dmap,1,1);
    MultiFab tracerPred(ba,dmap,1,1);
    MultiFab advFluxdivS(ba,dmap,1,1);

    // rhs_p GMRES solve
    MultiFab gmres_rhs_p(ba, dmap, 1, 0);
    gmres_rhs_p.setVal(0.);

    // pressure for GMRES solve
    MultiFab pres(ba,dmap,1,1);
    pres.setVal(0.);  // initial guess

    // rhs_u GMRES solve
    std::array< MultiFab, AMREX_SPACEDIM > gmres_rhs_u;
    AMREX_D_TERM(gmres_rhs_u[0].define(convert(ba,nodal_flag_x), dmap, 1, 1);,
                 gmres_rhs_u[1].define(convert(ba,nodal_flag_y), dmap, 1, 1);,
                 gmres_rhs_u[2].define(convert(ba,nodal_flag_z), dmap, 1, 1););
    
    // staggered velocities
    std::array< MultiFab, AMREX_SPACEDIM > umac;
    AMREX_D_TERM(umac[0].define(convert(ba,nodal_flag_x), dmap, 1, 1);,
                 umac[1].define(convert(ba,nodal_flag_y), dmap, 1, 1);,
                 umac[2].define(convert(ba,nodal_flag_z), dmap, 1, 1););

    std::array< MultiFab, AMREX_SPACEDIM > umacNew;
    AMREX_D_TERM(umacNew[0].define(convert(ba,nodal_flag_x), dmap, 1, 1);,
                 umacNew[1].define(convert(ba,nodal_flag_y), dmap, 1, 1);,
                 umacNew[2].define(convert(ba,nodal_flag_z), dmap, 1, 1););

    std::array< MultiFab, AMREX_SPACEDIM > Lumac;
    AMREX_D_TERM(Lumac[0].define(convert(ba,nodal_flag_x), dmap, 1, 1);,
                 Lumac[1].define(convert(ba,nodal_flag_y), dmap, 1, 1);,
                 Lumac[2].define(convert(ba,nodal_flag_z), dmap, 1, 1););

    // advective terms
    std::array< MultiFab, AMREX_SPACEDIM > advFluxdiv;
    AMREX_D_TERM(advFluxdiv[0].define(convert(ba,nodal_flag_x), dmap, 1, 1);,
                 advFluxdiv[1].define(convert(ba,nodal_flag_y), dmap, 1, 1);,
                 advFluxdiv[2].define(convert(ba,nodal_flag_z), dmap, 1, 1););

    std::array< MultiFab, AMREX_SPACEDIM > advFluxdivPred;
    AMREX_D_TERM(advFluxdivPred[0].define(convert(ba,nodal_flag_x), dmap, 1, 1);,
                 advFluxdivPred[1].define(convert(ba,nodal_flag_y), dmap, 1, 1);,
                 advFluxdivPred[2].define(convert(ba,nodal_flag_z), dmap, 1, 1););

    // staggered momentum
    std::array< MultiFab, AMREX_SPACEDIM > uMom;
    AMREX_D_TERM(uMom[0].define(convert(ba,nodal_flag_x), dmap, 1, 1);,
                 uMom[1].define(convert(ba,nodal_flag_y), dmap, 1, 1);,
                 uMom[2].define(convert(ba,nodal_flag_z), dmap, 1, 1););

    const RealBox& realDomain = geom.ProbDomain();
    int dm;

    for ( MFIter mfi(beta); mfi.isValid(); ++mfi ) {
        const Box& bx = mfi.validbox();
        
        AMREX_D_TERM(dm=0; init_vel(BL_TO_FORTRAN_BOX(bx),
                                    BL_TO_FORTRAN_ANYD(umac[0][mfi]), geom.CellSize(),
                                    geom.ProbLo(), geom.ProbHi() ,&dm, 
                                    ZFILL(realDomain.lo()), ZFILL(realDomain.hi()));,
                     dm=1; init_vel(BL_TO_FORTRAN_BOX(bx),
                                    BL_TO_FORTRAN_ANYD(umac[1][mfi]), geom.CellSize(),
                                    geom.ProbLo(), geom.ProbHi() ,&dm, 
                                    ZFILL(realDomain.lo()), ZFILL(realDomain.hi()));,
                     dm=2; init_vel(BL_TO_FORTRAN_BOX(bx),
                                    BL_TO_FORTRAN_ANYD(umac[2][mfi]), geom.CellSize(),
                                    geom.ProbLo(), geom.ProbHi() ,&dm, 
                                    ZFILL(realDomain.lo()), ZFILL(realDomain.hi())););

	// initialize tracer
        // init_s_vel(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()), 
	// 	   BL_TO_FORTRAN_ANYD(tracer[mfi]), geom.CellSize(), 
	// 	   ZFILL(realDomain.lo()), ZFILL(realDomain.hi()));

        init_s_vel(BL_TO_FORTRAN_BOX(bx),
		   BL_TO_FORTRAN_ANYD(tracer[mfi]),
		   dx, ZFILL(realDomain.lo()), ZFILL(realDomain.hi()));

        // init_s_vel(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()), 
	// 	   BL_TO_FORTRAN_ANYD(tracer[mfi]), dx);
    }

    // initial guess for new solution
    AMREX_D_TERM(MultiFab::Copy(umacNew[0], umac[0], 0, 0, 1, 0);,
                 MultiFab::Copy(umacNew[1], umac[1], 0, 0, 1, 0);,
                 MultiFab::Copy(umacNew[2], umac[2], 0, 0, 1, 0););

    int step = 0;
    Real time = 0.;

    Real norm_pre_rhs;

    // write out initial state
    WritePlotFile(step,time,geom,umac,tracer,pres);

    //Time stepping loop
    for(step=1;step<=max_step;++step) {

        AMREX_D_TERM(umac[0].FillBoundary(geom.periodicity());,
                     umac[1].FillBoundary(geom.periodicity());,
                     umac[2].FillBoundary(geom.periodicity()););

	// Compute tracer:
	if (step != 1) {
	  tracer.FillBoundary(geom.periodicity());
	  MkAdvSFluxdiv(umac,tracer,advFluxdivS,dx,geom,0);
	  advFluxdivS.mult(dt, 1);

	  // compute predictor
	  MultiFab::Copy(tracerPred, tracer, 0, 0, 1, 0);
	  MultiFab::Add(tracerPred, advFluxdivS, 0, 0, 1, 0);
	  tracerPred.FillBoundary(geom.periodicity());
	  MkAdvSFluxdiv(umac,tracerPred,advFluxdivS,dx,geom,0);
	  advFluxdivS.mult(dt, 1);

	  // advance in time
	  MultiFab::Add(tracer, tracerPred, 0, 0, 1, 0);
	  MultiFab::Add(tracer, advFluxdivS, 0, 0, 1, 0);
	  tracer.mult(0.5, 1);

	  // amrex::Print() << "tracer L0 norm = " << tracer.norm0() << "\n";
	}

	// PREDICTOR STEP (trapezoidal rule)
	// compute advective term
        AMREX_D_TERM(MultiFab::Copy(uMom[0], umac[0], 0, 0, 1, 0);,
                     MultiFab::Copy(uMom[1], umac[1], 0, 0, 1, 0);,
                     MultiFab::Copy(uMom[2], umac[2], 0, 0, 1, 0););

	// let rho = 1
	for (int d=0; d<AMREX_SPACEDIM; d++) {
	  uMom[d].mult(1.0, 1);
	}

        AMREX_D_TERM(uMom[0].FillBoundary(geom.periodicity());,
                     uMom[1].FillBoundary(geom.periodicity());,
                     uMom[2].FillBoundary(geom.periodicity()););

	MkAdvMFluxdiv(umac,uMom,advFluxdiv,dx,0);

	AMREX_D_TERM(MultiFab::Copy(gmres_rhs_u[0], umac[0], 0, 0, 1, 0);,
                     MultiFab::Copy(gmres_rhs_u[1], umac[1], 0, 0, 1, 0);,
                     MultiFab::Copy(gmres_rhs_u[2], umac[2], 0, 0, 1, 0););
	for (int d=0; d<AMREX_SPACEDIM; d++) {
	  gmres_rhs_u[d].mult(dtinv, 1);
	}
	AMREX_D_TERM(MultiFab::Add(gmres_rhs_u[0], advFluxdiv[0], 0, 0, 1, 0);,
                     MultiFab::Add(gmres_rhs_u[1], advFluxdiv[1], 0, 0, 1, 0);,
                     MultiFab::Add(gmres_rhs_u[2], advFluxdiv[2], 0, 0, 1, 0););

        AMREX_D_TERM(gmres_rhs_u[0].FillBoundary(geom.periodicity());,
                     gmres_rhs_u[1].FillBoundary(geom.periodicity());,
                     gmres_rhs_u[2].FillBoundary(geom.periodicity()););

	// initial guess for new solution
	AMREX_D_TERM(MultiFab::Copy(umacNew[0], umac[0], 0, 0, 1, 0);,
		     MultiFab::Copy(umacNew[1], umac[1], 0, 0, 1, 0);,
		     MultiFab::Copy(umacNew[2], umac[2], 0, 0, 1, 0););
	pres.setVal(0.);  // initial guess

        // call GMRES to compute predictor
	GMRES(gmres_rhs_u,gmres_rhs_p,umacNew,pres,alpha_fc,beta,beta_ed,gamma,theta_alpha,geom,norm_pre_rhs);

	// Compute predictor advective term
        AMREX_D_TERM(umacNew[0].FillBoundary(geom.periodicity());,
                     umacNew[1].FillBoundary(geom.periodicity());,
                     umacNew[2].FillBoundary(geom.periodicity()););

        AMREX_D_TERM(MultiFab::Copy(uMom[0], umacNew[0], 0, 0, 1, 0);,
                     MultiFab::Copy(uMom[1], umacNew[1], 0, 0, 1, 0);,
                     MultiFab::Copy(uMom[2], umacNew[2], 0, 0, 1, 0););

	// let rho = 1
	for (int d=0; d<AMREX_SPACEDIM; d++) {
	  uMom[d].mult(1.0, 1);
	}

        AMREX_D_TERM(uMom[0].FillBoundary(geom.periodicity());,
                     uMom[1].FillBoundary(geom.periodicity());,
                     uMom[2].FillBoundary(geom.periodicity()););

	MkAdvMFluxdiv(umacNew,uMom,advFluxdivPred,dx,0);

	// ADVANCE STEP (crank nicolson + trapezoidal rule)

	/////////////// Hack /////////////////////////////
	// VisMF::Write(advFluxdiv[0],"a_advFluxdiv0");

	// StagApplyOp(beta,gamma,beta_ed,umac,advFluxdiv,alpha_fc_0,dx,theta_alpha);

	// VisMF::Write(advFluxdiv[0],"a_Lumac0");
	// Abort("Done with hack");
	// exit(0);
	//////////////////////////////////////////////////

	// Compute gmres_rhs

        // trapezoidal advective terms
	for (int d=0; d<AMREX_SPACEDIM; d++) {
	  advFluxdiv[d].mult(0.5, 1);
	  advFluxdivPred[d].mult(0.5, 1);
	}

        // crank-nicolson terms
        StagApplyOp(beta_neghlf,gamma_neghlf,beta_ed_neghlf,umac,Lumac,alpha_fc_0,dx,theta_alpha);

	AMREX_D_TERM(MultiFab::Copy(gmres_rhs_u[0], umac[0], 0, 0, 1, 0);,
                     MultiFab::Copy(gmres_rhs_u[1], umac[1], 0, 0, 1, 0);,
                     MultiFab::Copy(gmres_rhs_u[2], umac[2], 0, 0, 1, 0););
	for (int d=0; d<AMREX_SPACEDIM; d++) {
	  gmres_rhs_u[d].mult(dtinv, 1);
	}
	AMREX_D_TERM(MultiFab::Add(gmres_rhs_u[0], Lumac[0], 0, 0, 1, 0);,
                     MultiFab::Add(gmres_rhs_u[1], Lumac[1], 0, 0, 1, 0);,
                     MultiFab::Add(gmres_rhs_u[2], Lumac[2], 0, 0, 1, 0););
	AMREX_D_TERM(MultiFab::Add(gmres_rhs_u[0], advFluxdiv[0], 0, 0, 1, 0);,
                     MultiFab::Add(gmres_rhs_u[1], advFluxdiv[1], 0, 0, 1, 0);,
                     MultiFab::Add(gmres_rhs_u[2], advFluxdiv[2], 0, 0, 1, 0););
	AMREX_D_TERM(MultiFab::Add(gmres_rhs_u[0], advFluxdivPred[0], 0, 0, 1, 0);,
                     MultiFab::Add(gmres_rhs_u[1], advFluxdivPred[1], 0, 0, 1, 0);,
                     MultiFab::Add(gmres_rhs_u[2], advFluxdivPred[2], 0, 0, 1, 0););

        AMREX_D_TERM(gmres_rhs_u[0].FillBoundary(geom.periodicity());,
                     gmres_rhs_u[1].FillBoundary(geom.periodicity());,
                     gmres_rhs_u[2].FillBoundary(geom.periodicity()););

	// initial guess for new solution
	AMREX_D_TERM(MultiFab::Copy(umacNew[0], umac[0], 0, 0, 1, 0);,
		     MultiFab::Copy(umacNew[1], umac[1], 0, 0, 1, 0);,
		     MultiFab::Copy(umacNew[2], umac[2], 0, 0, 1, 0););
	pres.setVal(0.);  // initial guess

        // call GMRES here
	GMRES(gmres_rhs_u,gmres_rhs_p,umacNew,pres,alpha_fc,beta_hlf,beta_ed_hlf,gamma_hlf,theta_alpha,geom,norm_pre_rhs);
	// GMRES(gmres_rhs_u,gmres_rhs_p,umacNew,pres,alpha_fc,beta,beta_ed,gamma,theta_alpha,geom,norm_pre_rhs);

        AMREX_D_TERM(MultiFab::Copy(umac[0], umacNew[0], 0, 0, 1, 0);,
                     MultiFab::Copy(umac[1], umacNew[1], 0, 0, 1, 0);,
                     MultiFab::Copy(umac[2], umacNew[2], 0, 0, 1, 0););

        amrex::Print() << "Advanced step " << step << "\n";

        time = time + dt;

        if (plot_int > 0 && step%plot_int == 0) {
          // write out umac & pres to a plotfile
	  WritePlotFile(step,time,geom,umac,tracer,pres);
        }
    }

    // Call the timer again and compute the maximum difference between the start time 
    // and stop time over all processors
    Real stop_time = ParallelDescriptor::second() - strt_time;
    ParallelDescriptor::ReduceRealMax(stop_time);
    amrex::Print() << "Run time = " << stop_time << std::endl;
    amrex::Print() << "dt = " << dt << "\n";

}
