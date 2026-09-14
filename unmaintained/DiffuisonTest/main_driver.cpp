#include "DiffusionTest_functions.H"

#include "common_functions.H"

#include "gmres_functions.H"



#include <AMReX_VisMF.H>

using namespace amrex;

// argv contains the name of the inputs file entered at the command line
void main_driver(const char* argv)
{

    // store the current time so we can later compute total run time.
    Real strt_time = ParallelDescriptor::second();

    std::string inputs_file = argv;

    // read in parameters from inputs file into F90 modules
    // we use "+1" because of amrex_string_c_to_f expects a null char termination
    read_common_namelist(inputs_file.c_str(),inputs_file.size()+1);

    // copy contents of F90 modules to C++ namespaces
    InitializeCommonNamespace();
    InitializeGmresNamespace();

    // is the problem periodic?
    Vector<int> is_periodic(AMREX_SPACEDIM,0);  // set to 0 (not periodic) by default
    for (int i=0; i<AMREX_SPACEDIM; ++i) {
        if (bc_vel_lo[i] == -1 && bc_vel_hi[i] == -1) {
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

    // how boxes are distrubuted among MPI processes
    DistributionMapping dmap(ba);

    // beta cell centred
    MultiFab betaCC(ba, dmap, 1, 1);

    // gamma cell centred
    MultiFab gammaCC(ba, dmap, 1, 1);

    Real dt = fixed_dt;
    const Real* dx = geom.CellSize();
    if (fixed_dt < 0.) {
        dt = 0.5*(dx[0]*dx[0])/(2.0*AMREX_SPACEDIM*1.0);
    }

    // beta on nodes in 2d
    // beta on edges in 3d
    std::array< MultiFab, NUM_EDGE > betaEdge;
#if (AMREX_SPACEDIM == 2)
    betaEdge[0].define(convert(ba,nodal_flag), dmap, 1, 1);

    if ( algorithm_type == 1 || algorithm_type == 2) {
        betaEdge[0].setVal(visc_coef*dt);
    } else if (algorithm_type == 0) {
        betaEdge[0].setVal(-visc_coef*dt);
    } else {
        Print() << "Error: Invalid choice of algorithm_type\n";
    }

#elif (AMREX_SPACEDIM == 3)
    betaEdge[0].define(convert(ba,nodal_flag_xy), dmap, 1, 1);
    betaEdge[1].define(convert(ba,nodal_flag_xz), dmap, 1, 1);
    betaEdge[2].define(convert(ba,nodal_flag_yz), dmap, 1, 1);

    if ( algorithm_type == 1 || algorithm_type == 2) {
        betaEdge[0].setVal(visc_coef*dt);
        betaEdge[1].setVal(visc_coef*dt);
        betaEdge[2].setVal(visc_coef*dt);
    } else if (algorithm_type == 0) {
        betaEdge[0].setVal(-visc_coef*dt);
        betaEdge[1].setVal(-visc_coef*dt);
        betaEdge[2].setVal(-visc_coef*dt);
    } else {
        Print() << "Error: Invalid choice of algorithm_type\n";
    }

#endif

    // staggered velocities
    std::array< MultiFab, AMREX_SPACEDIM > umac;
    AMREX_D_TERM(umac[0].define(convert(ba,nodal_flag_x), dmap, 1, 1);,
                 umac[1].define(convert(ba,nodal_flag_y), dmap, 1, 1);,
                 umac[2].define(convert(ba,nodal_flag_z), dmap, 1, 1););

    // alpha arrays
    std::array< MultiFab, AMREX_SPACEDIM > alpha;
    AMREX_D_TERM(alpha[0].define(convert(ba,nodal_flag_x), dmap, 1, 1);,
                 alpha[1].define(convert(ba,nodal_flag_y), dmap, 1, 1);,
                 alpha[2].define(convert(ba,nodal_flag_z), dmap, 1, 1););


    if ( algorithm_type == 1 || algorithm_type == 2) {
        betaCC.setVal(visc_coef*dt);
        gammaCC.setVal(0.);
        AMREX_D_TERM(alpha[0].setVal(1.);,
                     alpha[1].setVal(1.);,
                     alpha[2].setVal(1.););
    } else if (algorithm_type == 0) {
        betaCC.setVal(-visc_coef*dt);
        gammaCC.setVal(0.);
        AMREX_D_TERM(alpha[0].setVal(0.);,
                     alpha[1].setVal(0.);,
                     alpha[2].setVal(0.););
    } else {
        Print() << "Error: Invalid choice of algorithm_type\n";
    }

    // For testing timestepping
    std::array< MultiFab, AMREX_SPACEDIM > umacNew;
    AMREX_D_TERM(umacNew[0].define(convert(ba,nodal_flag_x), dmap, 1, 1);,
                 umacNew[1].define(convert(ba,nodal_flag_y), dmap, 1, 1);,
                 umacNew[2].define(convert(ba,nodal_flag_z), dmap, 1, 1););

    const RealBox& realDomain = geom.ProbDomain();

    int dm = 0;
    for ( MFIter mfi(betaCC); mfi.isValid(); ++mfi ) {
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
    }

    int step = 0;
    Real time = 0.;

    // write out initial state
    WritePlotFile(step,time,geom,umac);

    // initial guess for updated velocity
    AMREX_D_TERM(MultiFab::Copy(umacNew[0], umac[0], 0, 0, 1, 0);,
                 MultiFab::Copy(umacNew[1], umac[1], 0, 0, 1, 0);,
                 MultiFab::Copy(umacNew[2], umac[2], 0, 0, 1, 0););

    // definitions for Crank-Nicolson:
    // right hand side of poisson eqn std::array< MultiFab, AMREX_SPACEDIM > rhs; AMREX_D_TERM(rhs[0].define(convert(ba,nodal_flag_x), dmap, 1, 1);, rhs[1].define(convert(ba,nodal_flag_y), dmap, 1, 1);, rhs[2].define(convert(ba,nodal_flag_z), dmap, 1, 1);); MultiFab betaCC_hlf(ba, dmap, 1, 1); MultiFab gammaCC_hlf(ba, dmap, 1, 1);
    std::array< MultiFab, NUM_EDGE > betaEdge_hlf;
#if (AMREX_SPACEDIM == 2)
    betaEdge_hlf[0].define(convert(ba,nodal_flag), dmap, 1, 1);
#elif (AMREX_SPACEDIM == 3)
    betaEdge_hlf[0].define(convert(ba,nodal_flag_xy), dmap, 1, 1);
    betaEdge_hlf[1].define(convert(ba,nodal_flag_xz), dmap, 1, 1);
    betaEdge_hlf[2].define(convert(ba,nodal_flag_yz), dmap, 1, 1);
#endif

    MultiFab betaCC_neghlf(ba, dmap, 1, 1);
    MultiFab gammaCC_neghlf(ba, dmap, 1, 1);
    std::array< MultiFab, NUM_EDGE > betaEdge_neghlf;
#if (AMREX_SPACEDIM == 2)
    betaEdge_neghlf[0].define(convert(ba,nodal_flag), dmap, 1, 1);
#elif (AMREX_SPACEDIM == 3)
    betaEdge_neghlf[0].define(convert(ba,nodal_flag_xy), dmap, 1, 1);
    betaEdge_neghlf[1].define(convert(ba,nodal_flag_xz), dmap, 1, 1);
    betaEdge_neghlf[2].define(convert(ba,nodal_flag_yz), dmap, 1, 1);
#endif

    if (algorithm_type == 2) {

        AMREX_D_TERM(umac[0].FillBoundary(geom.periodicity());,
                     umac[1].FillBoundary(geom.periodicity());,
                     umac[2].FillBoundary(geom.periodicity()););

        Real weight_lap;

        // multiply beta & gamma by weight_lap = -1/2 : to compute rhs
        weight_lap = -0.5;
        MultiFab::Copy( betaCC_neghlf, betaCC,0,0,1,1);
        MultiFab::Copy(gammaCC_neghlf,gammaCC,0,0,1,1);
        betaCC_neghlf.mult(weight_lap,0,1,1);
        gammaCC_neghlf.mult(weight_lap,0,1,1);
#if (AMREX_SPACEDIM == 2)
        MultiFab::Copy( betaEdge_neghlf[0],betaEdge[0],0,0,1,1);
        betaEdge_neghlf[0].mult(weight_lap,0,1,1);
#elif (AMREX_SPACEDIM == 3)
        MultiFab::Copy(betaEdge_neghlf[0],betaEdge[0],0,0,1,1);
        MultiFab::Copy(betaEdge_neghlf[1],betaEdge[1],0,0,1,1);
        MultiFab::Copy(betaEdge_neghlf[2],betaEdge[2],0,0,1,1);
        betaEdge_neghlf[0].mult(weight_lap,0,1,1);
        betaEdge_neghlf[1].mult(weight_lap,0,1,1);
        betaEdge_neghlf[2].mult(weight_lap,0,1,1);
#endif
        // Fill ghost cells for good measure
        betaCC_neghlf.FillBoundary(geom.periodicity());
        gammaCC_neghlf.FillBoundary(geom.periodicity());
        AMREX_D_TERM(betaEdge_neghlf[0].FillBoundary(geom.periodicity());,
                     betaEdge_neghlf[1].FillBoundary(geom.periodicity());,
                     betaEdge_neghlf[2].FillBoundary(geom.periodicity()););

        // multiply beta & gamma by weight_lap = 1/2 : to correct L-matrix
        weight_lap = 0.5;
        MultiFab::Copy( betaCC_hlf, betaCC,0,0,1,1);
        MultiFab::Copy(gammaCC_hlf,gammaCC,0,0,1,1);
        betaCC_hlf.mult(weight_lap,0,1,1);
        gammaCC_hlf.mult(weight_lap,0,1,1);
#if (AMREX_SPACEDIM == 2)
        MultiFab::Copy( betaEdge_hlf[0],betaEdge[0],0,0,1,1);
        betaEdge_hlf[0].mult(weight_lap,0,1,1);
#elif (AMREX_SPACEDIM == 3)
        MultiFab::Copy(betaEdge_hlf[0],betaEdge[0],0,0,1,1);
        MultiFab::Copy(betaEdge_hlf[1],betaEdge[1],0,0,1,1);
        MultiFab::Copy(betaEdge_hlf[2],betaEdge[2],0,0,1,1);
        betaEdge_hlf[0].mult(weight_lap,0,1,1);
        betaEdge_hlf[1].mult(weight_lap,0,1,1);
        betaEdge_hlf[2].mult(weight_lap,0,1,1);
#endif
        // Fill ghost cells for good measure
        betaCC_hlf.FillBoundary(geom.periodicity());
        gammaCC_hlf.FillBoundary(geom.periodicity());
        AMREX_D_TERM(betaEdge_hlf[0].FillBoundary(geom.periodicity());,
                     betaEdge_hlf[1].FillBoundary(geom.periodicity());,
                     betaEdge_hlf[2].FillBoundary(geom.periodicity()););

    }

    //Time stepping loop
    for(step=1;step<=max_step;++step) {

        AMREX_D_TERM(umac[0].FillBoundary(geom.periodicity());,
                     umac[1].FillBoundary(geom.periodicity());,
                     umac[2].FillBoundary(geom.periodicity()););

        if (algorithm_type == 2) {
            // compute rhs of poisson eqn
            StagApplyOp(geom,betaCC_neghlf,gammaCC_neghlf,betaEdge_neghlf,umac,rhs,alpha,dx,1.);

            AMREX_D_TERM(rhs[0].FillBoundary(geom.periodicity());,
                         rhs[1].FillBoundary(geom.periodicity());,
                         rhs[2].FillBoundary(geom.periodicity()););

            // VisMF::Write(rhs[0],"a_rhs0");
            // exit(0);

            StagMGSolver(alpha,betaCC_hlf,betaEdge_hlf,gammaCC_hlf,umacNew,rhs,1.0,geom);
        } else if (algorithm_type == 1) {
            StagMGSolver(alpha,betaCC,betaEdge,gammaCC,umacNew,umac,1.0,geom);
        } else if (algorithm_type == 0) {
            StagExpSolver(alpha,betaCC,betaEdge,gammaCC,umacNew,umac,1.0,geom);
        } else {
            Print() << "Error: Invalid choice of algorithm_type\n";
        }

        AMREX_D_TERM(MultiFab::Copy(umac[0], umacNew[0], 0, 0, 1, 0);,
                     MultiFab::Copy(umac[1], umacNew[1], 0, 0, 1, 0);,
                     MultiFab::Copy(umac[2], umacNew[2], 0, 0, 1, 0););

        amrex::Print() << "Advanced step " << step << "\n";

        time = time + dt;

        if (plot_int > 0 && step%plot_int == 0) {
            // write out umac to a plotfile
            WritePlotFile(step,time,geom,umac);
        }
    }

    // Call the timer again and compute the maximum difference between the start time
    // and stop time over all processors
    Real stop_time = ParallelDescriptor::second() - strt_time;
    ParallelDescriptor::ReduceRealMax(stop_time);
    amrex::Print() << "Run time = " << stop_time << std::endl;
    amrex::Print() << "dt = " << dt << "\n";

}