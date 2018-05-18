

#include "common_functions.H"
#include "gmres_functions.H"

#include "common_functions_F.H"
#include "common_namespace.H"
#include "common_namespace_declarations.H"

#include "gmres_functions_F.H"
#include "gmres_namespace.H"
#include "gmres_namespace_declarations.H"

#include "INS_functions.H"

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
        if (bc_lo[i] == 0 && bc_hi[i] == 0) {
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

    // total density
    MultiFab rhotot(ba, dmap, 1, 1);

    // divergence
    MultiFab div(ba, dmap, 1, 1);

    // potential
    MultiFab phi(ba, dmap, 1, 1);

    // beta cell centred
    MultiFab betaCC(ba, dmap, 1, 1);

    // gamma cell centred
    MultiFab gammaCC(ba, dmap, 1, 1);

    //Print() << nodal_flag_xy << "\n";
    //while(true);

    // beta and gamma on nodes in 2d
#if (AMREX_SPACEDIM == 2)
    MultiFab betaNodal(convert(ba,nodal_flag_xy), dmap, 1, 1);
    MultiFab gammaNodal(convert(ba,nodal_flag_xy), dmap, 1, 1);

    betaNodal.setVal(1.);
    gammaNodal.setVal(0);
#endif

    // beta and gamma on edges in 3d
#if (AMREX_SPACEDIM == 3)
    std::array< MultiFab, AMREX_SPACEDIM > betaEdge;
    AMREX_D_TERM(betaEdge[0].define(convert(ba,nodal_flag_xy), dmap, 1, 1);,
                 betaEdge[1].define(convert(ba,nodal_flag_xz), dmap, 1, 1);,
                 betaEdge[2].define(convert(ba,nodal_flag_yz), dmap, 1, 1););

    std::array< MultiFab, AMREX_SPACEDIM > gammaEdge;
    AMREX_D_TERM(gammaEdge[0].define(convert(ba,nodal_flag_xy), dmap, 1, 1);,
                 gammaEdge[1].define(convert(ba,nodal_flag_xz), dmap, 1, 1);,
                 gammaEdge[2].define(convert(ba,nodal_flag_yz), dmap, 1, 1););

    betaEdge[0].setVal(1.);  
    betaEdge[1].setVal(1.);
    betaEdge[2].setVal(1.);

    gammaEdge[0].setVal(0.);    
    gammaEdge[1].setVal(0.);
    gammaEdge[2].setVal(0.);

#endif

    //Replace with proper initialiser
    phi.setVal(100.);

    betaCC.setVal(1.);
    gammaCC.setVal(0);

    // temporary placeholder for potential gradients on cell faces
    std::array< MultiFab, AMREX_SPACEDIM > umacT;
    AMREX_D_TERM(umacT[0].define(convert(ba,nodal_flag_x), dmap, 1, 1);,
                 umacT[1].define(convert(ba,nodal_flag_y), dmap, 1, 1);,
                 umacT[2].define(convert(ba,nodal_flag_z), dmap, 1, 1););

    // set density to 1
    rhotot.setVal(1.);

    // staggered velocities
    std::array< MultiFab, AMREX_SPACEDIM > umac;
    AMREX_D_TERM(umac[0].define(convert(ba,nodal_flag_x), dmap, 1, 1);,
                 umac[1].define(convert(ba,nodal_flag_y), dmap, 1, 1);,
                 umac[2].define(convert(ba,nodal_flag_z), dmap, 1, 1););

    // alpha arrays
    std::array< MultiFab, AMREX_SPACEDIM > alpha;
    AMREX_D_TERM(alpha[0].define(convert(ba,nodal_flag_x), dmap, 1, 1);
                 alpha[0].setVal(0);,
                 alpha[1].define(convert(ba,nodal_flag_y), dmap, 1, 1);
                 alpha[1].setVal(0);,
                 alpha[2].define(convert(ba,nodal_flag_z), dmap, 1, 1);
                 alpha[2].setVal(0););

    // For testing timestepping
    std::array< MultiFab, AMREX_SPACEDIM > umacNew;
    AMREX_D_TERM(umacNew[0].define(convert(ba,nodal_flag_x), dmap, 1, 1);,
                 umacNew[1].define(convert(ba,nodal_flag_y), dmap, 1, 1);,
                 umacNew[2].define(convert(ba,nodal_flag_z), dmap, 1, 1););

    // For testing StagApplyOp
    std::array< MultiFab, AMREX_SPACEDIM > umacOut;
    AMREX_D_TERM(umacOut[0].define(convert(ba,nodal_flag_x), dmap, 1, 1);,
                 umacOut[1].define(convert(ba,nodal_flag_y), dmap, 1, 1);,
                 umacOut[2].define(convert(ba,nodal_flag_z), dmap, 1, 1););


    // ***REPLACE THIS WITH A FUNCTION THAT SETS THE INITIAL VELOCITY***
    // ***SETTING THESE TO DUMMY VALUES FOR NOW***
    //AMREX_D_TERM(umac[0].setVal(100.);,
    //             umac[1].setVal(100.);,
     //            umac[2].setVal(100.););

	int dm = 0;
	for ( MFIter mfi(rhotot); mfi.isValid(); ++mfi )
    {
        const Box& bx = mfi.validbox();

		AMREX_D_TERM(dm=0; init_vel(BL_TO_FORTRAN_BOX(bx),
							BL_TO_FORTRAN_ANYD(umac[0][mfi]), geom.CellSize(),
            				geom.ProbLo(), geom.ProbHi() ,&dm);,
					dm=1; init_vel(BL_TO_FORTRAN_BOX(bx),
            			    BL_TO_FORTRAN_ANYD(umac[1][mfi]), geom.CellSize(),
            			    geom.ProbLo(), geom.ProbHi() ,&dm);,
					dm=2; init_vel(BL_TO_FORTRAN_BOX(bx),
                  			BL_TO_FORTRAN_ANYD(umac[2][mfi]), geom.CellSize(),
                   			geom.ProbLo(), geom.ProbHi() ,&dm););
    }

    // compute the time step
    const Real* dx = geom.CellSize();
    Real dt = 0.9*dx[0]*dx[0] / (2.0*AMREX_SPACEDIM);

    int step = 0;
    Real time = 0.;

    //Particles!

    //ParticleContainer<2*AMREX_SPACEDIM,0> particles(geom, dmap, ba);

    // write out initial state
    WritePlotFile(step,time,geom,rhotot,umac,div);


    //Time stepping loop
    for(step=1;step<=max_step;++step)
    {
        AMREX_D_TERM(
        umac[0].FillBoundary(geom.periodicity());,
        umac[1].FillBoundary(geom.periodicity());,
        umac[2].FillBoundary(geom.periodicity()););

        eulerStep(betaCC, gammaCC, 
#if (AMREX_SPACEDIM == 2)
        betaNodal, gammaNodal,
#endif
#if (AMREX_SPACEDIM == 3)
        betaEdge, gammaEdge,
#endif
        umac, umacOut, umacNew, alpha, geom, visc_type, &dt);


        AMREX_D_TERM(
        MultiFab::Copy(umac[0], umacNew[0], 0, 0, 1, 0);,
        MultiFab::Copy(umac[1], umacNew[1], 0, 0, 1, 0);,
        MultiFab::Copy(umac[2], umacNew[2], 0, 0, 1, 0););


        amrex::Print() << "Advanced step " << step << "\n";

        time = time + dt;

        if (plot_int > 0 && step%plot_int == 0)
        {
            // write out rhotot and umac to a plotfile
            WritePlotFile(step,time,geom,rhotot,umac,div);
        }


    }


    //Compute divergence, last arguement is flag for increment
    
    //ComputeDiv(div,umac,geom,0);

    //Compute gradient
    //ComputeGrad(phi,umacT,geom);


    // Call the timer again and compute the maximum difference between the start time 
    // and stop time over all processors
    Real stop_time = ParallelDescriptor::second() - strt_time;
    ParallelDescriptor::ReduceRealMax(stop_time);
    amrex::Print() << "Run time = " << stop_time << std::endl;
}
