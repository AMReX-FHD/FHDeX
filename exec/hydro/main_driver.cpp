

#include "hydro_functions.H"

#include "common_functions.H"
#include "common_functions_F.H"

#include "gmres_functions.H"
#include "gmres_functions_F.H"

#include "common_namespace.H"
#include "common_namespace_declarations.H"

#include "gmres_namespace.H"
#include "gmres_namespace_declarations.H"

#include <AMReX_VisMF.H>

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

    Real dt = fixed_dt;

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

    // alpha_fc arrays
    std::array< MultiFab, AMREX_SPACEDIM > alpha_fc;
    AMREX_D_TERM(alpha_fc[0].define(convert(ba,nodal_flag_x), dmap, 1, 1);,
                 alpha_fc[1].define(convert(ba,nodal_flag_y), dmap, 1, 1);,
                 alpha_fc[2].define(convert(ba,nodal_flag_z), dmap, 1, 1););
    AMREX_D_TERM(alpha_fc[0].setVal(1.);,
                 alpha_fc[1].setVal(1.);,
                 alpha_fc[2].setVal(1.););

    // beta cell centred
    MultiFab beta(ba, dmap, 1, 1);
    beta.setVal(visc_coef*dt);

    // beta on nodes in 2d
    // beta on edges in 3d
    std::array< MultiFab, NUM_EDGE > beta_ed;
#if (AMREX_SPACEDIM == 2)
    beta_ed[0].define(convert(ba,nodal_flag), dmap, 1, 1);
    beta_ed[0].setVal(visc_coef*dt);
#elif (AMREX_SPACEDIM == 3)
    beta_ed[0].define(convert(ba,nodal_flag_xy), dmap, 1, 1);
    beta_ed[1].define(convert(ba,nodal_flag_xz), dmap, 1, 1);
    beta_ed[2].define(convert(ba,nodal_flag_yz), dmap, 1, 1);
    beta_ed[0].setVal(visc_coef*dt);  
    beta_ed[1].setVal(visc_coef*dt);
    beta_ed[2].setVal(visc_coef*dt);
#endif

    // cell-centered gamma
    MultiFab gamma(ba, dmap, 1, 1);
    gamma.setVal(0.);

    // rhs_p GMRES solve
    MultiFab gmres_rhs_p(ba, dmap, 1, 0);
    gmres_rhs_p.setVal(0.);

    // pressure for GMRES solve
    MultiFab pres(ba,dmap,1,1);
    pres.setVal(0.);  // initial guess

    // staggered velocities
    std::array< MultiFab, AMREX_SPACEDIM > umac;
    AMREX_D_TERM(umac[0].define(convert(ba,nodal_flag_x), dmap, 1, 1);,
                 umac[1].define(convert(ba,nodal_flag_y), dmap, 1, 1);,
                 umac[2].define(convert(ba,nodal_flag_z), dmap, 1, 1););

    std::array< MultiFab, AMREX_SPACEDIM > umacNew;
    AMREX_D_TERM(umacNew[0].define(convert(ba,nodal_flag_x), dmap, 1, 1);,
                 umacNew[1].define(convert(ba,nodal_flag_y), dmap, 1, 1);,
                 umacNew[2].define(convert(ba,nodal_flag_z), dmap, 1, 1););

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
    }

    // initial guess for new solution
    AMREX_D_TERM(MultiFab::Copy(umacNew[0], umac[0], 0, 0, 1, 0);,
                 MultiFab::Copy(umacNew[1], umac[1], 0, 0, 1, 0);,
                 MultiFab::Copy(umacNew[2], umac[2], 0, 0, 1, 0););

    int step = 0;
    Real time = 0.;

    Real norm_pre_rhs;

     // write out initial state
    WritePlotFile(step,time,geom,umac);

    //Time stepping loop
    for(step=1;step<=max_step;++step) {

        AMREX_D_TERM(umac[0].FillBoundary(geom.periodicity());,
                     umac[1].FillBoundary(geom.periodicity());,
                     umac[2].FillBoundary(geom.periodicity()););

        // call GMRES here
        GMRES(umac,gmres_rhs_p,umacNew,pres,alpha_fc,beta,beta_ed,gamma,1.,geom,norm_pre_rhs);

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
