#include "common_functions.H"
#include "gmres_functions.H"

#include "common_functions_F.H"
#include "common_namespace.H"
#include "common_namespace_declarations.H"

#include "gmres_functions_F.H"
#include "gmres_namespace.H"
#include "gmres_namespace_declarations.H"

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
  
    // how boxes are distrubuted among MPI processes
    DistributionMapping dmap(ba);

    std::array< MultiFab, AMREX_SPACEDIM > umac_exact;
    std::array< MultiFab, AMREX_SPACEDIM > umac;
    std::array< MultiFab, AMREX_SPACEDIM > umac_tmp;
    std::array< MultiFab, AMREX_SPACEDIM > rhs_u;
    std::array< MultiFab, AMREX_SPACEDIM > grad_pres;
    std::array< MultiFab, AMREX_SPACEDIM > alpha_fc;

#if (AMREX_SPACEDIM == 2)
    std::array< MultiFab, 1 > beta_ed;
    beta_ed[0].define(convert(ba,nodal_flag), dmap, 1, 0);
#elif (AMREX_SPACEDIM == 3)
    std::array< MultiFab, 3 > beta_ed;
    beta_ed[0].define(convert(ba,nodal_flag_xy), dmap, 1, 0);
    beta_ed[1].define(convert(ba,nodal_flag_xz), dmap, 1, 0);
    beta_ed[2].define(convert(ba,nodal_flag_yz), dmap, 1, 0);
#endif

    AMREX_D_TERM(umac_exact[0].define(convert(ba,nodal_flag_x), dmap, 1, 1);,
                 umac_exact[1].define(convert(ba,nodal_flag_y), dmap, 1, 1);,
                 umac_exact[2].define(convert(ba,nodal_flag_z), dmap, 1, 1););
    AMREX_D_TERM(umac[0]      .define(convert(ba,nodal_flag_x), dmap, 1, 1);,
                 umac[1]      .define(convert(ba,nodal_flag_y), dmap, 1, 1);,
                 umac[2]      .define(convert(ba,nodal_flag_z), dmap, 1, 1););
    AMREX_D_TERM(umac_tmp[0]  .define(convert(ba,nodal_flag_x), dmap, 1, 1);,
                 umac_tmp[1]  .define(convert(ba,nodal_flag_y), dmap, 1, 1);,
                 umac_tmp[2]  .define(convert(ba,nodal_flag_z), dmap, 1, 1););
    AMREX_D_TERM(rhs_u[0]     .define(convert(ba,nodal_flag_x), dmap, 1, 0);,
                 rhs_u[1]     .define(convert(ba,nodal_flag_y), dmap, 1, 0);,
                 rhs_u[2]     .define(convert(ba,nodal_flag_z), dmap, 1, 0););
    AMREX_D_TERM(grad_pres[0] .define(convert(ba,nodal_flag_x), dmap, 1, 0);,
                 grad_pres[1] .define(convert(ba,nodal_flag_y), dmap, 1, 0);,
                 grad_pres[2] .define(convert(ba,nodal_flag_z), dmap, 1, 0););
    AMREX_D_TERM(alpha_fc[0]  .define(convert(ba,nodal_flag_x), dmap, 1, 0);,
                 alpha_fc[1]  .define(convert(ba,nodal_flag_y), dmap, 1, 0);,
                 alpha_fc[2]  .define(convert(ba,nodal_flag_z), dmap, 1, 0););

    MultiFab rhs_p     (ba, dmap, 1, 0);
    MultiFab pres_exact(ba, dmap, 1, 1);
    MultiFab pres      (ba, dmap, 1, 1);
    MultiFab pres_tmp  (ba, dmap, 1, 1);
    MultiFab alpha     (ba, dmap, 1, 1);
    MultiFab beta      (ba, dmap, 1, 1);
    MultiFab gamma     (ba, dmap, 1, 1);






    // Call the timer again and compute the maximum difference between the start time 
    // and stop time over all processors
    Real stop_time = ParallelDescriptor::second() - strt_time;
    ParallelDescriptor::ReduceRealMax(stop_time);
    amrex::Print() << "Run time = " << stop_time << std::endl;
}
