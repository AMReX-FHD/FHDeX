#include <AMReX_Gpu.H>
#include <AMReX_Utility.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_BCRec.H>
// #include <AMReX_BCUtil.H>

#include "common_functions.H"
#include "rng_functions.H"

#include "myfunc.H"
#include "chrono"

using namespace amrex;
using namespace std::chrono;

int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);

    main_main();

    amrex::Finalize();
    return 0;
}

void main_main ()
{
    // What time is it now?  We'll use this to compute total run time.
    auto strt_time = ParallelDescriptor::second();

    // AMREX_SPACEDIM: number of dimensions
    int n_cell, max_grid_size, nsteps, plot_int;
    int alg_type;
    amrex::Real npts_scale;
    amrex::Real cfl;
    Vector<int> bc_lo(AMREX_SPACEDIM,0);
    Vector<int> bc_hi(AMREX_SPACEDIM,0);

    // inputs parameters
    {
        // ParmParse is way of reading inputs from the inputs file
        ParmParse pp;

        // We need to get n_cell from the inputs file - this is the number of cells on each side of
        //   a square (or cubic) domain.
        pp.get("n_cell",n_cell);

        // The domain is broken into boxes of size max_grid_size
        pp.get("max_grid_size",max_grid_size);

        // Default plot_int to -1, allow us to set it to something else in the inputs file
        //  If plot_int < 0 then no plot files will be written
        plot_int = -1;
        pp.query("plot_int",plot_int);

        // Default nsteps to 0, allow us to set it to something else in the inputs file
        nsteps = 10;
        pp.query("nsteps",nsteps);

        npts_scale = 1.;
        pp.query ("npts_scale",npts_scale);

        alg_type = 0;
        pp.query ("alg_type",alg_type);

        cfl=.9;
        pp.query ("cfl",cfl);


        // By default, the boundary conditions will be set to periodic, or bc_lo = bc_hi = 0.
        //Other options in this program include bc_lo, bc_hi = 2 for homogeneous Neumann, or
        //bc_lo, bc_hi = 3 for external Dirichlet boundary conditions.
        pp.queryarr("bc_lo", bc_lo);
        pp.queryarr("bc_hi", bc_hi);
    }

    Vector<int> is_periodic(AMREX_SPACEDIM,0);
    for (int idim=0; idim < AMREX_SPACEDIM; ++idim) {
        if (bc_lo[idim] == BCType::int_dir && bc_hi[idim] == BCType::int_dir){
            is_periodic[idim] = 1;
        }
    }

    // make BoxArray and Geometry
    BoxArray ba;
    Geometry geom;
    {
        IntVect dom_lo(AMREX_D_DECL(       0,        0,        0));
        IntVect dom_hi(AMREX_D_DECL(n_cell-1, n_cell-1, n_cell-1));
        Box domain(dom_lo, dom_hi);

        // Initialize the boxarray "ba" from the single box "bx"
        ba.define(domain);
        // Break up boxarray "ba" into chunks no larger than "max_grid_size" along a direction
        ba.maxSize(max_grid_size);

       // This defines the physical box, [-1,1] in each direction.
        RealBox real_box({AMREX_D_DECL( 0.0, 0.0, 0.0)},
                         {AMREX_D_DECL( 1.0, 1.0, 1.0)});

        // This defines a Geometry object
        geom.define(domain,&real_box,CoordSys::cartesian,is_periodic.data());
    }

    // Nghost = number of ghost cells for each array
    int Nghost = 1;

    // Ncomp = number of components for each array
    int Ncomp;
    if(alg_type == 0){
       Ncomp = 1;
    } else {
       Ncomp = 2;
    }

    // How Boxes are distrubuted among MPI processes
    DistributionMapping dm(ba);


    /////////////////////////////////////////
    //Initialise rngs
    /////////////////////////////////////////

    if (restart < 0) {

        if (seed > 0) {
            // initializes the seed for C++ random number calls
            InitRandom(seed+ParallelDescriptor::MyProc(),
                       ParallelDescriptor::NProcs(),
                       seed+ParallelDescriptor::MyProc());
        } else if (seed == 0) {
            // initializes the seed for C++ random number calls based on the clock
            auto now = time_point_cast<nanoseconds>(system_clock::now());
            int randSeed = now.time_since_epoch().count();
            // broadcast the same root seed to all processors
            ParallelDescriptor::Bcast(&randSeed,1,ParallelDescriptor::IOProcessorNumber());
            InitRandom(randSeed+ParallelDescriptor::MyProc(),
                       ParallelDescriptor::NProcs(),
                       randSeed+ParallelDescriptor::MyProc());
        } else {
            Abort("Must supply non-negative seed");
        }

    }


    // we allocate two phi multifabs; one will store the old state, the other the new.
    MultiFab phi_old(ba, dm, Ncomp, Nghost);
    MultiFab phi_new(ba, dm, Ncomp, Nghost);

    GpuArray<Real,AMREX_SPACEDIM> dx = geom.CellSizeArray();

    init_phi(phi_new, geom,npts_scale, Ncomp);

    //Boundary conditions are assigned to phi_old such that the ghost cells at the boundary will
    //be filled to satisfy those conditions.
    Vector<BCRec> bc(phi_old.nComp());
    for (int n = 0; n < phi_old.nComp(); ++n)
    {
        for(int idim = 0; idim < AMREX_SPACEDIM; ++idim)
        {
            //Internal Dirichlet Periodic Boundary conditions, or bc_lo = bc_hi = 0
            if (bc_lo[idim] == BCType::int_dir) {
                bc[n].setLo(idim, BCType::int_dir);
            }
            //First Order Extrapolation for Neumann boundary conditions or bc_lo, bc_hi = 2
            else if (bc_lo[idim] == BCType::foextrap) {
                bc[n].setLo(idim, BCType::foextrap);
            }
            //External Dirichlet Boundary Condition, or bc_lo, bc_hi = 3
            else if(bc_lo[idim] == BCType::ext_dir) {
                bc[n].setLo(idim, BCType::ext_dir);
            }
            else {
                amrex::Abort("Invalid bc_lo");
            }

            //Internal Dirichlet Periodic Boundary conditions, or bc_lo = bc_hi = 0
            if (bc_hi[idim] == BCType::int_dir) {
                bc[n].setHi(idim, BCType::int_dir);
            }
            //First Order Extrapolation for Neumann boundary conditions or bc_lo, bc_hi = 2
            else if (bc_hi[idim] == BCType::foextrap) {
                bc[n].setHi(idim, BCType::foextrap);
            }
            //External Dirichlet Boundary Condition, or bc_lo, bc_hi = 3
            else if(bc_hi[idim] == BCType::ext_dir) {
                bc[n].setHi(idim, BCType::ext_dir);
            }
            else {
                amrex::Abort("Invalid bc_hi");
            }
        }
    }

    Real coeff = AMREX_D_TERM(   2./(dx[0]*dx[0]),
                               + 2./(dx[1]*dx[1]),
                               + 2./(dx[2]*dx[2]) );
    Real dt = cfl/(2.0*coeff);

    amrex::Print() << "dt = " << dt << " dx = " << dx[0] << std::endl;

    // time = starting time in the simulation
    Real time = 0.0;

    // Write a plotfile of the initial data if plot_int > 0 (plot_int was defined in the inputs file)
    if (plot_int > 0)
    {
        int n = 0;
        const std::string& pltfile = amrex::Concatenate("plt",n,6);
        WriteSingleLevelPlotfile(pltfile, phi_new, {"phi"}, geom, time, 0);
    }

    // build the flux multifabs
    Array<MultiFab, AMREX_SPACEDIM> flux;
    std::array< MultiFab, AMREX_SPACEDIM > stochFlux;
    for (int dir = 0; dir < AMREX_SPACEDIM; dir++)
    {
        // flux(dir) has one component, zero ghost cells, and is nodal in direction dir
        BoxArray edge_ba = ba;
        edge_ba.surroundingNodes(dir);
        flux[dir].define(edge_ba, dm, Ncomp, 0);
        stochFlux[dir].define(edge_ba, dm, Ncomp, 0);
    }

    AMREX_D_TERM(stochFlux[0].setVal(0.0);,
                 stochFlux[1].setVal(0.0);,
                 stochFlux[2].setVal(0.0););


    for (int n = 1; n <= nsteps; ++n)
    {
        MultiFab::Copy(phi_old, phi_new, 0, 0, Ncomp, 0);

        // new_phi = old_phi + dt * (something)
        advance(phi_old, phi_new, flux, stochFlux, dt, npts_scale, geom, bc, Ncomp);
        time = time + dt;

        // Tell the I/O Processor to write out which step we're doing
        amrex::Print() << "Advanced step " << n << "\n";

        // Write a plotfile of the current data (plot_int was defined in the inputs file)
        if (plot_int > 0 && n%plot_int == 0)
        {
            const std::string& pltfile = amrex::Concatenate("plt",n,6);
            WriteSingleLevelPlotfile(pltfile, phi_new, {"phi"}, geom, time, n);

//            amrex::Real Ephi=0.;
//            amrex::Real Ephi2=0.;
            Vector<Real> Ephi(3,0.);
            Real Ephimin = npts_scale;
            for ( MFIter mfi(phi_old); mfi.isValid(); ++mfi )
            {
                const Box& vbx = mfi.validbox();
                const auto lo = amrex::lbound(vbx);
                const auto hi = amrex::ubound(vbx);

               auto const& phiNew = phi_new.array(mfi);

               for (auto k = lo.z; k <= hi.z; ++k) {
               for (auto j = lo.y; j <= hi.y; ++j) {
               for (auto i = lo.x; i <= hi.x; ++i) {
                   Ephi[0] += phiNew(i,j,k);
                   Ephi[1] += phiNew(i,j,k)*phiNew(i,j,k);
                   Ephimin = std::min(Ephimin,phiNew(i,j,k));
                   Ephi[2] += (phiNew(i,j,k) < 0) ? 1. : 0. ;
                }
                }
                }


            }

            const int IOProc = ParallelDescriptor::IOProcessorNumber();
            ParallelDescriptor::ReduceRealSum(Ephi.dataPtr(),3);
            ParallelDescriptor::ReduceRealMin(Ephimin);
            amrex::Real scale = n_cell*n_cell;
            amrex::Real scale2 =  AMREX_D_TERM( dx[0],
                               * dx[1],
                               * dx[2] );

            amrex::Print() << "phi variance = " << Ephi[1]/scale - (Ephi[0]*Ephi[0]
                            /(scale*scale)) << std::endl;
            amrex::Print() << "phi integral = " << Ephi[0]*scale2 << std::endl;
            amrex::Print() << "phi min = " << Ephimin << " Number of negative points " << Ephi[2] <<  std::endl;

        }


    }

    // Call the timer again and compute the maximum difference between the start time and stop time
    //   over all processors
    auto stop_time = ParallelDescriptor::second() - strt_time;
    const int IOProc = ParallelDescriptor::IOProcessorNumber();
    ParallelDescriptor::ReduceRealMax(stop_time,IOProc);

    // Tell the I/O Processor to write out the "run time"
    amrex::Print() << "Run time = " << stop_time << std::endl;
}
