
#include "common_functions.H"
#include "rng_functions.H"

#include <AMReX_VisMF.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_MultiFabUtil.H>

#include "chrono"

using namespace std::chrono;


// argv contains the name of the inputs file entered at the command line
void main_driver(const char* argv)
{
  
    BL_PROFILE_VAR("main_driver()",main_driver);

    // store the current time so we can later compute total run time.
    Real strt_time = ParallelDescriptor::second();

    //=============================================================
    // Initialization
    //=============================================================
    
    std::string inputs_file = argv;

    InitializeCommonNamespace();

    // copy objects from namelist to temporary variables
    // equilibrium height "h0" is stored in rho0
    // surface tension "gamma" is stored in h_bar
    Real h0 = rho0;
    Real gamma = h_bar;

    // algorithm_type = 0 (2D, default)
    // algorithm_type = 1 (1D)
    int do_1d = (algorithm_type == 1) ? 1 : 0;

    // set fac_1d = 1. for 2D
    // set fac_1d = 0. for 1D
    Real fac_1d = (do_1d == 1) ? 0. : 1.;
    
    /////////////////////////////////////////
    // Initialize random number seed on all processors/GPUs
    /////////////////////////////////////////

    int mySeed;
    
    if (seed > 0) {
        // use seed from inputs file
        mySeed = seed;

    } else if (seed == 0) {
        // use a clock-based seed
        auto now = time_point_cast<nanoseconds>(system_clock::now());
        mySeed = now.time_since_epoch().count();
        // broadcast the same root seed to all processors
        ParallelDescriptor::Bcast(&mySeed,1,ParallelDescriptor::IOProcessorNumber());
    } else {
        Abort("Must supply non-negative seed");
    }

    // initialize the seed for C++ random number calls
    InitRandom(mySeed+ParallelDescriptor::MyProc(),
               ParallelDescriptor::NProcs(),
               mySeed+ParallelDescriptor::MyProc());

    /////////////////////////////////////////
    
    // is the problem periodic?
    Vector<int> is_periodic(AMREX_SPACEDIM,1);  // set to 1 (periodic) by default
        
    // This defines the physical box, [-1,1] in each direction.
    RealBox real_box({AMREX_D_DECL(prob_lo[0],prob_lo[1],prob_lo[2])},
                     {AMREX_D_DECL(prob_hi[0],prob_hi[1],prob_hi[2])});
        
    IntVect dom_lo(AMREX_D_DECL(           0,            0,            0));
    IntVect dom_hi(AMREX_D_DECL(n_cells[0]-1, n_cells[1]-1, n_cells[2]-1));
    Box domain(dom_lo, dom_hi);

    Geometry geom(domain,&real_box,CoordSys::cartesian,is_periodic.data());

    const GpuArray<Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();

    Real dVol = dx[0]*dx[1];
    
    // make BoxArray and Geometry
    BoxArray ba;

    // how boxes are distrubuted among MPI processes
    DistributionMapping dmap;    

    // Initialize the boxarray "ba" from the single box "bx"
    ba.define(domain);

    // Break up boxarray "ba" into chunks no larger than "max_grid_size" along a direction
    // note we are converting "Vector<int> max_grid_size" to an IntVect
    ba.maxSize(IntVect(max_grid_size));

    // define DistributionMapping
    dmap.define(ba);

    MultiFab height(ba, dmap, 1, 1);
    MultiFab Laph  (ba, dmap, 1, 1);

    std::array< MultiFab, AMREX_SPACEDIM > hface;
    std::array< MultiFab, AMREX_SPACEDIM > gradLh;
    std::array< MultiFab, AMREX_SPACEDIM > flux;
    std::array< MultiFab, AMREX_SPACEDIM > randface;
    
    AMREX_D_TERM(hface[0]   .define(convert(ba,nodal_flag_x), dmap, 1, 0);,
                 hface[1]   .define(convert(ba,nodal_flag_y), dmap, 1, 0);,
                 hface[2]   .define(convert(ba,nodal_flag_z), dmap, 1, 0););
    AMREX_D_TERM(gradLh[0]  .define(convert(ba,nodal_flag_x), dmap, 1, 0);,
                 gradLh[1]  .define(convert(ba,nodal_flag_y), dmap, 1, 0);,
                 gradLh[2]  .define(convert(ba,nodal_flag_z), dmap, 1, 0););
    AMREX_D_TERM(flux[0]    .define(convert(ba,nodal_flag_x), dmap, 1, 0);,
                 flux[1]    .define(convert(ba,nodal_flag_y), dmap, 1, 0);,
                 flux[2]    .define(convert(ba,nodal_flag_z), dmap, 1, 0););
    AMREX_D_TERM(randface[0].define(convert(ba,nodal_flag_x), dmap, 1, 0);,
                 randface[1].define(convert(ba,nodal_flag_y), dmap, 1, 0);,
                 randface[2].define(convert(ba,nodal_flag_z), dmap, 1, 0););

    // initialize height
    height.setVal(h0);
    
    // Physical time constant for dimensional time
    Real t0 = 3.0*visc_coef*h0/gamma;

    // constant factor in noise term
    Real ConstNoise = 2.*k_B*T_init[0] / (3.*visc_coef);
    Real Const3dx = gamma / (3.*visc_coef);
    
    Real time = 0.;
    Real dt = 0.1 * (t0/std::pow(h0,4)) * std::pow(dx[0],4) / 16.;
    
    if (plot_int > 0)
    {
        int step = 0;
        const std::string& pltfile = amrex::Concatenate("plt",step,5);
        WriteSingleLevelPlotfile(pltfile, height, {"height"}, geom, time, 0);
    }
    
    // Time stepping loop
    for(int istep=1; istep<=max_step; ++istep) {

        // timer
        Real step_strt_time = ParallelDescriptor::second();
    
        // fill random numbers
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            MultiFabFillRandom(randface[d], 0, variance_coef_mass, geom);
        }

        // compute Laph
        for ( MFIter mfi(Laph,TilingIfNotGPU()); mfi.isValid(); ++mfi ) {
        
            const Box& bx = mfi.tilebox();
            const Array4<Real> & L = Laph.array(mfi);
            const Array4<Real> & h = height.array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                L(i,j,k) =          ( h(i-1,j,k) - 2.*h(i,j,k) + h(i+1,j,k) ) / (dx[0]*dx[0])
                         + fac_1d * ( h(i,j-1,k) - 2.*h(i,j,k) + h(i,j+1,k) ) / (dx[1]*dx[1]);
            });
        }
        Laph.FillBoundary(geom.periodicity());
    
        // compute hface and grad(Laph)
        for ( MFIter mfi(height,TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

            AMREX_D_TERM(const Box & bx_x = mfi.nodaltilebox(0);,
                         const Box & bx_y = mfi.nodaltilebox(1);,
                         const Box & bx_z = mfi.nodaltilebox(2););
        
            AMREX_D_TERM(const Array4<Real> & hfacex = hface[0].array(mfi);,
                         const Array4<Real> & hfacey = hface[1].array(mfi);,
                         const Array4<Real> & hfacez = hface[2].array(mfi););
        
            AMREX_D_TERM(const Array4<Real> & gradLhx = gradLh[0].array(mfi);,
                         const Array4<Real> & gradLhy = gradLh[1].array(mfi);,
                         const Array4<Real> & gradLhz = gradLh[2].array(mfi););
            
            const Array4<Real> & L = Laph.array(mfi);
            const Array4<Real> & h = height.array(mfi);

            amrex::ParallelFor(bx_x, bx_y,
                               [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                hfacex(i,j,k) = 0.5*( h(i-1,j,k) + h(i,j,k) );
                gradLhx(i,j,k) = ( L(i,j,k) - L(i-1,j,k) ) / dx[0];
            },
                               [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                hfacey(i,j,k) = 0.5*( h(i,j-1,k) + h(i,j,k) );
                gradLhy(i,j,k) = ( L(i,j,k) - L(i,j-1,k) ) / dx[1];
            });

        }

        // compute flux
    
        // compute hface and grad(Laph)
        for ( MFIter mfi(height,TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

            AMREX_D_TERM(const Box & bx_x = mfi.nodaltilebox(0);,
                         const Box & bx_y = mfi.nodaltilebox(1);,
                         const Box & bx_z = mfi.nodaltilebox(2););
        
            AMREX_D_TERM(const Array4<Real> & hfacex = hface[0].array(mfi);,
                         const Array4<Real> & hfacey = hface[1].array(mfi);,
                         const Array4<Real> & hfacez = hface[2].array(mfi););
        
            AMREX_D_TERM(const Array4<Real> & gradLhx = gradLh[0].array(mfi);,
                         const Array4<Real> & gradLhy = gradLh[1].array(mfi);,
                         const Array4<Real> & gradLhz = gradLh[2].array(mfi););
        
            AMREX_D_TERM(const Array4<Real> & fluxx = flux[0].array(mfi);,
                         const Array4<Real> & fluxy = flux[1].array(mfi);,
                         const Array4<Real> & fluxz = flux[2].array(mfi););
        
            AMREX_D_TERM(const Array4<Real> & randfacex = randface[0].array(mfi);,
                         const Array4<Real> & randfacey = randface[1].array(mfi);,
                         const Array4<Real> & randfacez = randface[2].array(mfi););

            amrex::ParallelFor(bx_x, bx_y,
                               [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                fluxx(i,j,k) = std::sqrt(ConstNoise*std::pow(hfacex(i,j,k),3.) / (dt*dVol)) * randfacex(i,j,k)
                    + Const3dx * std::pow(hfacex(i,j,k),3.)*gradLhx(i,j,k);
            },
                               [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                fluxy(i,j,k) = fac_1d * (std::sqrt(ConstNoise*std::pow(hfacey(i,j,k),3.) / (dt*dVol)) * randfacey(i,j,k)
                                         + Const3dx * std::pow(hfacey(i,j,k),3.)*gradLhy(i,j,k) );
            });

        }

        // update height using forward Euler
        for ( MFIter mfi(height,TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

            const Box& bx = mfi.tilebox();
        
            AMREX_D_TERM(const Array4<Real> & fluxx = flux[0].array(mfi);,
                         const Array4<Real> & fluxy = flux[1].array(mfi);,
                         const Array4<Real> & fluxz = flux[2].array(mfi););
            
            const Array4<Real> & h = height.array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                h(i,j,k) -= dt * ( (fluxx(i+1,j,k) - fluxx(i,j,k)) / dx[0]
                                  +(fluxy(i,j+1,k) - fluxy(i,j,k)) / dx[1] );
            });

        }
        height.FillBoundary(geom.periodicity());

        time += dt;
        
        if (plot_int > 0 && istep%plot_int == 0)
        {
            const std::string& pltfile = amrex::Concatenate("plt",istep,5);
            WriteSingleLevelPlotfile(pltfile, height, {"height"}, geom, time, 0);
        }
        
        Real step_stop_time = ParallelDescriptor::second() - step_strt_time;
        ParallelDescriptor::ReduceRealMax(step_stop_time);
        amrex::Print() << "Time step " << istep << " complted in " << step_stop_time
                       << " seconds " << std::endl;
    
        // MultiFab memory usage
        const int IOProc = ParallelDescriptor::IOProcessorNumber();

        amrex::Long min_fab_megabytes  = amrex::TotalBytesAllocatedInFabsHWM()/1048576;
        amrex::Long max_fab_megabytes  = min_fab_megabytes;

        ParallelDescriptor::ReduceLongMin(min_fab_megabytes, IOProc);
        ParallelDescriptor::ReduceLongMax(max_fab_megabytes, IOProc);

        amrex::Print() << "High-water FAB megabyte spread across MPI nodes: ["
                       << min_fab_megabytes << " ... " << max_fab_megabytes << "]\n";

        min_fab_megabytes  = amrex::TotalBytesAllocatedInFabs()/1048576;
        max_fab_megabytes  = min_fab_megabytes;

        ParallelDescriptor::ReduceLongMin(min_fab_megabytes, IOProc);
        ParallelDescriptor::ReduceLongMax(max_fab_megabytes, IOProc);

        amrex::Print() << "Curent     FAB megabyte spread across MPI nodes: ["
                       << min_fab_megabytes << " ... " << max_fab_megabytes << "]\n";
        
    }

    // Call the timer again and compute the maximum difference between the start time 
    // and stop time over all processors
    Real stop_time = ParallelDescriptor::second() - strt_time;
    ParallelDescriptor::ReduceRealMax(stop_time);
    amrex::Print() << "Run time = " << stop_time << " seconds " << std::endl;

}
