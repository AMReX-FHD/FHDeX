#include "thinfilm_functions.H"
#include "common_functions.H"
#include "rng_functions.H"

#include <AMReX_VisMF.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_MultiFabUtil.H>

#include "chrono"

#ifdef AMREX_USE_CUDA
#include <cufft.h>
#else
#include <fftw3.h>
#include <fftw3-mpi.h>
#endif

using namespace std::chrono;

// argv contains the name of the inputs file entered at the command line
void main_driver(const char* argv)
{

    BL_PROFILE_VAR("main_driver()",main_driver);

    if (AMREX_SPACEDIM != 2) {
        Abort("Must build/run with DIM=2");
    }

    // store the current time so we can later compute total run time.
    Real strt_time = ParallelDescriptor::second();

    //=============================================================
    // Initialization
    //=============================================================

    std::string inputs_file = argv;

    InitializeCommonNamespace();
    InitializeThinfilmNamespace();

    // algorithm_type = 0 (1D in x; each row is an independent trial)
    // algorithm_type = 1 (1D in y; each col is an independent trial)
    // algorithm_type = 2 (2D)
    int do_1d_x = (algorithm_type == 0) ? 1 : 0;
    int do_1d_y = (algorithm_type == 1) ? 1 : 0;

    Real y_flux_fac = (do_1d_x == 1) ? 0. : 1.; // 1D in x; set y fluxes to zero
    Real x_flux_fac = (do_1d_y == 1) ? 0. : 1.; // 1D in y; set x fluxes to zero

    if (do_1d_x) {
        if (n_cells[0] > max_grid_size[0]) {
            Abort("For 1D-x mode, max_grid_size[0] must be >= n_cells[0]");
        }
    }
    if (do_1d_y) {
        if (n_cells[1] > max_grid_size[1]) {
            Abort("For 1D-y mode, max_grid_size[1] must be >= n_cells[1]");
        }
    }

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
    Vector<int> is_periodic(AMREX_SPACEDIM,0);  // set to 0 (non-periodic) by default

    // boundary conditions controlled by bc_mass_lo/hi
    // -1 = periodic
    //  0 = 90 degree contact angle, dh/dn=0
    //  1 = pinned, h=thinfilm_h0
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        if (bc_mass_lo[d] == -1 && bc_mass_hi[d] != -1 ||
            bc_mass_lo[d] != -1 && bc_mass_hi[d] == -1) {
            Abort("Error: specified periodic and non-periodic in the same directionality");
        }
        // set periodic
        if (bc_mass_lo[d] == -1) {
            is_periodic[d] = 1;
        }
    }

    // This defines the physical box, [-1,1] in each direction.
    RealBox real_box({prob_lo[0],prob_lo[1]},
                     {prob_hi[0],prob_hi[1]});

    IntVect dom_lo(           0,            0);
    IntVect dom_hi(n_cells[0]-1, n_cells[1]-1);
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
    MultiFab disjoining(ba, dmap, 1, 1);
    Laph.setVal(0.); // prevent intermediate NaN calculations behind physical boundaries
    disjoining.setVal(0.);

    // for statsitics

    // delta h
    MultiFab dheight    (ba, dmap, 1, 0);

    // variance
    MultiFab dheight2sum(ba, dmap, 1, 0);
    MultiFab dheight2avg(ba, dmap, 1, 0);
    dheight2sum.setVal(0.);

    // spatial correlation
    MultiFab dheightstarsum(ba, dmap, 1, 0);
    MultiFab dheightstaravg(ba, dmap, 1, 0);
    dheightstarsum.setVal(0.);

    std::array< MultiFab, AMREX_SPACEDIM > hface;
    AMREX_D_TERM(hface[0]   .define(convert(ba,nodal_flag_x), dmap, 1, 0);,
                 hface[1]   .define(convert(ba,nodal_flag_y), dmap, 1, 0);,
                 hface[2]   .define(convert(ba,nodal_flag_z), dmap, 1, 0););

    std::array< MultiFab, AMREX_SPACEDIM > gradh;
    AMREX_D_TERM(gradh[0]   .define(convert(ba,nodal_flag_x), dmap, 1, 0);,
                 gradh[1]   .define(convert(ba,nodal_flag_y), dmap, 1, 0);,
                 gradh[2]   .define(convert(ba,nodal_flag_z), dmap, 1, 0););

    std::array< MultiFab, AMREX_SPACEDIM > gradLaph;
    AMREX_D_TERM(gradLaph[0].define(convert(ba,nodal_flag_x), dmap, 1, 0);,
                 gradLaph[1].define(convert(ba,nodal_flag_y), dmap, 1, 0);,
                 gradLaph[2].define(convert(ba,nodal_flag_z), dmap, 1, 0););

    std::array< MultiFab, AMREX_SPACEDIM > gradDisjoining;
    AMREX_D_TERM(gradDisjoining[0].define(convert(ba,nodal_flag_x), dmap, 1, 0);,
                 gradDisjoining[1].define(convert(ba,nodal_flag_y), dmap, 1, 0);,
                 gradDisjoining[2].define(convert(ba,nodal_flag_z), dmap, 1, 0););

    std::array< MultiFab, AMREX_SPACEDIM > flux;
    AMREX_D_TERM(flux[0]    .define(convert(ba,nodal_flag_x), dmap, 1, 0);,
                 flux[1]    .define(convert(ba,nodal_flag_y), dmap, 1, 0);,
                 flux[2]    .define(convert(ba,nodal_flag_z), dmap, 1, 0););

    std::array< MultiFab, AMREX_SPACEDIM > randface;
    AMREX_D_TERM(randface[0].define(convert(ba,nodal_flag_x), dmap, 1, 0);,
                 randface[1].define(convert(ba,nodal_flag_y), dmap, 1, 0);,
                 randface[2].define(convert(ba,nodal_flag_z), dmap, 1, 0););

    // initialize height
    height.setVal(thinfilm_h0);

    // update height using forward Euler
    if (thinfilm_pertamp > 0.) {
        for ( MFIter mfi(height,TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

            const Box& bx = mfi.tilebox();
            const Array4<Real> & h = height.array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                h(i,j,k) += thinfilm_pertamp * thinfilm_h0*std::cos(2.*M_PI*(i+.5)*dx[0]/(prob_hi[0]-prob_lo[0]));
            });

        }
    }

    height.FillBoundary(geom.periodicity());

    // Physical time constant for dimensional time
    Real t0 = 3.0*visc_coef*thinfilm_h0/thinfilm_gamma;

    // constant factor in noise term
    Real ConstNoise = 2.*k_B*T_init[0] / (3.*visc_coef);
    Real Const3dx = thinfilm_gamma / (3.*visc_coef);
    Real Const3dx_nogamma = 1. / (3.*visc_coef);

    Real time = 0.;

    Real dx_min;
    if (algorithm_type == 0) {
        dx_min = dx[0];
    } else if (algorithm_type == 1) {
        dx_min = dx[1];
    } else if (algorithm_type == 2) {
        dx_min = std::min(dx[0],dx[1]);
    }

    //////////////////////////////////////
    // for ASCII writeout and 1D ffts in 1D mode
    BoxArray ba_onegrid(domain);
    DistributionMapping dmap_onegrid(ba_onegrid);

    MultiFab height_onegrid(ba_onegrid, dmap_onegrid, 1, 0);

    // data structures to take 1D ffts in 1D mode
    MultiFab height_flat_onegrid;
    MultiFab fft_sum;
    MultiFab fft_avg;

    Geometry geom_fft;

    long npts_flat;

    if (algorithm_type == 0 || algorithm_type == 1) {
        // make the "1D" MultiFab x-based, regardless of whether problem is x or y-based
        IntVect dom_lo_flat(0,0);
        IntVect dom_hi_flat(n_cells[algorithm_type]-1,0);

        Box domain_flat(dom_lo_flat,dom_hi_flat);
        BoxArray ba_flat_onegrid(domain_flat);

        npts_flat = domain_flat.numPts();

        // use same dmap as height_onegrid so we can easily copy a strip into here
        height_flat_onegrid.define(ba_flat_onegrid,dmap_onegrid,1,0);

        // fft contains n_cell/2+1 points (index space 0 to n_cell/2
        IntVect dom_lo_fft(0,0);
        IntVect dom_hi_fft(n_cells[algorithm_type]/2,0);
        Box domain_fft(dom_lo_fft,dom_hi_fft);
        BoxArray ba_fft(domain_fft);

        // magnitude of fft, sum and average
        // use same dmap as height_onegrid
        fft_sum.define(ba_fft,dmap_onegrid,1,0);
        fft_avg.define(ba_fft,dmap_onegrid,1,0);

        fft_sum.setVal(0.);

        geom_fft.define(domain_fft,&real_box,CoordSys::cartesian,is_periodic.data());
    }

    //////////////////////////////////////

    // time step
    Real dt = cfl * (t0/std::pow(thinfilm_h0,4)) * std::pow(dx_min,4) / 16.;

    int stats_count = 0;

    if (plot_int > 0)
    {
        int step = 0;
        const std::string& pltfile = amrex::Concatenate("plt",step,10);
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

        // compute hface and gradh
        for ( MFIter mfi(height,TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

            AMREX_D_TERM(const Box & bx_x = mfi.nodaltilebox(0);,
                         const Box & bx_y = mfi.nodaltilebox(1);,
                         const Box & bx_z = mfi.nodaltilebox(2););

            AMREX_D_TERM(const Array4<Real> & hfacex = hface[0].array(mfi);,
                         const Array4<Real> & hfacey = hface[1].array(mfi);,
                         const Array4<Real> & hfacez = hface[2].array(mfi););

            AMREX_D_TERM(const Array4<Real> & gradhx = gradh[0].array(mfi);,
                         const Array4<Real> & gradhy = gradh[1].array(mfi);,
                         const Array4<Real> & gradhz = gradh[2].array(mfi););

            const Array4<Real> & h = height.array(mfi);

            amrex::ParallelFor(bx_x, bx_y,
                               [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                hfacex(i,j,k) = 0.5*( h(i-1,j,k) + h(i,j,k) );
                gradhx(i,j,k) = ( h(i,j,k) - h(i-1,j,k) ) / dx[0];
            },
                               [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                hfacey(i,j,k) = 0.5*( h(i,j-1,k) + h(i,j,k) );
                gradhy(i,j,k) = ( h(i,j,k) - h(i,j-1,k) ) / dx[1];
            });

            // fix gradh at physical boundaries
            // no need to fix hface since we zero the total fluxes later

            // dh/dx = 0, x-faces
            if (bc_mass_lo[0] == 0 || bc_mass_hi[0] == 0) {
                amrex::ParallelFor(bx_x, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    // lo
                    if (bc_mass_lo[0] == 0 && i==0)          gradhx(i,j,k) = 0.;
                    // hi
                    if (bc_mass_hi[0] == 0 && i==n_cells[0]) gradhx(i,j,k) = 0.;
                });
            }

            // h = thinfilm_h0, xfaces
            if (bc_mass_lo[0] == 1 || bc_mass_hi[0] == 1) {
                amrex::ParallelFor(bx_x, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    // lo
                    if (bc_mass_lo[0] == 1 && i==0)          gradhx(i,j,k) = (h(i,j,k) - thinfilm_h0) / (0.5*dx[0]);
                    // hi
                    if (bc_mass_hi[0] == 1 && i==n_cells[0]) gradhx(i,j,k) = (thinfilm_h0 - h(i-1,j,k)) / (0.5*dx[0]);
                });
            }

            // dh/dx = 0, y-faces
            if (bc_mass_lo[1] == 0 || bc_mass_hi[1] == 0) {
                amrex::ParallelFor(bx_y, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    // lo
                    if (bc_mass_lo[1] == 0 && j==0)          gradhy(i,j,k) = 0.;
                    // hi
                    if (bc_mass_hi[1] == 0 && j==n_cells[1]) gradhy(i,j,k) = 0.;
                });
            }

            // h = thinfilm_h0, yfaces
            if (bc_mass_lo[1] == 1 || bc_mass_hi[1] == 1) {
                amrex::ParallelFor(bx_y, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    // lo
                    if (bc_mass_lo[1] == 1 && j==0)          gradhy(i,j,k) = (h(i,j,k) - thinfilm_h0) / (0.5*dx[1]);
                    // hi
                    if (bc_mass_hi[1] == 1 && j==n_cells[1]) gradhy(i,j,k) = (thinfilm_h0 - h(i,j-1,k)) / (0.5*dx[1]);
                });
            }

        }

        // compute Laph and disjoining
        for ( MFIter mfi(Laph,TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

            const Box& bx = mfi.tilebox();

            const Array4<Real> & L = Laph.array(mfi);

            const Array4<Real> & h = height.array(mfi);

            const Array4<Real> & Disjoining = disjoining.array(mfi);

            AMREX_D_TERM(const Array4<Real> & gradhx = gradh[0].array(mfi);,
                         const Array4<Real> & gradhy = gradh[1].array(mfi);,
                         const Array4<Real> & gradhz = gradh[2].array(mfi););

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                L(i,j,k) = x_flux_fac * (gradhx(i+1,j,k) - gradhx(i,j,k)) / dx[0]
                         + y_flux_fac * (gradhy(i,j+1,k) - gradhy(i,j,k)) / dx[1];

                Disjoining(i,j,k) = thinfilm_hamaker / (6.*M_PI*std::pow(h(i,j,k),3.));
            });
        }
        Laph.FillBoundary(geom.periodicity());
        disjoining.FillBoundary(geom.periodicity());

        // compute gradLaph and gradDisjoining
        for ( MFIter mfi(height,TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

            AMREX_D_TERM(const Box & bx_x = mfi.nodaltilebox(0);,
                         const Box & bx_y = mfi.nodaltilebox(1);,
                         const Box & bx_z = mfi.nodaltilebox(2););

            AMREX_D_TERM(const Array4<Real> & gradLaphx = gradLaph[0].array(mfi);,
                         const Array4<Real> & gradLaphy = gradLaph[1].array(mfi);,
                         const Array4<Real> & gradLaphz = gradLaph[2].array(mfi););

            AMREX_D_TERM(const Array4<Real> & gradDisjoiningx = gradDisjoining[0].array(mfi);,
                         const Array4<Real> & gradDisjoiningy = gradDisjoining[1].array(mfi);,
                         const Array4<Real> & gradDisjoiningz = gradDisjoining[2].array(mfi););

            const Array4<Real> & L = Laph.array(mfi);

            const Array4<Real> & Disjoining = disjoining.array(mfi);

            amrex::ParallelFor(bx_x, bx_y,
                               [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                gradLaphx(i,j,k) = ( L(i,j,k) - L(i-1,j,k) ) / dx[0];
                gradDisjoiningx(i,j,k) = ( Disjoining(i,j,k) -Disjoining(i-1,j,k) ) / dx[0];
            },
                               [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                gradLaphy(i,j,k) = ( L(i,j,k) - L(i,j-1,k) ) / dx[1];
                gradDisjoiningy(i,j,k) = ( Disjoining(i,j,k) -Disjoining(i,j-1,k) ) / dx[1];
            });
        }

        // compute flux
        for ( MFIter mfi(height,TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

            AMREX_D_TERM(const Box & bx_x = mfi.nodaltilebox(0);,
                         const Box & bx_y = mfi.nodaltilebox(1);,
                         const Box & bx_z = mfi.nodaltilebox(2););

            AMREX_D_TERM(const Array4<Real> & hfacex = hface[0].array(mfi);,
                         const Array4<Real> & hfacey = hface[1].array(mfi);,
                         const Array4<Real> & hfacez = hface[2].array(mfi););

            AMREX_D_TERM(const Array4<Real> & gradLaphx = gradLaph[0].array(mfi);,
                         const Array4<Real> & gradLaphy = gradLaph[1].array(mfi);,
                         const Array4<Real> & gradLaphz = gradLaph[2].array(mfi););

            AMREX_D_TERM(const Array4<Real> & fluxx = flux[0].array(mfi);,
                         const Array4<Real> & fluxy = flux[1].array(mfi);,
                         const Array4<Real> & fluxz = flux[2].array(mfi););

            AMREX_D_TERM(const Array4<Real> & randfacex = randface[0].array(mfi);,
                         const Array4<Real> & randfacey = randface[1].array(mfi);,
                         const Array4<Real> & randfacez = randface[2].array(mfi););

            AMREX_D_TERM(const Array4<Real> & gradDisjoiningx = gradDisjoining[0].array(mfi);,
                         const Array4<Real> & gradDisjoiningy = gradDisjoining[1].array(mfi);,
                         const Array4<Real> & gradDisjoiningz = gradDisjoining[2].array(mfi););

            amrex::ParallelFor(bx_x, bx_y,
                               [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                fluxx(i,j,k) = x_flux_fac * (
                               std::sqrt(ConstNoise*std::pow(hfacex(i,j,k),3.) / (dt*dVol)) * randfacex(i,j,k)
                               + Const3dx         * std::pow(hfacex(i,j,k),3.)*gradLaphx(i,j,k)
                               + Const3dx_nogamma * std::pow(hfacex(i,j,k),3.)*gradDisjoiningx(i,j,k));
            },
                               [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                fluxy(i,j,k) = y_flux_fac * (
                               std::sqrt(ConstNoise*std::pow(hfacey(i,j,k),3.) / (dt*dVol)) * randfacey(i,j,k)
                               + Const3dx         * std::pow(hfacey(i,j,k),3.)*gradLaphy(i,j,k)
                               + Const3dx_nogamma * std::pow(hfacex(i,j,k),3.)*gradDisjoiningy(i,j,k) );
            });

            // lo x-faces
            if (bc_mass_lo[0] == 0 || bc_mass_lo[0] == 1) {
                amrex::ParallelFor(bx_x, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (i==0) fluxx(i,j,k) = 0.;
                });
            }

            // hi x-faces
            if (bc_mass_hi[0] == 0 || bc_mass_hi[0] == 1) {
                amrex::ParallelFor(bx_x, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (i==n_cells[0]) fluxx(i,j,k) = 0.;
                });
            }

            // lo y-faces
            if (bc_mass_lo[1] == 0 || bc_mass_lo[1] == 1) {
                amrex::ParallelFor(bx_y, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (j==0) fluxy(i,j,k) = 0.;
                });
            }

            // hi y-faces
            if (bc_mass_hi[1] == 0 || bc_mass_hi[1] == 1) {
                amrex::ParallelFor(bx_y, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (j==n_cells[1]) fluxy(i,j,k) = 0.;
                });
            }
        }

        // update height using forward Euler
        for ( MFIter mfi(height,TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

            const Box& bx = mfi.tilebox();

            AMREX_D_TERM(const Array4<Real> & fluxx = flux[0].array(mfi);,
                         const Array4<Real> & fluxy = flux[1].array(mfi);,
                         const Array4<Real> & fluxz = flux[2].array(mfi););

            const Array4<Real> & h = height.array(mfi);

            // amrex::Print{} << "HEIGHT " << time << " " << h(0,0,0) << " " << h(31,0,0) << std::endl;

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {

                h(i,j,k) -= dt * ( (fluxx(i+1,j,k) - fluxx(i,j,k)) / dx[0]
                                  +(fluxy(i,j+1,k) - fluxy(i,j,k)) / dx[1]);
            });

        }
        height.FillBoundary(geom.periodicity());

        time += dt;

        // plotfile
        if (plot_int > 0 && istep%plot_int == 0)
        {
            const std::string& pltfile = amrex::Concatenate("plt",istep,10);
            WriteSingleLevelPlotfile(pltfile, height, {"height"}, geom, time, 0);

            // ASCII

            // copy distributed data into 1D data
            height_onegrid.ParallelCopy(height, 0, 0, 1);

            amrex::Real sumh = 0.;

            for ( MFIter mfi(height_onegrid,false); mfi.isValid(); ++mfi ) {

                std::ofstream hstream;
                std::string heightBaseName = "height";

                std::string heightName = Concatenate(heightBaseName,istep,10);
                heightName += ".txt";

                hstream.open(heightName);

                const Box& bx = mfi.validbox();
                const auto lo = amrex::lbound(bx);
                const auto hi = amrex::ubound(bx);

                const Array4<Real>& mfdata = height_onegrid.array(mfi);

                for (auto j = lo.y; j <= hi.y; ++j) {
                for (auto i = lo.x; i <= hi.x; ++i) {
                    hstream << std::setprecision(15) << mfdata(i,j,0) << " ";
                    sumh += mfdata(i,j,0);
                    if(j==0 && i==0){
                       amrex::Print{} << "HEIGHT " << time << " " << mfdata(0,0,0) << " " << mfdata(31,0,0) << std::endl;
                    }
                }
                hstream << "\n";
                }
                amrex::Print{} << "SUM " << sumh << std::endl;

            } // end MFIter

        }

        // statistics
        if (istep > n_steps_skip) {

            ++stats_count;
            Real stats_count_inv = 1./stats_count;

            ///////////////////////
            // variance

            // compute delta h
            dheight.setVal(thinfilm_h0);
            MultiFab::Subtract(dheight, height, 0, 0, 1, 0);

            // compute delta h * delta h
            MultiFab::Multiply(dheight, dheight, 0, 0, 1, 0);

            // compute sum of delta h * delta h
            MultiFab::Add(dheight2sum, dheight, 0, 0, 1, 0);

            // write variance to plotfile
            if (plot_int > 0 && istep%plot_int == 0)
            {
                // compute average of delta h * delta h
                MultiFab::Copy(dheight2avg, dheight2sum, 0, 0, 1, 0);
                dheight2avg.mult( stats_count_inv, 0, 1, 0);

                const std::string& pltfile = amrex::Concatenate("var",istep,10);
                WriteSingleLevelPlotfile(pltfile, dheight2avg, {"var"}, geom, time, 0);
            }

            ///////////////////////
            // spatial correlation

            if (do_1d_x || do_1d_y) {
                // for 1D mode, find delta h at each row (column) at icorr (jcorr) for x-mode (y-mode)

                // increment delta h^* * delta h
                for ( MFIter mfi(height,TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

                    const Box& bx = mfi.tilebox();

                    const Array4<Real> & h = height.array(mfi);
                    const Array4<Real> & dhstar = dheightstarsum.array(mfi);

                    if (do_1d_x) {
                        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                        {
                            dhstar(i,j,k) += ( h(thinfilm_icorr,j,k)-thinfilm_h0 ) * ( h(i,j,k)-thinfilm_h0 );
                        });
                    } else {
                        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                        {
                            dhstar(i,j,k) += ( h(i,thinfilm_jcorr,k)-thinfilm_h0 ) * ( h(i,j,k)-thinfilm_h0 );
                        });
                    }
                }

            } else {
                // for 2D mode, find delta h at point icorr,jcorr and compute correlation
                Real h_local = 0.;

                for ( MFIter mfi(height,TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

                    const Box& bx = mfi.tilebox();

                    const Array4<Real> & h = height.array(mfi);

                    if (bx.strictly_contains(IntVect(thinfilm_icorr,thinfilm_jcorr))) {
                        h_local += h(thinfilm_icorr,thinfilm_jcorr,0);
                    }

                }
                ParallelDescriptor::ReduceRealSum(h_local);

                // compute delta h^*
                h_local -= thinfilm_h0;

                // compute delta h
                dheight.setVal(thinfilm_h0);
                MultiFab::Subtract(dheight, height, 0, 0, 1, 0);

                // compute delta h^* * delta h
                dheight.mult(h_local, 0, 1, 0);

                // compute sum of delta h^* * delta h
                MultiFab::Add(dheightstarsum, dheight, 0, 0, 1, 0);
            }

            // write correlations to plotfile
            if (plot_int > 0 && istep%plot_int == 0)
            {
                // compute average of delta h^* * delta h
                MultiFab::Copy(dheightstaravg, dheightstarsum, 0, 0, 1, 0);
                dheightstaravg.mult( stats_count_inv, 0, 1, 0);

                const std::string& pltfile = amrex::Concatenate("star",istep,10);
                WriteSingleLevelPlotfile(pltfile, dheightstaravg, {"star"}, geom, time, 0);
            }

            if (do_fft_diag == 1) {

                // take ffts of each strip, add to running sum
                if (do_1d_x || do_1d_y) {

                    // copy distributed data into 1D data
                    height_onegrid.ParallelCopy(height, 0, 0, 1);

                    int numstrips = do_1d_x ? n_cells[1] : n_cells[0];

                    // work on each strip separately
                    for (int strip=0; strip<numstrips; ++strip) {

                        // copy strip into flat multifab
                        for ( MFIter mfi(height_flat_onegrid,TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

                            const Array4<Real> & h = height_onegrid.array(mfi);
                            const Array4<Real> & h_flat = height_flat_onegrid.array(mfi);

                            const Box& bx = mfi.tilebox();

                            if (do_1d_x) {
                                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                                {
                                    h_flat(i,j,k) = h(i,strip,k);
                                });
                            } else {
                                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                                {
                                    h_flat(i,j,k) = h(strip,i,k);  // convert column to row
                                });
                            }

                        }

#ifdef AMREX_USE_CUDA
                        using FFTplan = cufftHandle;
                        using FFTcomplex = cuDoubleComplex;
#else
                        using FFTplan = fftw_plan;
                        using FFTcomplex = fftw_complex;
#endif

                        Vector<std::unique_ptr<BaseFab<GpuComplex<Real> > > > spectral_field;
                        Vector<FFTplan> forward_plan;

                        Real sqrtnpts_flat = std::sqrt(npts_flat);

                        // take fft of strip and add magnitude of result to fft_sum
                        for (MFIter mfi(height_flat_onegrid); mfi.isValid(); ++mfi) {

                            // grab a single box
                            Box realspace_bx = mfi.fabbox();

                            // size of box
                            IntVect fft_size = realspace_bx.length(); // This will be different for FFTs of complex data

                            // this is the size of the box, except the 0th component is 'halved plus 1'
                            IntVect spectral_bx_size = fft_size;
                            spectral_bx_size[0] = fft_size[0]/2 + 1;

                            // spectral box
                            Box spectral_bx = Box(IntVect(0), spectral_bx_size - IntVect(1));

                            spectral_field.emplace_back(new BaseFab<GpuComplex<Real> >(spectral_bx,1,
                                                                                       The_Device_Arena()));
                            spectral_field.back()->setVal<RunOn::Device>(0.0); // touch the memory

                            FFTplan fplan;

#ifdef AMREX_USE_CUDA

                            cufftResult result = cufftPlan1d(&fplan, fft_size[0], CUFFT_D2Z, 1);
                            if (result != CUFFT_SUCCESS) {
                                AllPrint() << " cufftplan1d forward failed! Error: "
                                           << cufftErrorToString(result) << "\n";
                            }

#else // host

                            fplan = fftw_plan_dft_r2c_1d(fft_size[0],
                                                         height_flat_onegrid[mfi].dataPtr(),
                                                         reinterpret_cast<FFTcomplex*>
                                                         (spectral_field.back()->dataPtr()),
                                                         FFTW_ESTIMATE);

#endif

                            forward_plan.push_back(fplan);
                        }

                        ParallelDescriptor::Barrier();

                        // ForwardTransform
                        for (MFIter mfi(height_flat_onegrid); mfi.isValid(); ++mfi) {
                            int i = mfi.LocalIndex();
#ifdef AMREX_USE_CUDA
                            cufftSetStream(forward_plan[i], Gpu::gpuStream());
                            cufftResult result = cufftExecD2Z(forward_plan[i],
                                                              height_flat_onegrid[mfi].dataPtr(),
                                                              reinterpret_cast<FFTcomplex*>
                                                              (spectral_field[i]->dataPtr()));
                            if (result != CUFFT_SUCCESS) {
                                AllPrint() << " forward transform using cufftExec failed! Error: "
                                           << cufftErrorToString(result) << "\n";
                            }
#else
                            fftw_execute(forward_plan[i]);
#endif
                        }

                        // add magnitude to sum
                        for (MFIter mfi(fft_sum); mfi.isValid(); ++mfi) {

                            Array4< GpuComplex<Real> > spectral = (*spectral_field[0]).array();

                            Array4<Real> const& sum = fft_sum.array(mfi);

                            Box bx = mfi.fabbox();

                            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                            {
                                sum(i,j,k) += ( spectral(i,j,k).real()*spectral(i,j,k).real() + spectral(i,j,k).imag()*spectral(i,j,k).imag() ) / sqrtnpts_flat;
                                if (i == 0) sum(i,j,k) = 0.;
                            });
                        }

                        // destroy fft plan
                        for (int i = 0; i < forward_plan.size(); ++i) {
#ifdef AMREX_USE_CUDA
                            cufftDestroy(forward_plan[i]);
#else
                            fftw_destroy_plan(forward_plan[i]);
#endif
                        }
                    }

                    // write spectrum to plotfile
                    if (plot_int > 0 && istep%plot_int == 0)
                    {
                        // compute fft_avg by dividing fft_sum by stats_count*numstrips
                        MultiFab::Copy(fft_avg, fft_sum, 0, 0, 1, 0);

                        Real stats_count_fft_inv = 1./(stats_count*numstrips);
                        fft_avg.mult( stats_count_fft_inv, 0, 1, 0);

                        const std::string& pltfile = amrex::Concatenate("fft",istep,10);
                        WriteSingleLevelPlotfile(pltfile, fft_avg, {"fft"}, geom_fft, time, 0);
                    }

                } else {
                    Abort("Power spectrum not implemented for 2D");
                }

            } // end if test for fft diagnostic

        } // end loop over time steps

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