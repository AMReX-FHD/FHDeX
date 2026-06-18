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
    int nskip, nstat, ncopies, icor, istat;
    int alg_type, avg_type;
    int seed;
    amrex::Real phileft, phiright;
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

        nskip = 0;
        pp.query("nskip",nskip);

        nstat = -1;
        pp.query("nstat",nstat);

        icor = n_cell/2;;
        pp.query("icor",icor);

        ncopies = 1;
        pp.query("ncopies",ncopies);

        istat = 0;

        npts_scale = 1.;
        pp.query ("npts_scale",npts_scale);

        alg_type = 0;
        pp.query ("alg_type",alg_type);

        avg_type = 0;
        pp.query ("avg_type",avg_type);

        phileft = 32;
        pp.query ("phileft",phileft);

        phiright = 32;
        pp.query ("phiright",phiright);

        seed = 0;
        pp.query("seed",seed);

        cfl=.9;
        pp.query ("cfl",cfl);


        // By default, the boundary conditions will be set to periodic, or bc_lo = bc_hi = 0.
        //Other options in this program include bc_lo, bc_hi = 2 for homogeneous Neumann, or
        //bc_lo, bc_hi = 3 for external Dirichlet boundary conditions.
        pp.queryarr("bc_lo", bc_lo);
        pp.queryarr("bc_hi", bc_hi);
    }

    amrex::Real copies = ncopies;

    Vector<int> is_periodic(AMREX_SPACEDIM,0);
    for (int idim=0; idim < AMREX_SPACEDIM; ++idim) {
        if (bc_lo[idim] == BCType::int_dir && bc_hi[idim] == BCType::int_dir){
            is_periodic[idim] = 1;
        }
    }

    // make BoxArray and Geometry
    BoxArray ba;
    Geometry geom;
    MultiFab stats_onegrid;
    {
        IntVect dom_lo(AMREX_D_DECL(       0,        0,        0));
        IntVect dom_hi(AMREX_D_DECL(n_cell-1, ncopies-1, n_cell-1));
        Box domain(dom_lo, dom_hi);

        // Initialize the boxarray "ba" from the single box "bx"
        ba.define(domain);
        // Break up boxarray "ba" into chunks no larger than "max_grid_size" along a direction
        ba.maxSize(max_grid_size);


        // This defines the physical box, [0,1] in each direction.
        RealBox real_box({AMREX_D_DECL( 0.0, 0.0, 0.0)},
                         {AMREX_D_DECL( 1.0, copies, 1.0)});

        // This defines a Geometry object
        geom.define(domain,&real_box,CoordSys::cartesian,is_periodic.data());

        BoxArray ba_onegrid;
        Vector<int> pmap(1);
        const int ioproc = ParallelDescriptor::IOProcessorNumber();
        if(ioproc != 0){
            amrex::Print() << "IO processor is " << ioproc << std::endl;
        }
        pmap[0] = ioproc;
        ba_onegrid.define(domain);
        //DistributionMapping dmap_onegrid(ba_onegrid,pmap);
        DistributionMapping dmap_onegrid(pmap);
        stats_onegrid.define(ba_onegrid,dmap_onegrid,4,0);

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

    int restart = -1;
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

    MultiFab stats(ba, dm, 4, 0);
    stats.setVal(0.);


    GpuArray<Real,AMREX_SPACEDIM> dx = geom.CellSizeArray();

    init_phi(phi_new, geom,npts_scale, Ncomp, phileft, phiright);

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

//    Real coeff = AMREX_D_TERM(   2./(dx[0]*dx[0]),
//                               + 2./(dx[1]*dx[1]),
//                               + 2./(dx[2]*dx[2]) );
    Real coeff = 2./(dx[0]*dx[0]);
    Real dt = cfl/(2.0*coeff);

    amrex::Print() << "dt = " << dt << " dx = " << dx[0] << std::endl;

    // time = starting time in the simulation
    Real time = 0.0;

    // Write a plotfile of the initial data if plot_int > 0 (plot_int was defined in the inputs file)
    if (plot_int > 0)
    {
        int n = 0;
        const std::string& pltfile = amrex::Concatenate("plt",n,7);
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

    flux[1].setVal(0.);
    AMREX_D_TERM(stochFlux[0].setVal(0.0);,
                 stochFlux[1].setVal(0.0);,
                 stochFlux[2].setVal(0.0););

    nsteps += nskip;

    for (int n = 1; n <= nsteps; ++n)
    {
        MultiFab::Copy(phi_old, phi_new, 0, 0, Ncomp, 0);

        // new_phi = old_phi + dt * (something)
        advance(phi_old, phi_new, flux, stochFlux, dt, npts_scale, geom, bc, Ncomp, phileft, phiright, avg_type);
        time = time + dt;

        // Tell the I/O Processor to write out which step we're doing
        amrex::Print() << "Advanced step " << n << "\n";

        // Write a plotfile of the current data (plot_int was defined in the inputs file)
        if (nstat > 0 && n > nskip && n%nstat == 0){

            istat += 1;

            for ( MFIter mfi(stats); mfi.isValid(); ++mfi )
            {
                const Box& vbx = mfi.validbox();
                const auto lo = amrex::lbound(vbx);
                const auto hi = amrex::ubound(vbx);

                auto const& phiNew = phi_new.array(mfi);
                auto const& stat_arr = stats.array(mfi);

                for (auto k = lo.z; k <= hi.z; ++k) {
                    for (auto j = lo.y; j <= hi.y; ++j) {
                        for (auto i = lo.x; i <= hi.x; ++i) {
                            stat_arr(i,j,k,0) += phiNew(i,j,k);
                            stat_arr(i,j,k,1) += phiNew(i,j,k)*phiNew(i,j,k);
                            stat_arr(i,j,k,2) += phiNew(i,j,k)*phiNew(icor,j,k);
                            stat_arr(i,j,k,3) += phiNew(icor,j,k);
                        }
                    }
                }

            }

        }
        if (plot_int > 0 && n%plot_int == 0)
        {
            const std::string& pltfile = amrex::Concatenate("plt",n,7);
            WriteSingleLevelPlotfile(pltfile, phi_new, {"phi"}, geom, time, n);

            if(nstat > 0 && istat > 0){

                stats_onegrid.ParallelCopy(stats,0,0,4);
                amrex::Real mean[n_cell];
                amrex::Real var[n_cell];
                amrex::Real cor[n_cell];
                amrex::Real mean_avg[n_cell];
                amrex::Real var_avg[n_cell];
                amrex::Real cor_avg[n_cell];

                amrex::Real num_stats = istat;

                for (MFIter mfi(stats_onegrid); mfi.isValid(); ++mfi ) {

                    const Array4<Real>& stat_arr = stats_onegrid.array(mfi);

                    for (int i=0; i<n_cell; i++){
                        mean_avg[i]=0.;
                        var_avg[i]=0.;
                        cor_avg[i]=0.;
                    }

                    amrex::Real sum_check;

                    for (int j=0; j<ncopies; j++){

                        for (int i=0; i<n_cell; i++){

                            mean[i] = stat_arr(i,j,0,0)/num_stats;
                            mean_avg[i] += mean[i];
                            var[i] = stat_arr(i,j,0,1)/num_stats - mean[i]*mean[i];
                            var_avg[i] += var[i];
                            cor[i] = stat_arr(i,j,0,2)/num_stats - mean[i]*stat_arr(icor,j,0,3)/num_stats;
                            cor_avg[i] += cor[i];

                        }

                    }
                    for (int i=0; i<n_cell; i++){
                        mean_avg[i]=mean_avg[i]/copies;
                        var_avg[i]=var_avg[i]/copies;
                        cor_avg[i]=cor_avg[i]/copies;
                    }

                }

                if (ParallelDescriptor::IOProcessor()) {

                    std::ofstream meandat;
                    std::string meanBaseName = "mean";
                    std::string meanName = Concatenate(meanBaseName,n,9);
                    meanName += ".dat";
                    meandat.open(meanName);

                    std::ofstream vardat;
                    std::string varBaseName = "var";
                    std::string varName = Concatenate(varBaseName,n,9);
                    varName += ".dat";
                    vardat.open(varName);

                    std::ofstream cordat;
                    std::string corBaseName = "cor";
                    std::string corName = Concatenate(corBaseName,n,9);
                    corName += ".dat";
                    cordat.open(corName);

                    for (int i=0; i<n_cell; i++) {
                        meandat << (i+0.5)*dx[0]  << " " << mean_avg[i] << std::endl;
                        vardat << (i+0.5)*dx[0]  << " " << var_avg[i] << std::endl;
                        cordat << (i+0.5)*dx[0]  << " " << cor_avg[i] << std::endl;
                    }

                    meandat.close();
                    vardat.close();
                    cordat.close();

                }
            }
        }

//        amrex::Real Ephi=0.;
//        amrex::Real Ephi2=0.;
/*
        Vector<Real> Ephi(2,0.);
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
                    }
                }
            }

        }

        const int IOProc = ParallelDescriptor::IOProcessorNumber();
        ParallelDescriptor::ReduceRealSum(Ephi.dataPtr(),2);
        ParallelDescriptor::ReduceRealMin(Ephimin);
        amrex::Real scale = n_cell*n_cell;
        amrex::Real scale2 =  AMREX_D_TERM( dx[0],
                               * dx[1],
                               * dx[2] );

        amrex::Print() << "phi variance = " << Ephi[1]/scale - (Ephi[0]*Ephi[0]
                /(scale*scale)) << std::endl;
        amrex::Print() << "phi integral = " << Ephi[0]*scale2 << std::endl;
        amrex::Print() << "phi min = " << Ephimin << std::endl;
*/

    }

    // Call the timer again and compute the maximum difference between the start time and stop time
    //   over all processors
    auto stop_time = ParallelDescriptor::second() - strt_time;
    const int IOProc = ParallelDescriptor::IOProcessorNumber();
    ParallelDescriptor::ReduceRealMax(stop_time,IOProc);

    // Tell the I/O Processor to write out the "run time"
    amrex::Print() << "Run time = " << stop_time << std::endl;
}