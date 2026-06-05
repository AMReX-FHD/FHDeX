#include "common_functions.H"
#include "FPU.H"
#include <AMReX_Vector.H>
#include "rng_functions.H"
#include "StructFact.H"
#include "chrono"

using namespace std::chrono;
using namespace amrex;

// argv contains the name of the inputs file entered at the command line
void main_driver(const char* argv)
{
    BL_PROFILE_VAR("main_driver()",main_driver);

    // store the current time so we can later compute total run time.
    Real strt_time = ParallelDescriptor::second();

    std::string inputs_file = argv;

    InitializeCommonNamespace();

    InitializeNamespace();

    int step_start, statsCount;
    amrex::Real time;

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

    // conserved quantaties
    MultiFab cu; // stretch, avg. mom, energy

    // staggered momentum
    MultiFab cumom;

    //statistics
    MultiFab cuMeans;
    MultiFab cuVars;

    MultiFab coVars;

    // Commenting out computations for higher moments    
//    MultiFab mom3;
//    MultiFab mom4;

    MultiFab cumomMeans;
    MultiFab cumomVars;

    // make BoxArray and Geometry
    BoxArray ba;
    Geometry geom;
    DistributionMapping dmap;

    IntVect dom_lo(AMREX_D_DECL(           0,            0,            0));
    IntVect dom_hi(AMREX_D_DECL(n_cells[0]-1, n_cells[1]-1, n_cells[2]-1));
    Box domain(dom_lo, dom_hi);

    // This defines the physical box, [-1,1] in each direction.
    RealBox real_box({AMREX_D_DECL(prob_lo[0],prob_lo[1],prob_lo[2])},
                     {AMREX_D_DECL(prob_hi[0],prob_hi[1],prob_hi[2])});

    // This defines a Geometry object
    Vector<int> is_periodic(AMREX_SPACEDIM,1);  // set to 1 (periodic by default)
    geom.define(domain,&real_box,CoordSys::cartesian,is_periodic.data());

    Real dt = fixed_dt;
    const Real* dx = geom.CellSize();
    const RealBox& realDomain = geom.ProbDomain();
    Real sys_volume = 1.0;
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        sys_volume *= (realDomain.hi(d) - realDomain.lo(d));
    }

    // MultiFabs to copy data into for snapshots for full 3D data
    MultiFab structFactConsMF;

    BoxArray ba_pencil;
    DistributionMapping dmap_pencil;

    // Structure factor for pencils
    // enabled if do_1D=1
    Vector < StructFact > structFactConsArray;
    int cnt = 0;
    int numvars;
    std::string x;

    // "conserved" variable structure factor will contain
    // stretch
    // mom (averaged)
    // energy
    int structVarsCons = 3;

    Vector< std::string > cons_var_names;
    cons_var_names.resize(structVarsCons);

    cnt = 0;

    cons_var_names[cnt] = "stretch";
    ++cnt;
    
    cons_var_names[cnt] = "mom";
    ++cnt;
    
    cons_var_names[cnt] = "energy";
    ++cnt;

    // scale SF results by inverse cell volume
    Vector<Real> var_scaling_cons;
    var_scaling_cons.resize(structVarsCons*(structVarsCons+1)/2);
    for (int d=0; d<var_scaling_cons.size(); ++d) {
        var_scaling_cons[d] = 1./(dx[0]*dx[1]*dx[2]);
    }

    /////////////////////////////////////////////
    // Initialize based on fresh start or restart
    /////////////////////////////////////////////
    if (restart > 0) {

        ReadCheckPoint(step_start, time, statsCount, geom, domain, cu, cuMeans, cuVars,
                       cumom, cumomMeans, cumomVars,
                       coVars, ba, dmap); // REDEFINE

        if (reset_stats == 1) statsCount = 1;

    } else {

        ///////////////////////////////////////////
        // Define geometry, box arrays and MFs
        ///////////////////////////////////////////

        // Initialize the boxarray "ba" from the single box "bx"
        ba.define(domain);

        // Break up boxarray "ba" into chunks no larger than "max_grid_size" along a direction
        // note we are converting "Vector<int> max_grid_size" to an IntVect
        ba.maxSize(IntVect(max_grid_size));

        // how boxes are distrubuted among MPI processes
        dmap.define(ba);

        // conserved quantaties
        // 0        (stretch)
        // 1        (momentum, avg. to cc)
        // 2        (energy)
        cu.define(ba,dmap,3,ngc);

        cumom.define(convert(ba,nodal_flag_dir[0]), dmap, 1, ngc);

        cuMeans.define(ba,dmap,3,ngc);
        cuVars.define(ba,dmap,3,ngc);
        cuMeans.setVal(0.0);
        cuVars.setVal(0.0);

//        if (plot_mom3) {
//            mom3.define(ba,dmap,nvars+1,0); // nvars (vel instead of momemtum) with temperature and pressure
//            mom3.setVal(0.0);
//        }
//
//        if (plot_mom4) {
//            mom4.define(ba,dmap,nvars+1,0); // nvars (vel instead of momemtum) with temperature and pressure
//            mom4.setVal(0.0);
//        }
//
        // List of covariances (all cell centered)
        // 0: <stretch mom>
        // 1: <stretch energy>
        // 2: <mom energy>
        coVars.define(ba,dmap,3,0);
        coVars.setVal(0.0);

        cumomMeans.define(convert(ba,nodal_flag_dir[0]), dmap, 1, 0);
        cumomVars.define(convert(ba,nodal_flag_dir[0]), dmap, 1, 0);
        cumomMeans.setVal(0.);
        cumomVars.setVal(0.);

        ///////////////////////////////////////////
        // Initialize everything
        ///////////////////////////////////////////

        // initialize conserved variables
        //InitConsVarStag(cu,cumom,geom); // REDEFINE

        // Set BC: 1) fill boundary 2) physical
        cu.FillBoundary(geom.periodicity());
        cumom.FillBoundary(geom.periodicity());

        if (plot_int > 0) {
            WritePlotFile(0, 0.0, geom, cu, cuMeans, cuVars, cumom, cumomMeans, cumomVars, coVars); // REDEFINE
        }

        step_start = 1;
        time = 0.;
        statsCount = 1;

    } // else restart/non-restart

    ///////////////////////////////////////////
    // Setup 1D Structure factor
    ///////////////////////////////////////////

    if (struct_fact_int > 0) {

        structFactConsMF.define(ba,dmap,structVarsCons,0);

        MultiFab pencil;

        // we are only calling ExtractXPencil here to obtain
        // a built version of pencil so can obtain what we need to build the
        // structure factor objects for pencil data
        ExtractXPencil(cu, pencil, 0, 0, 0, 1);

        ba_pencil = pencil.boxArray();
        dmap_pencil = pencil.DistributionMap();

        structFactConsArray.resize(n_cells[1]*n_cells[2]);

        for (int i=0; i<n_cells[1]*n_cells[2];  ++i) {
            structFactConsArray[i].define(ba_pencil,dmap_pencil,cons_var_names,var_scaling_cons);
        }

    }

    //fluxes (except momentum) at faces
    // need +4 to separate out heat, viscous heating (diagonal vs shear)  and Dufour contributions to the energy flux
    // stacked at the end (see below)
    // index: flux term
    // 0: flux for stretch: momentum
    // 1: flux for energy
    MultiFab faceflux;
    faceflux.define(convert(ba,nodal_flag_x), dmap, 2, 0);

    //momentum flux (center)
    MultiFab cenflux;
    cenflux.define(ba,dmap,1,1);

    /////////////////////////////////////////////////
    //Time stepping loop
    /////////////////////////////////////////////////

    for (int step=step_start;step<=max_step;++step) {

        // timer
        Real ts1 = ParallelDescriptor::second();
        
        // time step
        RK3step(cu, cumom, faceflux, cenflux, geom, dt, step); // REDEFINE

        // timer
        Real ts2 = ParallelDescriptor::second() - ts1;
        ParallelDescriptor::ReduceRealMax(ts2, ParallelDescriptor::IOProcessorNumber());
        if (step%100 == 0) {
            amrex::Print() << "Advanced step " << step << " in " << ts2 << " seconds\n";
        }

        // timer
        Real aux1 = ParallelDescriptor::second();

        // reset statistics after n_steps_skip
        // if n_steps_skip is negative, we use it as an interval
        if ((n_steps_skip > 0 && step == n_steps_skip) ||
            (n_steps_skip < 0 && step%amrex::Math::abs(n_steps_skip) == 0) ) {

            cuMeans.setVal(0.0);
            cuVars.setVal(0.0);

            cumomMeans.setVal(0.0);
            cumomVars.setVal(0.0);

            coVars.setVal(0.0);

//            if (plot_mom3) mom3.setVal(0.0);

//            if (plot_mom4) mom4.setVal(0.0);

            std::printf("Resetting stat collection.\n");

            statsCount = 1;
        }

        // Evaluate Statistics
        evaluateStats(cu, cuMeans, cuVars,
                            cumom, cumomMeans, cumomVars, coVars,
                            statsCount, geom); // REDEFINE
        statsCount++;
        if (step%100 == 0) {
            amrex::Print() << "Mean stretch: "  << ComputeSpatialMean(cu, 0)
                           << " Mean momentum:" << ComputeSpatialMean(cumom, 0)
                           << " Mean energy:"   << ComputeSpatialMean(cu, 1)
                           << "\n";
        }


        // write a plotfile
        bool writePlt = false;
        if (plot_int > 0) {
            if (n_steps_skip >= 0) { // for positive n_steps_skip, write out at plot_int
                writePlt = (step%plot_int == 0);
            }
            else if (n_steps_skip < 0) { // for negative n_steps_skip, write out at plot_int-1
                writePlt = ((step+1)%plot_int == 0);
            }
        }

        if (writePlt) {
            WritePlotFile(step, time, geom, cu, cuMeans, cuVars, cumom, cumomMeans, cumomVars, coVars); // REDEFINE
        }

        bool SF_snapshot_taken = false;

        // collect a snapshot for structure factor
        if (struct_fact_int > 0 &&
            step > amrex::Math::abs(n_steps_skip) &&
            step%struct_fact_int == 0) {

            SF_snapshot_taken = true;

            ////////////// structFactConsMF /////////////
            cnt = 0;

            // copy [stretch, mom, energy]
            numvars = 3;
            MultiFab::Copy(structFactConsMF, cu, 0, cnt, numvars, 0);
            cnt+=numvars;

            // copy momFace
            ShiftFaceToCC(cumom,0,structFactConsMF,cnt,1);
            ++cnt;
            ////////////////////////////////////////////////////


            for (int i=0; i<n_cells[1]*n_cells[2]; ++i) {
                {
                    MultiFab pencil;
                    ExtractXPencil(structFactConsMF, pencil, i/n_cells[1], i%n_cells[1], 
                                   0, structVarsCons);
                    structFactConsArray[i].FortStructure(pencil);
                }
            }

        } // logic for doing structure factor

        // write out structure factor
        if (struct_fact_int > 0 &&
            SF_snapshot_taken &&
            plot_int > 0 &&
            step%plot_int == 0) {

            MultiFab cons_mag, cons_realimag;

            cons_mag     .define(ba_pencil,dmap_pencil,  structFactConsArray[0].get_ncov(),0);
            cons_realimag.define(ba_pencil,dmap_pencil,2*structFactConsArray[0].get_ncov(),0);

            cons_mag.setVal(0.0);
            cons_realimag.setVal(0.0);

            for (int i=0; i<n_cells[1]*n_cells[2]; ++i) {
                structFactConsArray[i].AddToExternal(cons_mag,cons_realimag);
            }

            Real ncellsinv = 1.0/(n_cells[1]*n_cells[2]);
            cons_mag.mult(ncellsinv);
            cons_realimag.mult(ncellsinv);

            WritePlotFilesSF_1D(cons_mag,cons_realimag,step,time,
                                structFactConsArray[0].get_names(),"plt_SF_1D");

        }

        // write checkpoint file
        if (chk_int > 0 && step > 0 && step%chk_int == 0)
        {
            WriteCheckPoint(step, time, statsCount, geom, cu, cuMeans, cuVars,
                            cumom, cumomMeans, cumomVars,
                            coVars); // REDEFINE
        }

        // timer
        Real aux2 = ParallelDescriptor::second() - aux1;
        ParallelDescriptor::ReduceRealMax(aux2,  ParallelDescriptor::IOProcessorNumber());
        if (step%100 == 0) {
            amrex::Print() << "Aux time (stats, struct fac, plotfiles) " << aux2 << " seconds\n";
        }

        time = time + dt;

        // MultiFab memory usage
        const int IOProc = ParallelDescriptor::IOProcessorNumber();

        amrex::Long min_fab_megabytes  = amrex::TotalBytesAllocatedInFabsHWM()/1048576;
        amrex::Long max_fab_megabytes  = min_fab_megabytes;

        ParallelDescriptor::ReduceLongMin(min_fab_megabytes, IOProc);
        ParallelDescriptor::ReduceLongMax(max_fab_megabytes, IOProc);

        if (step%100 == 0) {
            amrex::Print() << "High-water FAB megabyte spread across MPI nodes: ["
                           << min_fab_megabytes << " ... " << max_fab_megabytes << "]\n";
        }

        min_fab_megabytes  = amrex::TotalBytesAllocatedInFabs()/1048576;
        max_fab_megabytes  = min_fab_megabytes;

        ParallelDescriptor::ReduceLongMin(min_fab_megabytes, IOProc);
        ParallelDescriptor::ReduceLongMax(max_fab_megabytes, IOProc);

        if (step%100 == 0) {
            amrex::Print() << "Curent     FAB megabyte spread across MPI nodes: ["
                           << min_fab_megabytes << " ... " << max_fab_megabytes << "]\n";
        }
    }

    // timer
    Real stop_time = ParallelDescriptor::second() - strt_time;
    ParallelDescriptor::ReduceRealMax(stop_time, ParallelDescriptor::IOProcessorNumber());
    amrex::Print() << "Run time = " << stop_time << std::endl;
}
