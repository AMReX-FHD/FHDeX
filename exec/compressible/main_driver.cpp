
//#include "hydro_test_functions.H"
//#include "hydro_test_functions_F.H"

//#include "hydro_functions.H"
//#include "hydro_functions_F.H"
//#include "StochMFlux.H"
//#include "StructFact.H"

#include "rng_functions_F.H"

#include "common_functions.H"
#include "common_functions_F.H"

#include "gmres_functions.H"
#include "gmres_functions_F.H"

#include "common_namespace.H"
#include "common_namespace_declarations.H"

#include "gmres_namespace.H"
#include "gmres_namespace_declarations.H"

#include "compressible_functions.H"
#include "compressible_functions_F.H"

#include "exec_functions.H"

#include <AMReX_VisMF.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_MultiFab.H>

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

    //if gas heat capacities are negative, calculate using dofs. This will only update the Fortran values.
    get_hc_gas();
  
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
    const Real* dx = geom.CellSize();

    // how boxes are distrubuted among MPI processes
    DistributionMapping dmap(ba);

    /////////////////////////////////////////
    //Initialise rngs
    /////////////////////////////////////////
    const int n_rngs = 1;

    int fhdSeed = 1;
    int particleSeed = 2;
    int selectorSeed = 3;
    int thetaSeed = 4;
    int phiSeed = 5;
    int generalSeed = 6;

    const int proc = ParallelDescriptor::MyProc();

    fhdSeed += 10000*proc;
    particleSeed += 20000*proc;
    selectorSeed += 30000*proc;
    thetaSeed += 40000*proc;
    phiSeed += 50000*proc;
    generalSeed += 60000*proc;

    //Initialise rngs
    rng_initialize(&fhdSeed,&particleSeed,&selectorSeed,&thetaSeed,&phiSeed,&generalSeed);
    /////////////////////////////////////////

    Real Cfl_max[4];
    int eval_dt = 0;

    //transport properties
    MultiFab eta  (ba,dmap,1,ngc);
    MultiFab zeta (ba,dmap,1,ngc);
    MultiFab kappa(ba,dmap,1,ngc);
    MultiFab chi(ba,dmap,1,ngc);
    MultiFab D(ba,dmap,1,ngc);

    //conserved quantaties
    MultiFab cu  (ba,dmap,nvars,ngc);
    MultiFab cup (ba,dmap,nvars,ngc);
    MultiFab cup2(ba,dmap,nvars,ngc);
    MultiFab cup3(ba,dmap,nvars,ngc);

    MultiFab cuMeans(ba,dmap,nvars,ngc);
    MultiFab cuVars(ba,dmap,nvars,ngc);

    MultiFab cuMeansAv(ba,dmap,nvars,ngc);
    MultiFab cuVarsAv(ba,dmap,nvars,ngc);

    //primative quantaties
    MultiFab prim(ba,dmap,nprimvars,ngc);

    MultiFab primMeans(ba,dmap,nprimvars,ngc);
    MultiFab primVars(ba,dmap,nprimvars + 5,ngc);

    MultiFab primMeansAv(ba,dmap,nprimvars,ngc);
    MultiFab primVarsAv(ba,dmap,nprimvars + 5,ngc);

    MultiFab spatialCross(ba,dmap,3,ngc);

    MultiFab spatialCrossAv(ba,dmap,3,ngc);

    cuMeans.setVal(0.0);
    cuVars.setVal(0.0);

    primMeans.setVal(0.0);
    primVars.setVal(0.0);

    spatialCross.setVal(0.0);

    //possibly for later
    MultiFab source(ba,dmap,nprimvars,ngc);

    source.setVal(0.0);

    //Initialize physical parameters from input vals

    prim.setVal(rho0,0,1,ngc);
    prim.setVal(0,1,1,ngc);
    prim.setVal(0,2,1,ngc);
    prim.setVal(0,3,1,ngc);
    prim.setVal(T_init[0],4,1,ngc);

    double massvec[nspecies];
    double intEnergy, T0;

    T0 = T_init[0];

    for(int i=0;i<nspecies;i++)
    {
        prim.setVal(rhobar[i],6+i,1,ngc);
        cu.setVal(rhobar[i],5+i,1,ngc);

        massvec[i] = rhobar[i]*rho0;
    }

    get_energy(&intEnergy, massvec, &T0);

    cu.setVal(rho0,0,1,ngc);
    cu.setVal(0,1,1,ngc);
    cu.setVal(0,2,1,ngc);
    cu.setVal(0,3,1,ngc);
    cu.setVal(intEnergy,4,1,ngc);

    //Print() << intEnergy << "\n";

    //while(true);

    //fluxes
    std::array< MultiFab, AMREX_SPACEDIM > flux;
    AMREX_D_TERM(flux[0].define(convert(ba,nodal_flag_x), dmap, nvars, 0);,
                 flux[1].define(convert(ba,nodal_flag_y), dmap, nvars, 0);,
                 flux[2].define(convert(ba,nodal_flag_z), dmap, nvars, 0););

    //stochastic fluxes
    std::array< MultiFab, AMREX_SPACEDIM > stochFlux;
    AMREX_D_TERM(stochFlux[0].define(convert(ba,nodal_flag_x), dmap, nvars, 0);,
                 stochFlux[1].define(convert(ba,nodal_flag_y), dmap, nvars, 0);,
                 stochFlux[2].define(convert(ba,nodal_flag_z), dmap, nvars, 0););

    AMREX_D_TERM(stochFlux[0].setVal(0.0);,
                 stochFlux[1].setVal(0.0);,
                 stochFlux[2].setVal(0.0););

    Real time = 0;

    int step, statsCount;


    //Initialise everything
    calculateTransportCoeffs(prim, eta, zeta, kappa);

    conservedToPrimitive(prim, cu);
    cu.FillBoundary(geom.periodicity());
    prim.FillBoundary(geom.periodicity());

    calculateTransportCoeffs(prim, eta, zeta, kappa);

    eta.FillBoundary(geom.periodicity());
    zeta.FillBoundary(geom.periodicity());
    kappa.FillBoundary(geom.periodicity());

    setBC(prim, cu, eta, zeta, kappa);

    calculateFlux(cu, prim, eta, zeta, kappa, flux, stochFlux, geom, dx, dt);
    statsCount = 1;

    //double* test;

    //test = new double[nvars];

    //test[10] = 1.0;    

    //Time stepping loop
    for(step=1;step<=max_step;++step)
    {

        RK3step(cu, cup, cup2, cup3, prim, source, eta, zeta, kappa, flux, stochFlux, geom, dx, dt);

        if(step == n_steps_skip)
        {
            cuMeans.setVal(0.0);
            cuVars.setVal(0.0);

            primMeans.setVal(0.0);
            primVars.setVal(0.0);

            spatialCross.setVal(0.0);

            statsCount = 1;

            dt = 2.0*dt;

        }

        evaluateStats(cu, cuMeans, cuVars, prim, primMeans, primVars, spatialCross, statsCount, dx);

        statsCount++;

        if(step%5000 == 0)
        {    
                amrex::Print() << "Advanced step " << step << "\n";
        }

        if (plot_int > 0 && step > 0 && step%plot_int == 0)
        {

           yzAverage(cuMeans, cuVars, primMeans, primVars, spatialCross, cuMeansAv, cuVarsAv, primMeansAv, primVarsAv, spatialCrossAv);

           WritePlotFile(step, time, geom, cu, cuMeansAv, cuVarsAv, prim, primMeansAv, primVarsAv, spatialCross);
        }

        time = time + dt;
    }

    Real stop_time = ParallelDescriptor::second() - strt_time;
    ParallelDescriptor::ReduceRealMax(stop_time);
    amrex::Print() << "Run time = " << stop_time << std::endl;

    amrex::Finalize();

}

