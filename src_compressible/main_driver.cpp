#include "common_functions.H"
#include "compressible_functions.H"


#include "rng_functions.H"

#include "StructFact.H"

#include "chemistry_functions.H"

#include "MFsurfchem_functions.H"

#include "chrono"

#ifdef MUI
#include "surfchem_mui_functions.H"
#endif

using namespace std::chrono;
using namespace amrex;

// argv contains the name of the inputs file entered at the command line
void main_driver(const char* argv)
{
    BL_PROFILE_VAR("main_driver()",main_driver);

    // store the current time so we can later compute total run time.
    Real strt_time = ParallelDescriptor::second();

    std::string inputs_file = argv;

    // read the inputs file for common
    InitializeCommonNamespace();

    // read the inputs file for compressible
    InitializeCompressibleNamespace();

    // read the inputs file for chemistry
    InitializeChemistryNamespace();

    // read the inputs file for MFsurfchem
    InitializeMFSurfchemNamespace();

#ifdef MUI
    // read the inputs file for surfchem_mui
    InitializeSurfChemMUINamespace();

    if (n_ads_spec>0) {
        Abort("MFsurfchem cannot be used in compressible_mui");
    }
#endif

    if (nvars != AMREX_SPACEDIM + 2 + nspecies) {
        Abort("nvars must be equal to AMREX_SPACEDIM + 2 + nspecies");
    }

    if (nprimvars != AMREX_SPACEDIM + 3 + 2*nspecies) {
        Abort("nprimvars must be equal to AMREX_SPACEDIM + 3 + 2*nspecies");
    }

    // if gas heat capacities in the namelist are negative, calculate them using using dofs.
    // This will only update the Fortran values.
    GetHcGas();

    // check bc_vel_lo/hi to determine the periodicity
    Vector<int> is_periodic(AMREX_SPACEDIM,0);  // set to 0 (not periodic) by default
    for (int i=0; i<AMREX_SPACEDIM; ++i) {
        // lo/hi consistency check
        if (bc_vel_lo[i] == -1 || bc_vel_hi[i] == -1) {
            if (bc_vel_lo[i] != bc_vel_hi[i]) {
                Abort("Inconsistent periodicity definition in bc_vel_lo/hi");
            }
            else {
                is_periodic[i] = 1;
            }
        }
    }

    // for each direction, if bc_vel_lo/hi is periodic, then
    // set the corresponding bc_mass_lo/hi and bc_therm_lo/hi to periodic
    SetupBC();

    // if multispecies
    if (algorithm_type == 2) {
        // compute wall concentrations if BCs call for it
        SetupCWall();
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
    const Real* dx = geom.CellSize();
    const RealBox& realDomain = geom.ProbDomain();

    // how boxes are distrubuted among MPI processes
    DistributionMapping dmap(ba);

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

    /////////////////////////////////////////

    // transport properties
    /*
      Referring to K. Balakrishnan et al.,
      "Fluctuating hydrodynamics of multispecies nonreactive mixtures",
      Phys. Rev. E, 89, 1, 2014

      "kappa" and "zeta" in the code have opposite meanings from what they
      represent in the paper.  So kappa in the paper is bulk viscosity (see
      the equation for Pi immediately after (28)), but in the code it's zeta.
      Zeta is a thermodiffusion coefficient (see the equation for Q'
      immediately before (25)), but in the code it's kappa... and furthermore
      I believe kappa in the code is actually zeta/T^2.
    */
    MultiFab eta  (ba,dmap,1,ngc);
    MultiFab zeta (ba,dmap,1,ngc);
    MultiFab kappa(ba,dmap,1,ngc);
    MultiFab chi(ba,dmap,nspecies,ngc);
    MultiFab D(ba,dmap,nspecies*nspecies,ngc);

    eta.setVal(1.0,0,1,ngc);
    zeta.setVal(1.0,0,1,ngc);
    kappa.setVal(1.0,0,1,ngc);
    chi.setVal(1.0,0,nspecies,ngc);
    D.setVal(1.0,0,nspecies*nspecies,ngc);

    // conserved quantaties
    // in C++ indexing (add +1 for F90)
    // 0        (rho;     density)
    // 1-3      (j;       momentum)
    // 4        (rho*E;   total energy)
    // 5:5+ns-1 (rho*Yk;  mass densities)
    MultiFab cu  (ba,dmap,nvars,ngc);
    MultiFab cup (ba,dmap,nvars,ngc);
    MultiFab cup2(ba,dmap,nvars,ngc);
    MultiFab cup3(ba,dmap,nvars,ngc);

    //primative quantaties
    // in C++ indexing (add +1 for F90)
    // 0            (rho; density)
    // 1-3          (vel; velocity)
    // 4            (T;   temperature)
    // 5            (p;   pressure)
    // 6:6+ns-1     (Yk;  mass fractions)
    // 6+ns:6+2ns-1 (Xk;  mole fractions)
    MultiFab prim(ba,dmap,nprimvars,ngc);

    MultiFab surfcov;
    MultiFab dNadsdes;
    MultiFab dNads;
    MultiFab dNdes;
    if (n_ads_spec>0) {
        surfcov.define(ba,dmap,n_ads_spec,ngc);
        dNadsdes.define(ba,dmap,n_ads_spec,ngc);
        dNads.define(ba,dmap,n_ads_spec,ngc);
        dNdes.define(ba,dmap,n_ads_spec,ngc);
    }

    //statistics
    MultiFab cuMeans  (ba,dmap,nvars,ngc);
    MultiFab cuVars   (ba,dmap,nvars,ngc);
    MultiFab cuMeansAv(ba,dmap,nvars,ngc);
    MultiFab cuVarsAv (ba,dmap,nvars,ngc);

    cuMeans.setVal(0.0);
    cuVars.setVal(0.0);

    MultiFab primMeans  (ba,dmap,nprimvars  ,ngc);
    MultiFab primVars   (ba,dmap,nprimvars+5,ngc);
    MultiFab primMeansAv(ba,dmap,nprimvars  ,ngc);
    MultiFab primVarsAv (ba,dmap,nprimvars+5,ngc);
    primMeans.setVal(0.0);
    primVars.setVal(0.0);

    //miscStats & miscVals (only used internally -- not outputted)
    MultiFab miscStats(ba,dmap,10,ngc);
    Real miscVals[20];

    // spatial cross correlations
    // 0: slice average of mean temperature
    // 1: yz- average of mean temperature
    // 2: <T*T>
    // 3: <delta T* delta T>
    // 4: <delta T* delta rho>
    // 5: <delta u* delta rho>
    // 6: <delta jx* delta rho>
    // 7: <delta rho* delta rhoE>
    MultiFab spatialCross(ba,dmap,8,ngc);
    MultiFab spatialCrossAv(ba,dmap,8,ngc);

    miscStats.setVal(0.0);
    spatialCross.setVal(0.0);
    spatialCrossAv.setVal(0.0);

    // external source term - currently only chemistry source considered for nreaction>0
    MultiFab source(ba,dmap,nvars,ngc);
    source.setVal(0.0);

    MultiFab ranchem;
    if (nreaction>0) ranchem.define(ba,dmap,nreaction,ngc);

    //fluxes
    // need +4 to separate out heat, viscous heating (diagonal vs shear)  and Dufour contributions to the energy flux
    // stacked at the end (see below)
    // index: flux term
    // 0: density
    // 1: x-momentum
    // 2: y-momentum
    // 3: z-momentum
    // 4: total energy
    // 5:nvars-1: species flux (nvars = nspecies+5)
    // nvars: heat flux
    // nvars + 1: viscous heating (diagonal)
    // nvars + 2: viscous heating (shear)
    // nvars + 3: Dufour effect
    std::array< MultiFab, AMREX_SPACEDIM > flux;
    AMREX_D_TERM(flux[0].define(convert(ba,nodal_flag_x), dmap, nvars+4, 0);,
                 flux[1].define(convert(ba,nodal_flag_y), dmap, nvars+4, 0);,
                 flux[2].define(convert(ba,nodal_flag_z), dmap, nvars+4, 0););

    //stochastic fluxes
    std::array< MultiFab, AMREX_SPACEDIM > stochFlux;
    AMREX_D_TERM(stochFlux[0].define(convert(ba,nodal_flag_x), dmap, nvars, 0);,
                 stochFlux[1].define(convert(ba,nodal_flag_y), dmap, nvars, 0);,
                 stochFlux[2].define(convert(ba,nodal_flag_z), dmap, nvars, 0););

    AMREX_D_TERM(stochFlux[0].setVal(0.0);,
                 stochFlux[1].setVal(0.0);,
                 stochFlux[2].setVal(0.0););

    MultiFab rancorn;
    rancorn.define(convert(ba,nodal_flag), dmap, 1, 0);
    rancorn.setVal(0.0);

    //nodal arrays used for calculating viscous stress
    std::array< MultiFab, AMREX_SPACEDIM > cornx;
    AMREX_D_TERM(cornx[0].define(convert(ba,nodal_flag), dmap, 1, 0);,
                 cornx[1].define(convert(ba,nodal_flag), dmap, 1, 0);,
                 cornx[2].define(convert(ba,nodal_flag), dmap, 1, 0););

    std::array< MultiFab, AMREX_SPACEDIM > corny;
    AMREX_D_TERM(corny[0].define(convert(ba,nodal_flag), dmap, 1, 0);,
                 corny[1].define(convert(ba,nodal_flag), dmap, 1, 0);,
                 corny[2].define(convert(ba,nodal_flag), dmap, 1, 0););

    std::array< MultiFab, AMREX_SPACEDIM > cornz;
    AMREX_D_TERM(cornz[0].define(convert(ba,nodal_flag), dmap, 1, 0);,
                 cornz[1].define(convert(ba,nodal_flag), dmap, 1, 0);,
                 cornz[2].define(convert(ba,nodal_flag), dmap, 1, 0););

    MultiFab visccorn;
    visccorn.define(convert(ba,nodal_flag), dmap, 1, 0);

    // for energy dissipation calculation
    // [du/dx, dv/dy, dw/dz]
    MultiFab gradU;
    MultiFab rhoscaled_gradU;
    // Laplacian
    MultiFab LapU;
    // [du/dx du/dy du/dz dv/dx dv/dy dv/dz dw/dx dw/dy dw/dx]
    MultiFab gradUtensor;
    MultiFab rhoscaled_gradUtensor;
    // temporary storage
    MultiFab ccTemp;

    std::array< MultiFab, AMREX_SPACEDIM > gradUtensor_fc;
    std::array< MultiFab, AMREX_SPACEDIM > rhoscaled_gradUtensor_fc;
    std::array< MultiFab, AMREX_SPACEDIM > temp_fc;
    std::array< MultiFab, AMREX_SPACEDIM > rho_fc;

    Real dProb;
    Real p0;
    Real rho0;
    Real nu0;

    if (turbForcing == 1) {
        gradU.define(ba,dmap,AMREX_SPACEDIM,0);
        rhoscaled_gradU.define(ba,dmap,AMREX_SPACEDIM,0);
        LapU.define(ba,dmap,AMREX_SPACEDIM,0);
        gradUtensor.define(ba,dmap,AMREX_SPACEDIM*AMREX_SPACEDIM,0);
        rhoscaled_gradUtensor.define(ba,dmap,AMREX_SPACEDIM*AMREX_SPACEDIM,0);
        ccTemp.define(ba,dmap,1,0);

        AMREX_D_TERM(gradUtensor_fc[0].define(convert(ba,nodal_flag_x), dmap, AMREX_SPACEDIM, 0);,
                     gradUtensor_fc[1].define(convert(ba,nodal_flag_y), dmap, AMREX_SPACEDIM, 0);,
                     gradUtensor_fc[2].define(convert(ba,nodal_flag_z), dmap, AMREX_SPACEDIM, 0););
        AMREX_D_TERM(rhoscaled_gradUtensor_fc[0].define(convert(ba,nodal_flag_x), dmap, AMREX_SPACEDIM, 0);,
                     rhoscaled_gradUtensor_fc[1].define(convert(ba,nodal_flag_y), dmap, AMREX_SPACEDIM, 0);,
                     rhoscaled_gradUtensor_fc[2].define(convert(ba,nodal_flag_z), dmap, AMREX_SPACEDIM, 0););
        AMREX_D_TERM(temp_fc[0].define(convert(ba,nodal_flag_x), dmap, 1, 0);,
                     temp_fc[1].define(convert(ba,nodal_flag_y), dmap, 1, 0);,
                     temp_fc[2].define(convert(ba,nodal_flag_z), dmap, 1, 0););
        AMREX_D_TERM(rho_fc[0].define(convert(ba,nodal_flag_x), dmap, 1, 0);,
                     rho_fc[1].define(convert(ba,nodal_flag_y), dmap, 1, 0);,
                     rho_fc[2].define(convert(ba,nodal_flag_z), dmap, 1, 0););

        p0 = 884.147e3;
        dProb = (AMREX_SPACEDIM==2) ? 1./(n_cells[0]*n_cells[1]) : 1./(n_cells[0]*n_cells[1]*n_cells[2]);
        rho0 = molmass[0] / avogadro * p0 / (k_B * T_init[0]);
        nu0 = 0.185;
    }

    Real time = 0;

    int step;
    int statsCount = 1;

    ///////////////////////////////////////////
    // Structure factor:
    ///////////////////////////////////////////

    // "primitive" variable structure factor will contain
    // rho
    // vel
    // T
    // pressure
    // Yk
    int structVarsPrim = AMREX_SPACEDIM+nspecies+3;

    Vector< std::string > prim_var_names;
    prim_var_names.resize(structVarsPrim);

    int cnt = 0;
    std::string x;

    // rho
    prim_var_names[cnt++] = "rho";

    // velx, vely, velz
    for (int d=0; d<AMREX_SPACEDIM; d++) {
        x = "vel";
        x += (120+d);
        prim_var_names[cnt++] = x;
    }

    // Temp
    prim_var_names[cnt++] = "Temp";

    // Pressure
    prim_var_names[cnt++] = "Pressure";

    // Yk
    for (int d=0; d<nspecies; d++) {
        x = "Y";
        x += (49+d);
        prim_var_names[cnt++] = x;
    }

    MultiFab structFactPrimMF;
    structFactPrimMF.define(ba, dmap, structVarsPrim, 0);
    structFactPrimMF.setVal(0.0);

    // scale SF results by inverse cell volume
    Vector<Real> var_scaling_prim;
    var_scaling_prim.resize(structVarsPrim*(structVarsPrim+1)/2);
    for (int d=0; d<var_scaling_prim.size(); ++d) {
        var_scaling_prim[d] = 1./(dx[0]*dx[1]*dx[2]);
    }

    // compute all pairs
    // note: StructFactPrim option to compute only speicified pairs not written yet
    StructFact structFactPrim(ba,dmap,prim_var_names,var_scaling_prim);

    ///////////////////////////////////////////

    // structure factor class for flattened dataset
    StructFact structFactPrimFlattened;

    //////////////////////////////////////////////

    // "conserved" variable structure factor will contain
    // rho
    // j
    // rho*E
    // rho*Yk
    // Temperature (not in the conserved array; will have to copy it in)
    int structVarsCons = AMREX_SPACEDIM+nspecies+3;

    Vector< std::string > cons_var_names;
    cons_var_names.resize(structVarsCons);

    cnt = 0;

    // rho
    cons_var_names[cnt++] = "rho";

    // velx, vely, velz
    for (int d=0; d<AMREX_SPACEDIM; d++) {
        x = "j";
        x += (120+d);
        cons_var_names[cnt++] = x;
    }

    // rho*E
    cons_var_names[cnt++] = "rhoE";

    // rho*Yk
    for (int d=0; d<nspecies; d++) {
        x = "rhoY";
        x += (49+d);
        cons_var_names[cnt++] = x;
    }

    // Temp
    cons_var_names[cnt++] = "Temp";

    MultiFab structFactConsMF;
    structFactConsMF.define(ba, dmap, structVarsCons, 0);
    structFactConsMF.setVal(0.0);

    // scale SF results by inverse cell volume
    Vector<Real> var_scaling_cons;
    var_scaling_cons.resize(structVarsCons*(structVarsCons+1)/2);
    for (int d=0; d<var_scaling_cons.size(); ++d) {
        var_scaling_cons[d] = 1./(dx[0]*dx[1]*dx[2]);
    }

    // compute all pairs
    // note: StructFactCons option to compute only speicified pairs not written yet
    StructFact structFactCons(ba,dmap,cons_var_names,var_scaling_cons);

    //////////////////////////////////////////////

    // structure factor class for flattened dataset
    StructFact structFactConsFlattened;

    //////////////////////////////////////////////

    if(project_dir >= 0){
        MultiFab Flattened;  // flattened multifab defined below

        // we are only calling ComputeVerticalAverage or ExtractSlice here to obtain
        // a built version of primFlattened so can obtain what we need to build the
        // structure factor and geometry objects for flattened data
        if (slicepoint < 0) {
            ComputeVerticalAverage(structFactPrimMF, Flattened, project_dir, 0, 1);
        } else {
            ExtractSlice(structFactPrimMF, Flattened, project_dir, slicepoint, 0, 1);
        }
        BoxArray ba_flat = Flattened.boxArray();
        const DistributionMapping& dmap_flat = Flattened.DistributionMap();

        structFactPrimFlattened.define(ba_flat,dmap_flat,prim_var_names,var_scaling_prim);
        structFactConsFlattened.define(ba_flat,dmap_flat,cons_var_names,var_scaling_cons);
    }

    ///////////////////////////////////////////
    // Structure factor object to help compute tubulent energy spectra
    ///////////////////////////////////////////

    // option to compute only specified pairs
    amrex::Vector< int > s_pairA(AMREX_SPACEDIM);
    amrex::Vector< int > s_pairB(AMREX_SPACEDIM);

    // need to use dVol for scaling
    Real dVol = (AMREX_SPACEDIM==2) ? dx[0]*dx[1]*cell_depth : dx[0]*dx[1]*dx[2];

    Vector< std::string > var_names;
    var_names.resize(AMREX_SPACEDIM);

    var_names[0] = "xvel";
    var_names[1] = "yvel";
#if (AMREX_SPACEDIM == 3)
    var_names[2] = "zvel";
#endif

    Vector<Real> var_scaling(AMREX_SPACEDIM);
    for (int d=0; d<var_scaling.size(); ++d) {
        var_scaling[d] = 1./dVol;
    }

    // Select which variable pairs to include in structure factor:
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        s_pairA[d] = d;
        s_pairB[d] = d;
    }
    StructFact turbStructFact;
    MultiFab structFactMF;
    if (turbForcing == 1) {
        turbStructFact.define(ba,dmap,var_names,var_scaling,s_pairA,s_pairB);
        structFactMF.define(ba, dmap, AMREX_SPACEDIM, 0);
    }

    //////////////////////////////////////////////

    // initialize conserved variables
    InitConsVar(cu,geom);

    // initialize primitive variables
    conservedToPrimitive(prim, cu);

    if (n_ads_spec>0) init_surfcov(surfcov, geom);

    // Set BC: 1) fill boundary 2) physical
    cu.FillBoundary(geom.periodicity());
    prim.FillBoundary(geom.periodicity());
    setBC(prim, cu);

    if (plot_int > 0) {
        WritePlotFile(0, 0.0, geom, cu, cuMeans, cuVars,
                      prim, primMeans, primVars, spatialCross, eta, kappa);
    }

#ifdef MUI
    // MUI setting
    mui::uniface2d uniface( "mpi://FHD-side/FHD-KMC-coupling" );

    mui_announce_send_recv_span(uniface,cu,dx);
#endif

    //Time stepping loop
    for(step=1;step<=max_step;++step) {

        if (restart > 0 && step==1) {
            ReadCheckPoint(step, time, statsCount, geom, cu, cuMeans, cuVars, prim,
                           primMeans, primVars, spatialCross, miscStats, eta, kappa);
        }

        // timer
        Real ts1 = ParallelDescriptor::second();

        // sample surface chemistry (via either surfchem_mui or MFsurfchem)
#ifdef MUI
        mui_push(cu, prim, dx, uniface, step);
#endif
        if (n_ads_spec>0) sample_MFsurfchem(cu, prim, surfcov, dNadsdes, dNads, dNdes, geom, dt);

        // FHD
        RK3step(cu, cup, cup2, cup3, prim, source, eta, zeta, kappa, chi, D, flux,
                stochFlux, cornx, corny, cornz, visccorn, rancorn, ranchem, geom, dt);

        // update surface chemistry (via either surfchem_mui or MFsurfchem)
#ifdef MUI
        mui_fetch(cu, prim, dx, uniface, step);

        conservedToPrimitive(prim, cu);

        // Set BC: 1) fill boundary 2) physical
        cu.FillBoundary(geom.periodicity());
        prim.FillBoundary(geom.periodicity());
        setBC(prim, cu);
#endif
        if (n_ads_spec>0) {

            update_MFsurfchem(cu, prim, surfcov, dNadsdes, dNads, dNdes, geom);

            conservedToPrimitive(prim, cu);

            // Set BC: 1) fill boundary 2) physical
            cu.FillBoundary(geom.periodicity());
            prim.FillBoundary(geom.periodicity());
            setBC(prim, cu);
        }

        // timer
        Real ts2 = ParallelDescriptor::second() - ts1;
        ParallelDescriptor::ReduceRealMax(ts2);
        if (step%100 == 0) {
            amrex::Print() << "Advanced step " << step << " in " << ts2 << " seconds\n";
        }

        // compute mean and variances
        if (step > n_steps_skip && stats_int > 0 && step%stats_int == 0) {

            // timer
            Real t1 = ParallelDescriptor::second();

            evaluateStats(cu, cuMeans, cuVars, prim, primMeans, primVars,
                          spatialCross, miscStats, miscVals, statsCount, dx);
            statsCount++;

            // timer
            Real t2 = ParallelDescriptor::second() - t1;
            ParallelDescriptor::ReduceRealMax(t2);
            if (step%100 == 0) {
                amrex::Print() << "evaluateStats time " << t2 << " seconds\n";
            }
        }
        if (step%100 == 0) {
            amrex::Print() << "Mean Momentum (x, y, z): " << ComputeSpatialMean(cu, 1) << " " << ComputeSpatialMean(cu, 2) << " " << ComputeSpatialMean(cu, 3) << "\n";
        }

        // write a plotfile
        if (plot_int > 0 && step%plot_int == 0) {

            // timer
            Real t1 = ParallelDescriptor::second();

            /*
              yzAverage(cuMeans, cuVars, primMeans, primVars, spatialCross,
                        cuMeansAv, cuVarsAv, primMeansAv, primVarsAv, spatialCrossAv);
              WritePlotFile(step, time, geom, cu, cuMeansAv, cuVarsAv,
                            prim, primMeansAv, primVarsAv, spatialCrossAv, eta, kappa);
            */
            WritePlotFile(step, time, geom, cu, cuMeans, cuVars,
                         prim, primMeans, primVars, spatialCross, eta, kappa);

#ifdef MUI
            // also horizontal average
            WriteHorizontalAverage(cu,2,0,5+nspecies,step,geom);
#endif
            if (n_ads_spec>0) WriteHorizontalAverage(cu,2,0,5+nspecies,step,geom);

            // timer
            Real t2 = ParallelDescriptor::second() - t1;
            ParallelDescriptor::ReduceRealMax(t2);
            amrex::Print() << "WritePlotFile time " << t2 << " seconds\n";

            // snapshot of instantaneous energy spectra
            if (turbForcing == 1) {

                // timer
                t1 = ParallelDescriptor::second();

                // copy velocities into structFactMF
                MultiFab::Copy(structFactMF, prim, 1, 0, AMREX_SPACEDIM, 0);

                // reset and compute structure factor
                turbStructFact.FortStructure(structFactMF,1);
                turbStructFact.CallFinalize();

                // integrate cov_mag over shells in k and write to file
                turbStructFact.IntegratekShells(step);

                // timer
                t2 = ParallelDescriptor::second() - t1;
                ParallelDescriptor::ReduceRealMax(t2);
                amrex::Print() << "Ek time " << t2 << " seconds\n";
            }
        }

        if (chk_int > 0 && step > 0 && step%chk_int == 0) {

            // timer
            Real t1 = ParallelDescriptor::second();

            WriteCheckPoint(step, time, statsCount, geom, cu, cuMeans,
                            cuVars, prim, primMeans, primVars, spatialCross, miscStats, eta, kappa);

            // timer
            Real t2 = ParallelDescriptor::second() - t1;
            ParallelDescriptor::ReduceRealMax(t2);
            amrex::Print() << "WriteCheckPoint time " << t2 << " seconds\n";
        }

        // collect a snapshot for structure factor
        if (step > n_steps_skip && struct_fact_int > 0 && (step-n_steps_skip)%struct_fact_int == 0) {

            // timer
            Real t1 = ParallelDescriptor::second();

            MultiFab::Copy(structFactPrimMF, prim, 0,                0,                structVarsPrim,   0);
            MultiFab::Copy(structFactConsMF, cu,   0,                0,                structVarsCons-1, 0);
            MultiFab::Copy(structFactConsMF, prim, AMREX_SPACEDIM+1, structVarsCons-1, 1,                0); // temperature too
            structFactPrim.FortStructure(structFactPrimMF);
            structFactCons.FortStructure(structFactConsMF);
            if(project_dir >= 0) {
                MultiFab primFlattened;  // flattened multifab defined below
                MultiFab consFlattened;  // flattened multifab defined below
                if (slicepoint < 0) {
                    ComputeVerticalAverage(structFactPrimMF, primFlattened, project_dir, 0, structVarsPrim);
                    ComputeVerticalAverage(structFactConsMF, consFlattened, project_dir, 0, structVarsCons);
                } else {
                    ExtractSlice(structFactPrimMF, primFlattened, project_dir, slicepoint, 0, structVarsPrim);
                    ExtractSlice(structFactConsMF, consFlattened, project_dir, slicepoint, 0, structVarsCons);
                }
                structFactPrimFlattened.FortStructure(primFlattened);
                structFactConsFlattened.FortStructure(consFlattened);
            }

            // timer
            Real t2 = ParallelDescriptor::second() - t1;
            ParallelDescriptor::ReduceRealMax(t2);
            amrex::Print() << "StructFact snapshot time " << t2 << " seconds\n";
        }

        // write out structure factor
        if (step > n_steps_skip && struct_fact_int > 0 && plot_int > 0 && step%plot_int == 0) {

            // timer
            Real t1 = ParallelDescriptor::second();

            Print() << "HERE1\n";

            structFactPrim.WritePlotFile(step,time,"plt_SF_prim");

            Print() << "HERE2\n";
            structFactCons.WritePlotFile(step,time,"plt_SF_cons");
            if(project_dir >= 0) {
                structFactPrimFlattened.WritePlotFile(step,time,"plt_SF_prim_Flattened");
                structFactConsFlattened.WritePlotFile(step,time,"plt_SF_cons_Flattened");
            }

            // timer
            Real t2 = ParallelDescriptor::second() - t1;
            ParallelDescriptor::ReduceRealMax(t2);
            amrex::Print() << "StructFact plotfile time " << t2 << " seconds\n";
        }

        // energy dissipation rate
        if (turbForcing == 1) {

            // timer
            Real t1 = ParallelDescriptor::second();

            // FORM 1: <rho/rho0 du_i/dx_i du_i/dx_i>

            // compute gradU = [du/dx, dv/dy, dw/dz]
            for (int d=0; d<AMREX_SPACEDIM; ++d) {
                ComputeCentredGradCompDir(prim,gradU,d,d+1,d,geom);
            }

            // create a copy of gradU scaled by rho
            MultiFab::Copy(rhoscaled_gradU, gradU, 0, 0, AMREX_SPACEDIM, 0);
            for (int d=0; d<AMREX_SPACEDIM; ++d) {
                MultiFab::Multiply(rhoscaled_gradU, prim, 0, d, 1, 0);
            }

            // compute <rho/rho0 du_i/dx_i du_i/dx_i>
            Vector<Real> gradUdotgradU(3);
            for (int d=0; d<AMREX_SPACEDIM; ++d) {
                CCInnerProd(gradU,d,rhoscaled_gradU,d,ccTemp,gradUdotgradU[d]);
            }

            Real FORM1 = dProb*(gradUdotgradU[0] + gradUdotgradU[1] + gradUdotgradU[2]) * (nu0 / rho0);

            // FORM 2: <-rho/rho0 u_j Lap(u_j)>

            // Lap(u_j)
            ComputeLap(prim,LapU,1,0,AMREX_SPACEDIM,geom);

            // <rho u_j Lap(u_j)>
            Vector<Real> rhoULapU(3);
            for (int d=0; d<AMREX_SPACEDIM; ++d) {
                CCInnerProd(cu,d+1,LapU,d,ccTemp,rhoULapU[d]);
            }

            // <rho/rho0 u_j Lap(u_j)>
            Real FORM2 = -dProb*(rhoULapU[0] + rhoULapU[1] + rhoULapU[2]) * (nu0 / rho0);

            // FORM 3: <rho/rho0 du_i/dx_j du_i/dx_j> using cell-centered gradients

            // compute [du/dx du/dy du/dz dv/dx dv/dy dv/dz dw/dx dw/dy dw/dx]
            for (int j=0; j<AMREX_SPACEDIM; ++j) {
                for (int i=0; i<AMREX_SPACEDIM; ++i) {
                    ComputeCentredGradCompDir(prim,gradUtensor,i,j+1,i+j*AMREX_SPACEDIM,geom);
                }
            }

            // create a copy of gradUtensor scaled by rho
            MultiFab::Copy(rhoscaled_gradUtensor, gradUtensor, 0, 0, AMREX_SPACEDIM*AMREX_SPACEDIM, 0);
            for (int d=0; d<AMREX_SPACEDIM*AMREX_SPACEDIM; ++d) {
                MultiFab::Multiply(rhoscaled_gradUtensor, prim, 0, d, 1, 0);
            }

            // compute <rho/rho0 du_i/dx_j du_i/dx_j>
            Vector<Real> gradUdotgradUtensor(9);
            for (int d=0; d<AMREX_SPACEDIM*AMREX_SPACEDIM; ++d) {
                CCInnerProd(gradUtensor,d,rhoscaled_gradUtensor,d,ccTemp,gradUdotgradUtensor[d]);
            }

            Real FORM3 = dProb*(  gradUdotgradUtensor[0] + gradUdotgradUtensor[1] + gradUdotgradUtensor[2]
                                  + gradUdotgradUtensor[3] + gradUdotgradUtensor[4] + gradUdotgradUtensor[5]
                                  + gradUdotgradUtensor[6] + gradUdotgradUtensor[7] + gradUdotgradUtensor[8]) * (nu0 / rho0);


            // FORM 4: <rho/rho0 du_i/dx_j du_i/dx_j> using face-centered gradients

            for (int d=0; d<AMREX_SPACEDIM; ++d) {
                // du/dx du/dy du/dz goes into component 0 of each element in gradUtensor_fc
                // dv/dx dv/dy dv/dz goes into component 1 of each element in gradUtensor_fc
                // dw/dx dw/dy dw/dz goes into component 2 of each element in gradUtensor_fc
                // only works for periodic (9999 is a fake bc_comp that makes everything INT_DIR)
                ComputeGrad(prim,gradUtensor_fc,d+1,d,1,9999,geom);
            }

            // compute rho at faces
            AverageCCToFace(prim,rho_fc,0,1,RHO_BC_COMP,geom);

            // create a copy of gradUtensor_fc scaled by rho
            for (int d=0; d<AMREX_SPACEDIM; ++d) {
                // copy all 3 components of gradUtensor_fc[d] into a temporary
                MultiFab::Copy(rhoscaled_gradUtensor_fc[d], gradUtensor_fc[d], 0, 0, AMREX_SPACEDIM, 0);
                for (int dd=0; dd<AMREX_SPACEDIM; ++dd) {
                    // scale each temporary value by rho
                    MultiFab::Multiply(rhoscaled_gradUtensor_fc[d], rho_fc[d], 0, dd, 1, 0);
                }
            }

            Vector<Real> gradUdotgradUx(3);
            Vector<Real> gradUdotgradUy(3);
            Vector<Real> gradUdotgradUz(3);
            StagInnerProd(rhoscaled_gradUtensor_fc,0,gradUtensor_fc,0,temp_fc,gradUdotgradUx);
            StagInnerProd(rhoscaled_gradUtensor_fc,1,gradUtensor_fc,1,temp_fc,gradUdotgradUy);
            StagInnerProd(rhoscaled_gradUtensor_fc,2,gradUtensor_fc,2,temp_fc,gradUdotgradUz);

            Real FORM4 = dProb*(  gradUdotgradUx[0] + gradUdotgradUx[1] + gradUdotgradUx[2]
                                  + gradUdotgradUy[0] + gradUdotgradUy[1] + gradUdotgradUy[2]
                                  + gradUdotgradUz[0] + gradUdotgradUz[1] + gradUdotgradUz[2]) * (nu0 / rho0);

            // FORM 5: <curl(V) dot (curl(V)> using cell-centered gradients

            // compute curlU (store in gradU)
            for (int d=0; d<AMREX_SPACEDIM; ++d) {
                ComputeCurlCC(prim,1,gradU,0,geom);
            }

            // create a copy of curlU scaled by rho
            MultiFab::Copy(rhoscaled_gradU, gradU, 0, 0, AMREX_SPACEDIM, 0);
            for (int d=0; d<AMREX_SPACEDIM; ++d) {
                MultiFab::Multiply(rhoscaled_gradU, prim, 0, d, 1, 0);
            }

            // compute <curl(V) dot (curl(V)>
            Vector<Real> curlUdotcurlU(3);
            for (int d=0; d<AMREX_SPACEDIM; ++d) {
                CCInnerProd(gradU,d,rhoscaled_gradU,d,ccTemp,curlUdotcurlU[d]);
            }

            Real FORM5 =  dProb*(curlUdotcurlU[0] + curlUdotcurlU[1] + curlUdotcurlU[2]) * (nu0 / rho0);

            curlUdotcurlU[0]=0;
            ComputeDivCC(prim,1,gradU,0,geom);
            MultiFab::Copy(rhoscaled_gradU, gradU, 0, 0, 1, 0);
            MultiFab::Multiply(rhoscaled_gradU, prim, 0, 0, 1, 0);
            CCInnerProd(gradU,0,rhoscaled_gradU,0,ccTemp,curlUdotcurlU[0]);

            Real DIVCOR = 2.* dProb*curlUdotcurlU[0] * nu0 / (3.* rho0);

            Print() << "Non-viscosity scaled energy dissipation "
                    << time << " "
                    << FORM1 << " "
                    << FORM2 << " "
                    << FORM3 << " "
                    << FORM4 << " "
                    << FORM5 << " "
                    << FORM3 - DIVCOR  << " "
                    << std::endl;

            // timer
            Real t2 = ParallelDescriptor::second() - t1;
            ParallelDescriptor::ReduceRealMax(t2);
            amrex::Print() << "Energy dissipation compute time " << t2 << " seconds\n";

        }  // end if (turbForcing == 1)

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
    ParallelDescriptor::ReduceRealMax(stop_time);
    amrex::Print() << "Run time = " << stop_time << std::endl;
}