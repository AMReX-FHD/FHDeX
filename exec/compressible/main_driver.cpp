#include "common_functions.H"
#include "compressible_functions.H"


#include "rng_functions.H"

#include "StructFact.H"

#include "chemistry_functions.H"

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

    
    // copy contents of F90 modules to C++ namespaces
    InitializeCommonNamespace();

    InitializeCompressibleNamespace();

    // read the inputs file for chemistry
    InitializeChemistryNamespace();

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

    //Initialize physical parameters from input vals

    double intEnergy, T0;

    T0 = T_init[0];

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

    Real time = 0;

    int step, statsCount;

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
    MultiFab primFlattenedRotMaster;

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
    MultiFab consFlattenedRotMaster;

    //////////////////////////////////////////////
    
    Geometry geom_flat;
    
    if(project_dir >= 0){
      MultiFab Flattened;  // flattened multifab defined below
      
      // we are only calling ComputeVerticalAverage or ExtractSlice here to obtain
      // a built version of primFlattened so can obtain what we need to build the
      // structure factor and geometry objects for flattened data
      if (slicepoint < 0) {
          ComputeVerticalAverage(structFactPrimMF, Flattened, geom, project_dir, 0, 1);
      } else {
          ExtractSlice(structFactPrimMF, Flattened, geom, project_dir, slicepoint, 0, 1);
      }
      // we rotate this flattened MultiFab to have normal in the z-direction since
      // our structure factor class assumes this for flattened
      MultiFab FlattenedRot = RotateFlattenedMF(Flattened);
      BoxArray ba_flat = FlattenedRot.boxArray();
      const DistributionMapping& dmap_flat = FlattenedRot.DistributionMap();
      primFlattenedRotMaster.define(ba_flat,dmap_flat,structVarsPrim,0);
      consFlattenedRotMaster.define(ba_flat,dmap_flat,structVarsCons,0);
      {
        IntVect dom_lo(AMREX_D_DECL(0,0,0));
        IntVect dom_hi;

        // yes you could simplify this code but for now
        // these are written out fully to better understand what is happening
        // we wanted dom_hi[AMREX_SPACEDIM-1] to be equal to 0
        // and need to transmute the other indices depending on project_dir
#if (AMREX_SPACEDIM == 2)
        if (project_dir == 0) {
            dom_hi[0] = n_cells[1]-1;
        }
        else if (project_dir == 1) {
            dom_hi[0] = n_cells[0]-1;
        }
        dom_hi[1] = 0;
#elif (AMREX_SPACEDIM == 3)
        if (project_dir == 0) {
            dom_hi[0] = n_cells[1]-1;
            dom_hi[1] = n_cells[2]-1;
        } else if (project_dir == 1) {
            dom_hi[0] = n_cells[0]-1;
            dom_hi[1] = n_cells[2]-1;
        } else if (project_dir == 2) {
            dom_hi[0] = n_cells[0]-1;
            dom_hi[1] = n_cells[1]-1;
        }
        dom_hi[2] = 0;
#endif
        Box domain(dom_lo, dom_hi);

        // This defines the physical box
        Vector<Real> projected_hi(AMREX_SPACEDIM);

        // yes you could simplify this code but for now
        // these are written out fully to better understand what is happening
        // we wanted projected_hi[AMREX_SPACEDIM-1] to be equal to dx[projected_dir]
        // and need to transmute the other indices depending on project_dir
#if (AMREX_SPACEDIM == 2)
        if (project_dir == 0) {
            projected_hi[0] = prob_hi[1];
        } else if (project_dir == 1) {
            projected_hi[0] = prob_hi[0];
        }
        projected_hi[1] = prob_hi[project_dir] / n_cells[project_dir];
#elif (AMREX_SPACEDIM == 3)
        if (project_dir == 0) {
            projected_hi[0] = prob_hi[1];
            projected_hi[1] = prob_hi[2];
        } else if (project_dir == 1) {
            projected_hi[0] = prob_hi[0];
            projected_hi[1] = prob_hi[2];
        } else if (project_dir == 2) {
            projected_hi[0] = prob_hi[0];
            projected_hi[1] = prob_hi[1];
        }
        projected_hi[2] = prob_hi[project_dir] / n_cells[project_dir];
#endif

        RealBox real_box({AMREX_D_DECL(     prob_lo[0],     prob_lo[1],     prob_lo[2])},
                         {AMREX_D_DECL(projected_hi[0],projected_hi[1],projected_hi[2])});
        
        // This defines a Geometry object
        geom_flat.define(domain,&real_box,CoordSys::cartesian,is_periodic.data());
      }

      structFactPrimFlattened.define(ba_flat,dmap_flat,prim_var_names,var_scaling_prim);
      structFactConsFlattened.define(ba_flat,dmap_flat,cons_var_names,var_scaling_cons);
    }
    
    //////////////////////////////////////////////

    // Initialize everything
    
    prim.setVal(0.0,0,nprimvars,ngc);
    prim.setVal(rho0,0,1,ngc);      // density
    prim.setVal(0.,1,3,ngc);        // x/y/z velocity
    prim.setVal(T_init[0],4,1,ngc); // temperature
                                    // pressure computed later in conservedToPrimitive
    for(int i=0;i<nspecies;i++) {
        prim.setVal(rhobar[i],6+i,1,ngc);    // mass fractions
    }

    // compute internal energy
    GpuArray<Real,MAX_SPECIES> massvec;
    for(int i=0;i<nspecies;i++) {
        massvec[i] = rhobar[i];
    }
    GetEnergy(intEnergy, massvec, T0);

    cu.setVal(0.0,0,nvars,ngc);
    cu.setVal(rho0,0,1,ngc);           // density
    cu.setVal(0,1,3,ngc);              // x/y/z momentum
    cu.setVal(rho0*intEnergy,4,1,ngc); // total energy
    for(int i=0;i<nspecies;i++) {
        cu.setVal(rho0*rhobar[i],5+i,1,ngc); // mass densities
    }

    // RK stage storage
    cup.setVal(0.0,0,nvars,ngc);
    cup2.setVal(0.0,0,nvars,ngc);
    cup3.setVal(0.0,0,nvars,ngc);

    // set density
    cup.setVal(rho0,0,1,ngc);
    cup2.setVal(rho0,0,1,ngc);
    cup3.setVal(rho0,0,1,ngc);

    // initialize conserved variables
    if (prob_type > 1) {
        InitConsVar(cu,prim,geom);
    }

    statsCount = 1;

    // Write initial plotfile
    conservedToPrimitive(prim, cu);

    // Set BC: 1) fill boundary 2) physical
    cu.FillBoundary(geom.periodicity());
    prim.FillBoundary(geom.periodicity());
    setBC(prim, cu);
    
    if (plot_int > 0) {
        WritePlotFile(0, 0.0, geom, cu, cuMeans, cuVars,
                      prim, primMeans, primVars, spatialCross, eta, kappa);
    }

    //Time stepping loop
    for(step=1;step<=max_step;++step) {

        if (restart > 0 && step==1) {
            ReadCheckPoint(step, time, statsCount, geom, cu, cuMeans, cuVars, prim,
                           primMeans, primVars, spatialCross, miscStats, eta, kappa);
        }

        // timer
        Real ts1 = ParallelDescriptor::second();

        RK3step(cu, cup, cup2, cup3, prim, source, eta, zeta, kappa, chi, D, flux,
                stochFlux, cornx, corny, cornz, visccorn, rancorn, ranchem, geom, dt);

        // timer
        Real ts2 = ParallelDescriptor::second() - ts1;
        ParallelDescriptor::ReduceRealMax(ts2);
    	amrex::Print() << "Advanced step " << step << " in " << ts2 << " seconds\n";

        // timer
        Real aux1 = ParallelDescriptor::second();
        
        // compute mean and variances
        if (step > n_steps_skip) {
            evaluateStats(cu, cuMeans, cuVars, prim, primMeans, primVars,
                          spatialCross, miscStats, miscVals, statsCount, dx);
            statsCount++;
        }

        // write a plotfile
        if (plot_int > 0 && step > 0 && step%plot_int == 0) {
        /*
           yzAverage(cuMeans, cuVars, primMeans, primVars, spatialCross,
                     cuMeansAv, cuVarsAv, primMeansAv, primVarsAv, spatialCrossAv);
           WritePlotFile(step, time, geom, cu, cuMeansAv, cuVarsAv,
                         prim, primMeansAv, primVarsAv, spatialCrossAv, eta, kappa);
        */
           WritePlotFile(step, time, geom, cu, cuMeans, cuVars,
                         prim, primMeans, primVars, spatialCross, eta, kappa);
        }

        if (chk_int > 0 && step > 0 && step%chk_int == 0)
        {
           WriteCheckPoint(step, time, statsCount, geom, cu, cuMeans,
                           cuVars, prim, primMeans, primVars, spatialCross, miscStats, eta, kappa);
        }

	// collect a snapshot for structure factor
	if (step > n_steps_skip && struct_fact_int > 0 && (step-n_steps_skip)%struct_fact_int == 0) {
            MultiFab::Copy(structFactPrimMF, prim, 0,                0,                structVarsPrim,   0);
            MultiFab::Copy(structFactConsMF, cu,   0,                0,                structVarsCons-1, 0);
            MultiFab::Copy(structFactConsMF, prim, AMREX_SPACEDIM+1, structVarsCons-1, 1,                0); // temperature too
            structFactPrim.FortStructure(structFactPrimMF,geom);
            structFactCons.FortStructure(structFactConsMF,geom);
            if(project_dir >= 0) {
                MultiFab primFlattened;  // flattened multifab defined below
                MultiFab consFlattened;  // flattened multifab defined below
                if (slicepoint < 0) {
                    ComputeVerticalAverage(structFactPrimMF, primFlattened, geom, project_dir, 0, structVarsPrim);
                    ComputeVerticalAverage(structFactConsMF, consFlattened, geom, project_dir, 0, structVarsCons);
                } else {
                    ExtractSlice(structFactPrimMF, primFlattened, geom, project_dir, slicepoint, 0, structVarsPrim);
                    ExtractSlice(structFactConsMF, consFlattened, geom, project_dir, slicepoint, 0, structVarsCons);
                }
                // we rotate this flattened MultiFab to have normal in the z-direction since
                // our structure factor class assumes this for flattened
                MultiFab primFlattenedRot = RotateFlattenedMF(primFlattened);
                primFlattenedRotMaster.ParallelCopy(primFlattenedRot,0,0,structVarsPrim);
                structFactPrimFlattened.FortStructure(primFlattenedRotMaster,geom_flat);

                MultiFab consFlattenedRot = RotateFlattenedMF(consFlattened);
                consFlattenedRotMaster.ParallelCopy(consFlattenedRot,0,0,structVarsCons);
                structFactConsFlattened.FortStructure(consFlattenedRotMaster,geom_flat);
            }
        }

        // write out structure factor
        if (step > n_steps_skip && struct_fact_int > 0 && plot_int > 0 && step%plot_int == 0) {
            structFactPrim.WritePlotFile(step,time,geom,"plt_SF_prim");
            structFactCons.WritePlotFile(step,time,geom,"plt_SF_cons");
            if(project_dir >= 0) {
                structFactPrimFlattened.WritePlotFile(step,time,geom_flat,"plt_SF_prim_Flattened");
                structFactConsFlattened.WritePlotFile(step,time,geom_flat,"plt_SF_cons_Flattened");
            }
        }
        
        // timer
        Real aux2 = ParallelDescriptor::second() - aux1;
        ParallelDescriptor::ReduceRealMax(aux2);
        amrex::Print() << "Aux time (stats, struct fac, plotfiles) " << aux2 << " seconds\n";
        
        time = time + dt;

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

    // timer
    Real stop_time = ParallelDescriptor::second() - strt_time;
    ParallelDescriptor::ReduceRealMax(stop_time);
    amrex::Print() << "Run time = " << stop_time << std::endl;
}
