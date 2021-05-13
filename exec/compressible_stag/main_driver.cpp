#include "common_functions.H"
#include "compressible_functions.H"
#include "compressible_functions_stag.H"

#include "common_namespace_declarations.H"
#include "compressible_namespace_declarations.H"

#include "rng_functions.H"

#include "StructFact.H"

using namespace amrex;

// argv contains the name of the inputs file entered at the command line
void main_driver(const char* argv)
{
    BL_PROFILE_VAR("main_driver()",main_driver);
    
    // store the current time so we can later compute total run time.
    Real strt_time = ParallelDescriptor::second();

    std::string inputs_file = argv;

    // read in parameters from inputs file into F90 modules
    // we use "+1" because of amrex_string_c_to_f expects a null char termination
    read_common_namelist(inputs_file.c_str(),inputs_file.size()+1);
    read_compressible_namelist(inputs_file.c_str(),inputs_file.size()+1);
    
    // copy contents of F90 modules to C++ namespaces
    InitializeCommonNamespace();
    InitializeCompressibleNamespace();

    int step_start, statsCount;
    amrex::Real time;

    // if gas heat capacities in the namelist are negative, calculate them using using dofs.
    // This will only update the Fortran values.
    get_hc_gas();
    // now update C++ values
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
    setup_bc(); // do the same in the fortran namelist

    // if multispecies
    if (algorithm_type == 2) {
        // compute wall concentrations if BCs call for it
        setup_cwall();
        SetupCWall();
    }

    /////////////////////////////////////////
    //Initialise rngs
    /////////////////////////////////////////

    if (restart < 0) {
        // NOTE: only fhdSeed is used currently
        // zero is a clock-based seed
        int fhdSeed      = seed;
        int particleSeed = seed;
        int selectorSeed = seed;
        int thetaSeed    = seed;
        int phiSeed      = seed;
        int generalSeed  = seed;

        // "seed" controls all of them and gives distinct seeds to each physical process over each MPI process
        // this should be fixed so each physical process has its own seed control
        if (seed > 0) {
            fhdSeed      = 6*ParallelDescriptor::MyProc() + seed;
            particleSeed = 6*ParallelDescriptor::MyProc() + seed + 1;
            selectorSeed = 6*ParallelDescriptor::MyProc() + seed + 2;
            thetaSeed    = 6*ParallelDescriptor::MyProc() + seed + 3;
            phiSeed      = 6*ParallelDescriptor::MyProc() + seed + 4;
            generalSeed  = 6*ParallelDescriptor::MyProc() + seed + 5;
        }

        // Initialise rngs
        rng_initialize(&fhdSeed,&particleSeed,&selectorSeed,&thetaSeed,&phiSeed,&generalSeed);

        // initializes the seed for C++ random number calls
        InitRandom(seed+ParallelDescriptor::MyProc());
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
    MultiFab eta;
    MultiFab zeta;
    MultiFab kappa;
    MultiFab chi;
    MultiFab D;

    // conserved quantaties
    MultiFab cu;

    // staggered momentum
    std::array< MultiFab, AMREX_SPACEDIM > vel;
    std::array< MultiFab, AMREX_SPACEDIM > cumom;
  
    //primative quantaties
    MultiFab prim;

    //statistics    
    MultiFab cuMeans;
    MultiFab cuVars;
    
    MultiFab primMeans;
    MultiFab primVars;

    MultiFab coVars;

    std::array< MultiFab, AMREX_SPACEDIM > velMeans;
    std::array< MultiFab, AMREX_SPACEDIM > velVars;
    std::array< MultiFab, AMREX_SPACEDIM > cumomMeans;
    std::array< MultiFab, AMREX_SPACEDIM > cumomVars;
    
    if ((plot_cross) and ((cross_cell <= 0) or (cross_cell >= n_cells[0]-1))) {
        Abort("Cross cell needs to be within the domain: 0 < cross_cell < n_cells[0] - 1");
    }

    // contains yz-averaged running & instantaneous averages of conserved variables at cross cell + four primitive variables [vx, vy, vz, T]: 2*nvars + 2*4
    Vector<Real> yzAvMeans_cross(2*nvars+8, 0.0); 
    
    // 1: <delrho*delrho>
    // 2: <delrhoE*delrhoE>
    // 3: <deljx*deljx>
    // 4: <deljy*deljy>
    // 5: <deljz*deljz>
    // 6: <deljx*delrho>
    // 7: <delG*delG>
    // 8: <delG*delK>
    // 9: <delK*delG>
    // 10: <delrho*delK>
    // 11: <delK*delrho>
    // 12: <delrho*delG>
    // 13: <delG*delrho>
    // 14: <delT*delT>
    // 15: <delT*delrho>
    // 16: <delux*delrho>
    // nspecies: <delrhoYk*delrhoYk>
    // can add more -- change main_driver, statsStag, writeplotfilestag, and Checkpoint
    int ncross = 16+nspecies;
    Vector<Real> spatialCross(n_cells[0]*ncross, 0.0); 
    
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
    geom.define(domain,&real_box,CoordSys::cartesian,is_periodic.data());

    Real dt = fixed_dt;
    const Real* dx = geom.CellSize();
    const RealBox& realDomain = geom.ProbDomain();

    std::string filename = "crossMeans";
    std::ofstream outfile;

    /////////////////////////////////////////////
    // Setup Structure factor variables & scaling
    /////////////////////////////////////////////
    MultiFab structFactPrimMF;
    MultiFab structFactConsMF;

    StructFact structFactPrim;
    StructFact structFactCons;
    StructFact structFactPrimVerticalAverage;
    StructFact structFactConsVerticalAverage;
    
    Geometry geom_flat;

    // "primitive" variable structure factor will contain
    // rho
    // vel (shifted)
    // T
    // Yk
    // vel (averaged)
    int structVarsPrim = 2*AMREX_SPACEDIM+nspecies+2;

    Vector< std::string > prim_var_names;
    prim_var_names.resize(structVarsPrim);

    int cnt = 0;
    std::string x;

    // rho
    prim_var_names[cnt++] = "rho";

    // velx, vely, velz
    for (int d=0; d<AMREX_SPACEDIM; d++) {
        x = "velCC";
        x += (120+d);
        prim_var_names[cnt++] = x;
    }

    // Temp
    prim_var_names[cnt++] = "Temp";

    // Yk
    for (int d=0; d<nspecies; d++) {
        x = "Y";
        x += (49+d);
        prim_var_names[cnt++] = x;
    }

    // velx, vely, velz
    for (int d=0; d<AMREX_SPACEDIM; d++) {
        x = "velFACE";
        x += (120+d);
        prim_var_names[cnt++] = x;
    }

    // "conserved" variable structure factor will contain
    // rho
    // j (averaged)
    // rho*E
    // rho*Yk
    // Temperature (not in the conserved array; will have to copy it in)
    // j (shifted)
    int structVarsCons = 2*AMREX_SPACEDIM+nspecies+3;

    Vector< std::string > cons_var_names;
    cons_var_names.resize(structVarsCons);

    cnt = 0;

    // rho
    cons_var_names[cnt++] = "rho";

    // jx, jy, jz
    for (int d=0; d<AMREX_SPACEDIM; d++) {
        x = "jCC";
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

    // jx, jy, jz
    for (int d=0; d<AMREX_SPACEDIM; d++) {
        x = "jFACE";
        x += (120+d);
        cons_var_names[cnt++] = x;
    }

    // scale SF results by inverse cell volume
    Vector<Real> var_scaling_prim;
    var_scaling_prim.resize(structVarsPrim*(structVarsPrim+1)/2);
    for (int d=0; d<var_scaling_prim.size(); ++d) {
        var_scaling_prim[d] = 1./(dx[0]*dx[1]*dx[2]);
    }
    Vector<Real> var_scaling_cons;
    // scale SF results by inverse cell volume
    var_scaling_cons.resize(structVarsCons*(structVarsCons+1)/2);
    for (int d=0; d<var_scaling_cons.size(); ++d) {
        var_scaling_cons[d] = 1./(dx[0]*dx[1]*dx[2]);
    }


    /////////////////////////////////////////////
    // Initialize based on fresh start or restart
    /////////////////////////////////////////////

    if (restart > 0) {
        ReadCheckPoint(step_start, time, statsCount, geom, cu, cuMeans, cuVars, prim,
                       primMeans, primVars, cumom, cumomMeans, cumomVars, 
                       vel, velMeans, velVars, coVars, spatialCross,
                       ba, dmap);

        // transport properties
        eta.define(ba,dmap,1,ngc);
        zeta.define(ba,dmap,1,ngc);
        kappa.define(ba,dmap,1,ngc);
        chi.define(ba,dmap,nspecies,ngc);
        D.define(ba,dmap,nspecies*nspecies,ngc);

        eta.setVal(1.0,0,1,ngc);
        zeta.setVal(1.0,0,1,ngc);
        kappa.setVal(1.0,0,1,ngc);
        chi.setVal(1.0,0,nspecies,ngc);
        D.setVal(1.0,0,nspecies*nspecies,ngc);

        if (plot_cross) {
            if (ParallelDescriptor::IOProcessor()) outfile.open(filename, std::ios::app);
        }

        ///////////////////////////////////////////
        // Setup Structure factor
        ///////////////////////////////////////////

        structFactPrimMF.define(ba, dmap, structVarsPrim, 0);

        // compute all pairs
        // note: StructFactPrim option to compute only speicified pairs not written yet
        structFactPrim.define(ba,dmap,prim_var_names,var_scaling_prim);
        
        //////////////////////////////////////////////
        structFactConsMF.define(ba, dmap, structVarsCons, 0);

        // compute all pairs
        // note: StructFactCons option to compute only speicified pairs not written yet
        structFactCons.define(ba,dmap,cons_var_names,var_scaling_cons);
        
        //////////////////////////////////////////////
        
        // structure factor class for vertically-averaged dataset
        if(project_dir >= 0){
            // need a more elegant way of doing this
            MultiFab primVertAvg;  // flattened multifab defined below
            MultiFab prim_temp;
            prim_temp.define(ba,dmap,nprimvars,ngc);
            prim_temp.setVal(0.0);
            ComputeVerticalAverage(prim_temp, primVertAvg, geom, project_dir, 0, structVarsPrim);
            MultiFab primVertAvgRot = RotateFlattenedMF(primVertAvg);
            BoxArray ba_flat = primVertAvgRot.boxArray();
            const DistributionMapping& dmap_flat = primVertAvg.DistributionMap();
            {
                IntVect dom_lo(AMREX_D_DECL(0,0,0));
                IntVect dom_hi;
#if (AMREX_SPACEDIM == 2)
                if (project_dir == 0) {
                    dom_hi[0] = n_cells[1]-1;
                    dom_hi[1] = 0;
                }
                else if (project_dir == 1) {
                    dom_hi[0] = n_cells[0]-1;
                    dom_hi[1] = 0;
                }
#elif (AMREX_SPACEDIM == 3)
                if (project_dir == 0) {
                    dom_hi[0] = n_cells[1]-1;
                    dom_hi[1] = n_cells[2]-1;
                    dom_hi[2] = 0;
                } else if (project_dir == 1) {
                    dom_hi[0] = n_cells[0]-1;
                    dom_hi[1] = n_cells[2]-1;
                    dom_hi[2] = 0;
                } else if (project_dir == 2) {
                    dom_hi[0] = n_cells[0]-1;
                    dom_hi[1] = n_cells[1]-1;
                    dom_hi[2] = 0;
                }
#endif
                Box domain(dom_lo, dom_hi);

                // This defines the physical box
                Vector<Real> projected_hi(AMREX_SPACEDIM);
                for (int d=0; d<AMREX_SPACEDIM; d++) {
                    projected_hi[d] = prob_hi[d];
                }
#if (AMREX_SPACEDIM == 2)
                if (project_dir == 0) {
                    projected_hi[0] = prob_hi[1];
                }
#elif (AMREX_SPACEDIM == 3)
                if (project_dir == 0) {
                    projected_hi[0] = prob_hi[1];
                    projected_hi[1] = prob_hi[2];
                } else if (project_dir == 1) {
                    projected_hi[1] = prob_hi[2];
                }
#endif
        
                projected_hi[AMREX_SPACEDIM-1] = prob_hi[project_dir] / n_cells[project_dir];

                RealBox real_box({AMREX_D_DECL(     prob_lo[0],     prob_lo[1],     prob_lo[2])},
                                 {AMREX_D_DECL(projected_hi[0],projected_hi[1],projected_hi[2])});
          
                // This defines a Geometry object
                geom_flat.define(domain,&real_box,CoordSys::cartesian,is_periodic.data());
            }

            structFactPrimVerticalAverage.define(ba_flat,dmap_flat,prim_var_names,var_scaling_prim);
            structFactConsVerticalAverage.define(ba_flat,dmap_flat,cons_var_names,var_scaling_cons);
            //structFactPrimVerticalAverage.~StructFact(); // destruct
            //new(&structFactPrimVerticalAverage) StructFact(ba_flat,dmap_flat,prim_var_names,var_scaling); // reconstruct
    
        }

    }

    else {

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

        // transport properties
        eta.define(ba,dmap,1,ngc);
        zeta.define(ba,dmap,1,ngc);
        kappa.define(ba,dmap,1,ngc);
        chi.define(ba,dmap,nspecies,ngc);
        D.define(ba,dmap,nspecies*nspecies,ngc);

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
        cu.define(ba,dmap,nvars,ngc);

        for (int d=0; d<AMREX_SPACEDIM; d++) {
            vel[d].define(convert(ba,nodal_flag_dir[d]), dmap, 1, ngc);
            cumom[d].define(convert(ba,nodal_flag_dir[d]), dmap, 1, ngc);
        }

        //primative quantaties
        // in C++ indexing (add +1 for F90)
        // 0            (rho; density)
        // 1-3          (vel; velocity)
        // 4            (T;   temperature)
        // 5            (p;   pressure)
        // 6:6+ns-1     (Yk;  mass fractions)
        // 6+ns:6+2ns-1 (Xk;  mole fractions)
        prim.define(ba,dmap,nprimvars,ngc);

        cuMeans.define(ba,dmap,nvars,ngc);
        cuVars.define(ba,dmap,nvars,ngc);
        cuMeans.setVal(0.0);
        cuVars.setVal(0.0);
        
        primMeans.define(ba,dmap,nprimvars,ngc);
        primVars.define(ba,dmap,nprimvars+5,ngc);
        primMeans.setVal(0.0);
        primVars.setVal(0.0);

        // List of covariances (all cell centered)
        // 0: <rho jx>
        // 1: <rho jy>
        // 2: <rho jz>
        // 3: <jx jy>
        // 4: <jy jz>
        // 5: <jx jz>
        // 6: <rho rhoE>
        // 7: <rhoE jx>
        // 8: <rhoE jy>
        // 9: <rhoE jz>
        // 10: <rhoYklightest rhoYkheaviest>
        // 11: <rho vx>
        // 12: <rho vy>
        // 13: <rho vz>
        // 14: <vx vy>
        // 15: <vy vz>
        // 16: <vx vz>
        // 17: <rho T>
        // 18: <vx T>
        // 19: <vy T>
        // 20: <vz T>
        coVars.define(ba,dmap,21,0);
        coVars.setVal(0.0);

        for (int d=0; d<AMREX_SPACEDIM; d++) {
            velMeans[d].define(convert(ba,nodal_flag_dir[d]), dmap, 1, 0);
            cumomMeans[d].define(convert(ba,nodal_flag_dir[d]), dmap, 1, 0);
            velVars[d].define(convert(ba,nodal_flag_dir[d]), dmap, 1, 0);
            cumomVars[d].define(convert(ba,nodal_flag_dir[d]), dmap, 1, 0);
            velMeans[d].setVal(0.);
            velVars[d].setVal(0.);
            cumomMeans[d].setVal(0.);
            cumomVars[d].setVal(0.);
        }

        ///////////////////////////////////////////
        // Setup Structure factor
        ///////////////////////////////////////////

        structFactPrimMF.define(ba, dmap, structVarsPrim, 0);

        // compute all pairs
        // note: StructFactPrim option to compute only speicified pairs not written yet
        structFactPrim.define(ba,dmap,prim_var_names,var_scaling_prim);
        
        //////////////////////////////////////////////
        structFactConsMF.define(ba, dmap, structVarsCons, 0);

        // compute all pairs
        // note: StructFactCons option to compute only speicified pairs not written yet
        structFactCons.define(ba,dmap,cons_var_names,var_scaling_cons);
        
        //////////////////////////////////////////////
        
        // structure factor class for vertically-averaged dataset
        if(project_dir >= 0){
            MultiFab primVertAvg;  // flattened multifab defined below
            prim.setVal(0.0);
            ComputeVerticalAverage(prim, primVertAvg, geom, project_dir, 0, structVarsPrim);
            MultiFab primVertAvgRot = RotateFlattenedMF(primVertAvg);
            BoxArray ba_flat = primVertAvgRot.boxArray();
            const DistributionMapping& dmap_flat = primVertAvg.DistributionMap();
            {
                IntVect dom_lo(AMREX_D_DECL(0,0,0));
                IntVect dom_hi;
#if (AMREX_SPACEDIM == 2)
                if (project_dir == 0) {
                    dom_hi[0] = n_cells[1]-1;
                    dom_hi[1] = 0;
                }
                else if (project_dir == 1) {
                    dom_hi[0] = n_cells[0]-1;
                    dom_hi[1] = 0;
                }
#elif (AMREX_SPACEDIM == 3)
                if (project_dir == 0) {
                    dom_hi[0] = n_cells[1]-1;
                    dom_hi[1] = n_cells[2]-1;
                    dom_hi[2] = 0;
                } else if (project_dir == 1) {
                    dom_hi[0] = n_cells[0]-1;
                    dom_hi[1] = n_cells[2]-1;
                    dom_hi[2] = 0;
                } else if (project_dir == 2) {
                    dom_hi[0] = n_cells[0]-1;
                    dom_hi[1] = n_cells[1]-1;
                    dom_hi[2] = 0;
                }
#endif
                Box domain(dom_lo, dom_hi);

                // This defines the physical box
                Vector<Real> projected_hi(AMREX_SPACEDIM);
                for (int d=0; d<AMREX_SPACEDIM; d++) {
                    projected_hi[d] = prob_hi[d];
                }
#if (AMREX_SPACEDIM == 2)
                if (project_dir == 0) {
                    projected_hi[0] = prob_hi[1];
                }
#elif (AMREX_SPACEDIM == 3)
                if (project_dir == 0) {
                    projected_hi[0] = prob_hi[1];
                    projected_hi[1] = prob_hi[2];
                } else if (project_dir == 1) {
                    projected_hi[1] = prob_hi[2];
                }
#endif
        
                projected_hi[AMREX_SPACEDIM-1] = prob_hi[project_dir] / n_cells[project_dir];

                RealBox real_box({AMREX_D_DECL(     prob_lo[0],     prob_lo[1],     prob_lo[2])},
                                 {AMREX_D_DECL(projected_hi[0],projected_hi[1],projected_hi[2])});
          
                // This defines a Geometry object
                geom_flat.define(domain,&real_box,CoordSys::cartesian,is_periodic.data());
            }

            structFactPrimVerticalAverage.define(ba_flat,dmap_flat,prim_var_names,var_scaling_prim);
            structFactConsVerticalAverage.define(ba_flat,dmap_flat,cons_var_names,var_scaling_cons);
            //structFactPrimVerticalAverage.~StructFact(); // destruct
            //new(&structFactPrimVerticalAverage) StructFact(ba_flat,dmap_flat,prim_var_names,var_scaling); // reconstruct
    
        }

        ///////////////////////////////////////////
        // Initialize everything
        ///////////////////////////////////////////
        
        prim.setVal(0.0,0,nprimvars,ngc);
        prim.setVal(rho0,0,1,ngc);      // density
        prim.setVal(0.,1,3,ngc);        // x/y/z velocity
        prim.setVal(T_init[0],4,1,ngc); // temperature
                                        // pressure computed later in conservedToPrimitive
        for(int i=0;i<nspecies;i++) {
            prim.setVal(rhobar[i],6+i,1,ngc);    // mass fractions
        }

        //Initialize physical parameters from input vals
        double intEnergy, T0;
        T0 = T_init[0];

        // compute internal energy
        double massvec[nspecies];
        for(int i=0;i<nspecies;i++) {
            massvec[i] = rhobar[i];
        }
        get_energy(&intEnergy, massvec, &T0);

        cu.setVal(0.0,0,nvars,ngc);
        cu.setVal(rho0,0,1,ngc);           // density
        cu.setVal(0,1,3,ngc);              // x/y/z momentum
        cu.setVal(rho0*intEnergy,4,1,ngc); // total energy
        for(int i=0;i<nspecies;i++) {
            cu.setVal(rho0*rhobar[i],5+i,1,ngc); // mass densities
        }

        for (int d=0; d<AMREX_SPACEDIM; d++) { // staggered momentum & velocities
          cumom[d].setVal(0.);
          vel[d].setVal(0.);
        }

        // initialize conserved variables
        if (prob_type > 1) {
            InitConsVarStag(cu,cumom,prim,geom); // Need to add for staggered -- Ishan
        }

        // Write initial plotfile
        conservedToPrimitiveStag(prim, vel, cu, cumom);

        // Set BC: 1) fill boundary 2) physical (How to do for staggered? -- Ishan)
        cu.FillBoundary(geom.periodicity());
        prim.FillBoundary(geom.periodicity());
        for (int d=0; d<AMREX_SPACEDIM; d++) {
            cumom[d].FillBoundary(geom.periodicity());
            vel[d].FillBoundary(geom.periodicity());
        }

        setBCStag(prim, cu, cumom, vel, geom);
        
        if (plot_int > 0) {
          WritePlotFileStag(0, 0.0, geom, cu, cuMeans, cuVars, cumom, cumomMeans, cumomVars, 
                          prim, primMeans, primVars, vel, velMeans, velVars, coVars, eta, kappa);

          if (plot_cross) WriteSpatialCross(spatialCross, 0, dx);
        }

        if (plot_cross) {
            if (ParallelDescriptor::IOProcessor()) outfile.open(filename);
        }

        step_start = 1;
        time = 0.;
        statsCount = 1;

    } // end t=0 setup

    /////////////////////////////////////////////////
    // Initialize Fluxes and Sources
    /////////////////////////////////////////////////
    
    // external source term - possibly for later
    MultiFab source(ba,dmap,nprimvars,ngc);
    source.setVal(0.0);

    //fluxes (except momentum) at faces
    std::array< MultiFab, AMREX_SPACEDIM > faceflux;
    AMREX_D_TERM(faceflux[0].define(convert(ba,nodal_flag_x), dmap, nvars, 0);,
                 faceflux[1].define(convert(ba,nodal_flag_y), dmap, nvars, 0);,
                 faceflux[2].define(convert(ba,nodal_flag_z), dmap, nvars, 0););

    //momentum flux (edge + center)
#if (AMREX_SPACEDIM == 3)
    std::array< MultiFab, 2 > edgeflux_x; // divide by dx
    std::array< MultiFab, 2 > edgeflux_y; // divide by dy
    std::array< MultiFab, 2 > edgeflux_z; // divide by dz

    edgeflux_x[0].define(convert(ba,nodal_flag_xy), dmap, 1, 0); // 0-2: rhoU, rhoV, rhoW
    edgeflux_x[1].define(convert(ba,nodal_flag_xz), dmap, 1, 0);
                 
    edgeflux_y[0].define(convert(ba,nodal_flag_xy), dmap, 1, 0);
    edgeflux_y[1].define(convert(ba,nodal_flag_yz), dmap, 1, 0);

    edgeflux_z[0].define(convert(ba,nodal_flag_xz), dmap, 1, 0);
    edgeflux_z[1].define(convert(ba,nodal_flag_yz), dmap, 1, 0);

#elif (AMREX_SPACEDIM == 2)
    Abort("Currently requires AMREX_SPACEDIM=3");
#endif

    std::array< MultiFab, AMREX_SPACEDIM > cenflux;
    AMREX_D_TERM(cenflux[0].define(ba,dmap,1,1);, // 0-2: rhoU, rhoV, rhoW
                 cenflux[1].define(ba,dmap,1,1);,
                 cenflux[2].define(ba,dmap,1,1););
                
    /////////////////////////////////////////////////
    //Time stepping loop
    /////////////////////////////////////////////////
    
    for (int step=step_start;step<=max_step;++step) {

        // timer
        Real ts1 = ParallelDescriptor::second();
    
        RK3stepStag(cu, cumom, prim, vel, source, eta, zeta, kappa, chi, D, 
            faceflux, edgeflux_x, edgeflux_y, edgeflux_z, cenflux, geom, dt, step);

        // timer
        Real ts2 = ParallelDescriptor::second() - ts1;
        ParallelDescriptor::ReduceRealMax(ts2);
        amrex::Print() << "Advanced step " << step << " in " << ts2 << " seconds\n";

        // timer
        Real aux1 = ParallelDescriptor::second();
        
        // reset statistics after n_steps_skip
        // if n_steps_skip is negative, we use it as an interval
        if ((n_steps_skip > 0 && step == n_steps_skip) ||
            (n_steps_skip < 0 && step%n_steps_skip == 0) ) {

            cuMeans.setVal(0.0);
            cuVars.setVal(0.0);
            primMeans.setVal(0.0);
            primVars.setVal(0.0);

            for (int d=0; d<AMREX_SPACEDIM; d++) {
                velMeans[d].setVal(0.0);
                velVars[d].setVal(0.0);
                cumomMeans[d].setVal(0.0);
                cumomVars[d].setVal(0.0);
            }

            coVars.setVal(0.0);
            spatialCross.assign(spatialCross.size(), 0.0);

            std::printf("Resetting stat collection.\n");

            statsCount = 1;

        }

        // Evaluate Statistics
        evaluateStatsStag(cu, cuMeans, cuVars, prim, primMeans, primVars, vel, 
                          velMeans, velVars, cumom, cumomMeans, cumomVars, coVars,
                          yzAvMeans_cross, spatialCross, statsCount, dx);
        statsCount++;

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
             //yzAverage(cuMeans, cuVars, primMeans, primVars, spatialCross,
             //          cuMeansAv, cuVarsAv, primMeansAv, primVarsAv, spatialCrossAv);
            WritePlotFileStag(step, time, geom, cu, cuMeans, cuVars, cumom, cumomMeans, cumomVars,
                          prim, primMeans, primVars, vel, velMeans, velVars, coVars, eta, kappa);

            if (plot_cross) {
                WriteSpatialCross(spatialCross, step, dx);
                if (ParallelDescriptor::IOProcessor()) {
                    outfile << step << " ";
                    for (auto l=0; l<2*nvars+8; ++l) {
                        outfile << yzAvMeans_cross[l] << " ";
                    }
                    outfile << std::endl;
                }
            }
        }

        // collect a snapshot for structure factor
        if (step > amrex::Math::abs(n_steps_skip) && 
            struct_fact_int > 0 && 
            (step-amrex::Math::abs(n_steps_skip))%struct_fact_int == 0) {

            MultiFab::Copy(structFactPrimMF, prim, 0, 0, structVarsPrim-AMREX_SPACEDIM  , 0);
            MultiFab::Copy(structFactConsMF, cu,   0, 0, structVarsCons-AMREX_SPACEDIM-1, 0);
            // temperature too
            MultiFab::Copy(structFactConsMF, prim, AMREX_SPACEDIM+1, structVarsCons-1-AMREX_SPACEDIM, 1, 0);

            // append the shifted momentum and velocities to the end
            for (int d=0; d<AMREX_SPACEDIM; ++d) {
                ShiftFaceToCC(vel[d]  ,0,structFactPrimMF,d+structVarsPrim-AMREX_SPACEDIM,1);
                ShiftFaceToCC(cumom[d],0,structFactConsMF,d+structVarsCons-AMREX_SPACEDIM,1);
            }
            
            structFactPrim.FortStructure(structFactPrimMF,geom);
            structFactCons.FortStructure(structFactConsMF,geom);

            if(project_dir >= 0) {
                MultiFab primVertAvg;  // flattened multifab defined below
                MultiFab consVertAvg;  // flattened multifab defined below
                ComputeVerticalAverage(structFactPrimMF, primVertAvg, geom, project_dir, 0, structVarsPrim);
                ComputeVerticalAverage(structFactConsMF, consVertAvg, geom, project_dir, 0, structVarsCons);
                MultiFab primVertAvgRot = RotateFlattenedMF(primVertAvg);
                MultiFab consVertAvgRot = RotateFlattenedMF(consVertAvg);
                structFactPrimVerticalAverage.FortStructure(primVertAvg,geom_flat);
                structFactConsVerticalAverage.FortStructure(consVertAvg,geom_flat);
            }
        }

        // write out structure factor
        if (step > amrex::Math::abs(n_steps_skip) && 
            struct_fact_int > 0 && plot_int > 0 && 
            step%plot_int == 0) {

            structFactPrim.WritePlotFile(step,time,geom,"plt_SF_prim");
            structFactCons.WritePlotFile(step,time,geom,"plt_SF_cons");
            if(project_dir >= 0) {
                structFactPrimVerticalAverage.WritePlotFile(step,time,geom_flat,"plt_SF_prim_VerticalAverage");
                structFactConsVerticalAverage.WritePlotFile(step,time,geom_flat,"plt_SF_cons_VerticalAverage");
            }
        }
        

        // write checkpoint file
        if (chk_int > 0 && step > 0 && step%chk_int == 0)
        {
            WriteCheckPoint(step, time, statsCount, geom, cu, cuMeans, cuVars, prim,
                           primMeans, primVars, cumom, cumomMeans, cumomVars, 
                           vel, velMeans, velVars, coVars, spatialCross);
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

    if (ParallelDescriptor::IOProcessor()) outfile.close();

    // timer
    Real stop_time = ParallelDescriptor::second() - strt_time;
    ParallelDescriptor::ReduceRealMax(stop_time);
    amrex::Print() << "Run time = " << stop_time << std::endl;
}
