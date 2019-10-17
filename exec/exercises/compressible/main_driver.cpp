#include "common_functions.H"
#include "common_functions_F.H"
#include "common_namespace.H"
#include "common_namespace_declarations.H"

#include "compressible_functions.H"
#include "compressible_functions_F.H"
#include "compressible_namespace.H"
#include "compressible_namespace_declarations.H"

#include "rng_functions_F.H"

#include "StructFact.H"

using namespace amrex;
using namespace common;
using namespace compressible;

// argv contains the name of the inputs file entered at the command line
void main_driver(const char* argv)
{
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

    // if gas heat capacities in the namelist are negative, calculate them using using dofs.
    // This will only update the Fortran values.
    get_hc_gas();
  
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
    for (int i=0; i<AMREX_SPACEDIM; ++i) {
        if (bc_vel_lo[i] == -1) {
            bc_mass_lo[i] = -1;
            bc_mass_hi[i] = -1;
            bc_therm_lo[i] = -1;
            bc_therm_hi[i] = -1;
        }
    }
    setup_bc(); // do the same in the fortran namelist

    // if multispecies
    if (algorithm_type == 2) {
        // compute wall concentrations if BCs call for it
        setup_cwall();
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

    // NOTE: only fhdSeed is used currently
        // zero is a clock-based seed
    int fhdSeed      = 0;
    int particleSeed = 0;
    int selectorSeed = 0;
    int thetaSeed    = 0;
    int phiSeed      = 0;
    int generalSeed  = 0;

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
    cuMeans.setVal(0.0);
    cuVars.setVal(0.0);
    
    MultiFab cuVertAvg;  // flattened multifab defined below

    MultiFab primMeans  (ba,dmap,nprimvars  ,ngc);
    MultiFab primVars   (ba,dmap,nprimvars+5,ngc);
    primMeans.setVal(0.0);
    primVars.setVal(0.0);
    
    // external source term - possibly for later
    MultiFab source(ba,dmap,nprimvars,ngc);
    source.setVal(0.0);

    //Initialize physical parameters from input vals

    double intEnergy, T0;

    T0 = T_init[0];

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
    
    ////////////////////////////////
    // create equilibrium covariance matrix
    Real molmix, avgmolmass;
    Real P_bar;
    Real c_v, c_v2;
    Real dVol = dx[0]*dx[1];
    molmix = 0.0;
    avgmolmass = 0.0;
    for(int i=0; i<nspecies; i++) {
      molmix     += rhobar[i]/molmass[i];
      avgmolmass += rhobar[i]*molmass[i];
    }
    molmix = 1.0/molmix;                // molar mass of mixture
    P_bar = rho0*(Runiv/molmix)*T0;     // eqm pressure
    c_v  = 1.5*(Runiv/molmix);          // Assuming all gases are monoatomic (dof=3)
    c_v2 = c_v*c_v;
    // calc cell volume
    if (AMREX_SPACEDIM == 2) {
	dVol *= cell_depth;
    } else if (AMREX_SPACEDIM == 3) {
	dVol *= dx[2];
    }
    // calc momentum variance
    Real Jeqmvar = rho0*k_B*T0/dVol;
    Real Meqmvar = (rho0/P_bar)*Jeqmvar;
    Real Eeqmvar = (c_v2*T0*T0)*Meqmvar + (c_v*T0)*Jeqmvar;
    Real MEeqmcovar = (rho0/P_bar)*(c_v*T0)*Jeqmvar;

    // setup covariance eqm. variance matrix
    int cnt = 0;
    int nvar_sf = nvars;
    int ncov_sf = nvar_sf*(nvar_sf+1)/2;
    int nb_sf = 4;
    int nbcov_sf = nb_sf*(nb_sf+1)/2;

    // struct. fact. eqm. variances
    Vector< Real > eqmvars(ncov_sf);
    Vector< Real > beqmvars(nbcov_sf);
    Vector< int > blocks(nb_sf);

    // layout of covariance matrix:
    //
    // <cons_i,cons_j> =
    // | <rho,rho>         ...      ...       ...            ...      ... |
    // | <Jx,rho>       <Jx,Jx>     ...       ...            ...      ... |
    // |    ...            ...      ...       ...            ...      ... |
    // | <rhoE,rho>     <rhoE,Jx>   ...   <rhoE,rhoE>        ...      ... |
    // | <rho1,rho>     <rho1,Jx>   ...   <rho1,rhoE>    <rho1,rho1>  ... |
    // |    ...            ...      ...       ...            ...      ... |

    // covariance matrix block sizes
    cnt = 0;
    blocks[cnt++] = 1;
    blocks[cnt++] = AMREX_SPACEDIM;
    blocks[cnt++] = 1;
    blocks[cnt++] = nspecies;

    // covariance matrix block scaling
    cnt = 0;
    beqmvars[cnt++] = Meqmvar*avgmolmass/molmix;  // rho,rho
    beqmvars[cnt++] = sqrt(Meqmvar*Jeqmvar);      // J_i,rho      - cheat scale
    beqmvars[cnt++] = MEeqmcovar;                 // rhoE,rho
    beqmvars[cnt++] = Meqmvar;                    // rho_k,rho    - scaled by mass fracs later
    beqmvars[cnt++] = Jeqmvar;                    // J_i,J_j
    beqmvars[cnt++] = sqrt(Eeqmvar*Jeqmvar);      // rhoE,J_j     - cheat scale
    beqmvars[cnt++] = sqrt(Meqmvar*Jeqmvar);      // rho_k,J_j    - cheat scale
    beqmvars[cnt++] = Eeqmvar;                    // rhoE,rhoE
    beqmvars[cnt++] = MEeqmcovar;                 // rho_k,rhoE   - scaled by mass fracs later
    beqmvars[cnt++] = Meqmvar;                    // rho_k,rho_l  - scaled by mass fracs later
    
//    // loop over lower triangular block matrix
//    cnt = 0;
//    int bcnt = 0;
//    int ig, jg = 0;
//    // loop over matrix blocks
//    for(int jb=0; jb<nb_sf; jb++) {
//      ig = jg;
//      for(int ib=jb; ib<nb_sf; ib++) {
//	// loop within blocks
//	for(int j=0; j<blocks[jb]; j++) {
//	  int low_ind;
//	  if(ib==jb){      // if block lies on diagonal...
//	    low_ind=j;
//	  } else {
//	    low_ind=0;
//	  }
//	  for(int i=low_ind; i<blocks[ib]; i++) {
//	    cnt = (ig+i)+nvar_sf*(jg+j)-(jg+j)*(jg+j+1)/2;
//	    eqmvars[cnt] = beqmvars[bcnt];

//	    // fix scale for individual species
//	    if(ib != 2 && jb != 2) { // not for energy

//	      if(blocks[ib]==nspecies && blocks[jb]==nspecies) {
//	    	eqmvars[cnt] *= sqrt(rhobar[i]*molmass[i]/molmix);
//	    	eqmvars[cnt] *= sqrt(rhobar[j]*molmass[j]/molmix);
//	      } else if (blocks[ib]==nspecies) {
//	    	eqmvars[cnt] *= rhobar[i]*molmass[i]/molmix;
//	      } else if (blocks[jb]==nspecies) {
//	    	eqmvars[cnt] *= rhobar[j]*molmass[j]/molmix;
//	      }

//	    } else {
//	      
//	      // if rho_k & energy, only scale by Yk
//	      if(blocks[ib]==nspecies && jb==2) {
//	    	eqmvars[cnt] *= rhobar[i];
//	      }
//	      
//	    }

//	  }
//	}
//	bcnt++;
//	ig += blocks[ib];
//      }
//      jg += blocks[jb];
//    }

    ////////////////////////////////

    // set variable names
    cnt = 0;
    Vector< std::string > var_names;
    var_names.resize(nvar_sf);
    std::string x;
    var_names[cnt++] = "rho";
    for (int d=0; d<AMREX_SPACEDIM; d++) {
      x = "J";
      x += (120+d);
      var_names[cnt++] = x;
    }
    var_names[cnt++] = "rhoE";
    for (int d=0; d<nspecies; d++) {
      x = "rho";
      x += (49+d);
      var_names[cnt++] = x;
    }

    MultiFab struct_in_cc;
    struct_in_cc.define(ba, dmap, nvar_sf, 0);

    amrex::Vector< int > s_pairA(nvar_sf);
    amrex::Vector< int > s_pairB(nvar_sf);

    // Select which variable pairs to include in structure factor:
    for (int d=0; d<nvar_sf; d++) {
      s_pairA[d] = d;
      s_pairB[d] = d;
    }

    // structure factor class for full dataset
    StructFact structFact(ba,dmap,var_names,eqmvars);
    
    // structure factor class for vertically-averaged dataset
    StructFact structFactVA;
    
    Geometry geom_flat;

    if(project_dir >= 0){
      cu.setVal(0.0);
      ComputeVerticalAverage(cu, cuVertAvg, geom, project_dir, 0, nvars);
      BoxArray ba_flat = cuVertAvg.boxArray();
      const DistributionMapping& dmap_flat = cuVertAvg.DistributionMap();
      {
        IntVect dom_lo(AMREX_D_DECL(           0,            0,            0));
        IntVect dom_hi(AMREX_D_DECL(n_cells[0]-1, n_cells[1]-1, n_cells[2]-1));
    	dom_hi[project_dir] = 0;
        Box domain(dom_lo, dom_hi);
	
    	// This defines the physical box
    	Vector<Real> projected_hi(AMREX_SPACEDIM);
    	for (int d=0; d<AMREX_SPACEDIM; d++) {
    	  projected_hi[d] = prob_hi[d];
    	}
    	projected_hi[project_dir] = prob_hi[project_dir]/n_cells[project_dir];
        RealBox real_box({AMREX_D_DECL(     prob_lo[0],     prob_lo[1],     prob_lo[2])},
                         {AMREX_D_DECL(projected_hi[0],projected_hi[1],projected_hi[2])});

        // This defines a Geometry object
        geom_flat.define(domain,&real_box,CoordSys::cartesian,is_periodic.data());
      }

      structFactVA.~StructFact(); // destruct
      new(&structFactVA) StructFact(ba_flat,dmap_flat,var_names,eqmvars); // reconstruct
    
    }

    ///////////////////////////////////////////

    // Initialize everything
    
    prim.setVal(0.0,0,nprimvars,ngc);
    prim.setVal(rho0,0,1,ngc);      // density
    prim.setVal(0.,1,3,ngc);        // x/y/z velocity
    prim.setVal(T_init[0],4,1,ngc); // temperature
                                    // pressure computed later in cons_to_prim
    for(int i=0;i<nspecies;i++) {
        prim.setVal(rhobar[i],6+i,1,ngc);    // mass fractions
    }

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
        for ( MFIter mfi(cu); mfi.isValid(); ++mfi ) {
            const Box& bx = mfi.validbox();
            init_consvar(BL_TO_FORTRAN_BOX(bx),
                         BL_TO_FORTRAN_ANYD(cu[mfi]),
                         BL_TO_FORTRAN_ANYD(prim[mfi]),
                         dx, ZFILL(realDomain.lo()), ZFILL(realDomain.hi()));
        }
    }

    statsCount = 1;

    // Write initial plotfile
    conservedToPrimitive(prim, cu);
    if (plot_int > 0) {
	WritePlotFile(0, 0.0, geom, cu, cuMeans, cuVars,
                      prim, primMeans, primVars, eta, kappa);
    }

    //Time stepping loop
    for(step=1;step<=max_step;++step) {

        // timer
        Real ts1 = ParallelDescriptor::second();
    
        RK3step(cu, cup, cup2, cup3, prim, source, eta, zeta, kappa, chi, D, flux,
                stochFlux, cornx, corny, cornz, visccorn, rancorn, geom, dx, dt);

        // timer
        Real ts2 = ParallelDescriptor::second() - ts1;
        ParallelDescriptor::ReduceRealMax(ts2);
    	amrex::Print() << "Advanced step " << step << " in " << ts2 << " seconds\n";

        // timer
        Real aux1 = ParallelDescriptor::second();
        
        // compute mean and variances
	if (step > n_steps_skip) {
            evaluateStats(cu, cuMeans, cuVars, prim, primMeans, primVars, statsCount, dx);
            statsCount++;
	}

        // write a plotfile
        if (plot_int > 0 && step > 0 && step%plot_int == 0) {
            WritePlotFile(step, time, geom, cu, cuMeans, cuVars,
                          prim, primMeans, primVars, eta, kappa);
        }
 
	// collect a snapshot for structure factor
	if (step > n_steps_skip && struct_fact_int > 0 && (step-n_steps_skip)%struct_fact_int == 0) {
            MultiFab::Copy(struct_in_cc, cu, 0, 0, nvar_sf, 0);
            structFact.FortStructure(struct_in_cc,geom);
            if(project_dir >= 0) {
                ComputeVerticalAverage(cu, cuVertAvg, geom, project_dir, 0, nvars);
                structFactVA.FortStructure(cuVertAvg,geom_flat);
            }
        }

        // write out structure factor
        if (step > n_steps_skip && struct_fact_int > 0 && plot_int > 0 && step%plot_int == 0) {
            structFact.WritePlotFile(step,time,geom,"plt_SF");
            if(project_dir >= 0) {
                structFactVA.WritePlotFile(step,time,geom_flat,"plt_SF_VA");
            }
        }
        
        // timer
        Real aux2 = ParallelDescriptor::second() - aux1;
        ParallelDescriptor::ReduceRealMax(aux2);
        amrex::Print() << "Aux time (stats, struct fac, plotfiles) " << aux2 << " seconds\n";
        
        time = time + dt;
    }

    // timer
    Real stop_time = ParallelDescriptor::second() - strt_time;
    ParallelDescriptor::ReduceRealMax(stop_time);
    amrex::Print() << "Run time = " << stop_time << std::endl;
}

