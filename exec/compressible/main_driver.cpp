
#include "rng_functions.H"
#include "rng_functions_F.H"

#include "common_functions.H"
#include "common_functions_F.H"

#include "common_namespace.H"
#include "common_namespace_declarations.H"

#include "compressible_functions.H"
#include "compressible_functions_F.H"

#include "exec_functions.H"

#include "StructFact.H"

#include <AMReX_VisMF.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_MultiFab.H>

#include "compressible_test_functions_F.H"

using namespace amrex;
using namespace common;

// argv contains the name of the inputs file entered at the command line
void main_driver(const char* argv)
{

    // store the current time so we can later compute total run time.
    Real strt_time = ParallelDescriptor::second();

    std::string inputs_file = argv;

    // read in parameters from inputs file into F90 modules
    // we use "+1" because of amrex_string_c_to_f expects a null char termination
    read_common_namelist(inputs_file.c_str(),inputs_file.size()+1);
    // read_gmres_namelist(inputs_file.c_str(),inputs_file.size()+1);

    // copy contents of F90 modules to C++ namespaces
    InitializeCommonNamespace();
    // InitializeGmresNamespace();

    //if gas heat capacities are negative, calculate using dofs. This will only update the Fortran values.
    get_hc_gas();
  
    // is the problem periodic?
    Vector<int> is_periodic(AMREX_SPACEDIM,0);  // set to 0 (not periodic) by default
    for (int i=0; i<AMREX_SPACEDIM; ++i) {
        if (bc_lo[i] == -1 && bc_hi[i] == -1) {
            is_periodic[i] = 1;
            Print() << "Periodic: " << is_periodic[i] << "\n";
        }
        //is_periodic[i] = 0;
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

    // Print() << "Hack: boxarray = " << ba << "\n";

    Real dt = fixed_dt;
    Real dtinv = 1.0/dt;
    const Real* dx = geom.CellSize();
    const RealBox& realDomain = geom.ProbDomain();

    // how boxes are distrubuted among MPI processes
    DistributionMapping dmap(ba);

    /////////////////////////////////////////
    //Initialise rngs
    /////////////////////////////////////////
    const int n_rngs = 1;

    const int proc = ParallelDescriptor::MyProc();

    int fhdSeed = 0;
    int particleSeed = 2;
    int selectorSeed = 3;
    int thetaSeed = 4;
    int phiSeed = 5;
    int generalSeed = 0;

    //fhdSeed += 10000*proc;
    particleSeed += 20000*proc;
    selectorSeed += 30000*proc;
    thetaSeed += 40000*proc;
    phiSeed += 50000*proc;
    //generalSeed += 60000*proc;

    //Initialise rngs
    rng_initialize(&fhdSeed,&particleSeed,&selectorSeed,&thetaSeed,&phiSeed,&generalSeed);
    /////////////////////////////////////////

    Real Cfl_max[4];
    int eval_dt = 0;

    //transport properties
    MultiFab eta  (ba,dmap,1,ngc);
    MultiFab zeta (ba,dmap,1,ngc);
    MultiFab kappa(ba,dmap,1,ngc);
    MultiFab chi(ba,dmap,nspecies,ngc);
    MultiFab D(ba,dmap,nspecies*nspecies,ngc);

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

    MultiFab etaMean(ba,dmap,1,ngc);
    MultiFab kappaMean(ba,dmap,1,ngc);

    MultiFab etaMeanAv(ba,dmap,1,ngc);
    MultiFab kappaMeanAv(ba,dmap,1,ngc);


    Real delHolder1[n_cells[1]*n_cells[2]];
    Real delHolder2[n_cells[1]*n_cells[2]];
    Real delHolder3[n_cells[1]*n_cells[2]];
    Real delHolder4[n_cells[1]*n_cells[2]];
    Real delHolder5[n_cells[1]*n_cells[2]];
    Real delHolder6[n_cells[1]*n_cells[2]];

    MultiFab spatialCross(ba,dmap,6,ngc);

    MultiFab spatialCrossAv(ba,dmap,6,ngc);

    cuMeans.setVal(0.0);
    cuVars.setVal(0.0);

    primMeans.setVal(0.0);
    primVars.setVal(0.0);

    etaMean.setVal(0.0);
    kappaMean.setVal(0.0);

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

    // amrex::Print() << "Hack: T_init = " << T_init[0] << "\n";

    double massvec[nspecies];
    double intEnergy, T0;

    T0 = T_init[0];

    for(int i=0;i<nspecies;i++)
    {
        prim.setVal(rhobar[i],6+i,1,ngc);
        cu.setVal(rho0*rhobar[i],5+i,1,ngc);

        massvec[i] = rhobar[i];
    }

    get_energy(&intEnergy, massvec, &T0);

    cu.setVal(rho0,0,1,ngc);
    cu.setVal(0,1,1,ngc);
    cu.setVal(0,2,1,ngc);
    cu.setVal(0,3,1,ngc);
    cu.setVal(rho0*intEnergy,4,1,ngc);

    cup.setVal(rho0,0,1,ngc);
    cup2.setVal(rho0,0,1,ngc);
    cup3.setVal(rho0,0,1,ngc);
    
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
    // structure factor:
    ///////////////////////////////////////////
    
    Real molmix;
    Real P_bar;
    Real dVol = dx[0]*dx[1];
    Real c_v;
    int tot_n_cells = n_cells[0]*n_cells[1];
    // calc eqm pressure
    molmix = 0.0;
    for(int i=0; i<nspecies; i++) {
	molmix += rhobar[i]/molmass[i];
    }
    molmix = 1.0/molmix;
    P_bar = rho0*(Runiv/molmix)*T_init[0];
    c_v = 1.5*Runiv/molmix;
    // calc cell volume
    if (AMREX_SPACEDIM == 2) {
	dVol *= cell_depth;
    } else if (AMREX_SPACEDIM == 3) {
	dVol *= dx[2];
	tot_n_cells = n_cells[2]*tot_n_cells;
    }
    // calc momentum variance
    Real Jeqvar = rho0*k_B*T_init[0]/dVol;

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
    beqmvars[cnt++] = rho0/P_bar*Jeqvar;                      // rho,rho      - not scaled by mass fracs
    beqmvars[cnt++] = sqrt(rho0/P_bar)*Jeqvar;                // J_i,rho      - dimensionless only
    beqmvars[cnt++] = sqrt(rho0/P_bar*T_init[0]*c_v)*Jeqvar;  // rhoE,rho     - dimensionless only
    beqmvars[cnt++] = rho0/P_bar*Jeqvar;                      // rho_k,rho    - not scaled by mass fracs
    beqmvars[cnt++] = Jeqvar;                                 // J_i,J_j
    beqmvars[cnt++] = sqrt(T_init[0]*c_v)*Jeqvar;             // rhoE,J_j     - dimensionless only
    beqmvars[cnt++] = sqrt(rho0/P_bar)*Jeqvar;                // rho_k,J_j    - dimensionless only
    beqmvars[cnt++] = T_init[0]*c_v*Jeqvar;                   // rhoE,rhoE
    beqmvars[cnt++] = sqrt(rho0/P_bar*T_init[0]*c_v)*Jeqvar;  // rho_k,rhoE   - not scaled by mass fracs
    beqmvars[cnt++] = rho0/P_bar*Jeqvar;                      // rho_k,rho_l  - not scaled by mass fracs

    // for (int d=0; d<nbcov_sf; d++) {
    //   beqmvars[d] = 1.0;
    // }

    cnt = 0;
    int bcnt = 0;
    int ig, jg = 0;
    // loop over matrix blocks
    for(int jb=0; jb<nb_sf; jb++) {
      ig = jg;
      for(int ib=jb; ib<nb_sf; ib++) {
	// loop within blocks
	for(int j=0; j<blocks[jb]; j++) {
	  int low_ind;
	  if(ib==jb){      // if block lies on diagonal...
	    low_ind=j;
	  } else {
	    low_ind=0;
	  }
	  for(int i=low_ind; i<blocks[ib]; i++) {
	    cnt = (ig+i)+nvar_sf*(jg+j)-(jg+j)*(jg+j+1)/2;
	    eqmvars[cnt] = beqmvars[bcnt];
	    // Print() << "Hack: " << ig+i << " " << jg+j << " " << cnt << " " << eqmvars[cnt] << std::endl;
	  }
	}
	bcnt++;
	ig += blocks[ib];
      }
      jg += blocks[jb];
    }
    
    // for (int d=0; d<ncov_sf; d++) {
    //   Print() << "Hack: scaling = " << eqmvars[d] << " " << d << std::endl;
    // }
    // exit(0);

    // set variable names
    cnt = 0;
    Vector< std::string > var_names;
    var_names.resize(nvar_sf);
    std::string x;
    var_names[cnt++] = "rho";
    for (int d=0; d<AMREX_SPACEDIM; d++) {
      x = "mom";
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

    StructFact structFact(ba,dmap,var_names,eqmvars);
    // StructFact structFact(ba,dmap,var_names,eqmvars,s_pairA,s_pairB);

    ///////////////////////////////////////////

    //Initialize everything
    
    // chi.setVal(0.0);
    // D.setVal(0.0);
    // calculateTransportCoeffs(prim, eta, zeta, kappa, chi, D);

    // conservedToPrimitive(prim, cu);
    // cu.FillBoundary(geom.periodicity());
    // prim.FillBoundary(geom.periodicity());

    // calculateTransportCoeffs(prim, eta, zeta, kappa, chi, D);

    // eta.FillBoundary(geom.periodicity());
    // zeta.FillBoundary(geom.periodicity());
    // kappa.FillBoundary(geom.periodicity());
    // chi.FillBoundary(geom.periodicity());
    // D.FillBoundary(geom.periodicity());

    // setBC(prim, cu, eta, zeta, kappa);

    // calculateFlux(cu, prim, eta, zeta, kappa, flux, stochFlux, cornx, corny, cornz, visccorn, rancorn, geom, dx, dt);

    statsCount = 1;

    //Time stepping loop
    for(step=1;step<=max_step;++step)
    {

        RK3step(cu, cup, cup2, cup3, prim, source, eta, zeta, kappa, chi, D, flux, stochFlux, cornx, corny, cornz, visccorn, rancorn, geom, dx, dt);

        if(step == n_steps_skip)
        {
            cuMeans.setVal(0.0);
            cuVars.setVal(0.0);

            primMeans.setVal(0.0);
            primVars.setVal(0.0);

            spatialCross.setVal(0.0);

            statsCount = 1;

            //dt = 2.0*dt;
        }

//        if(step == (int)floor((double)n_steps_skip/2.0))
//        {
//            cuMeans.setVal(0.0);
//            cuVars.setVal(0.0);

//            primMeans.setVal(0.0);
//            primVars.setVal(0.0);

//            spatialCross.setVal(0.0);

//            statsCount = 1;

//            dt = 2.0*dt;
//        }

	// for ( MFIter mfi(cu); mfi.isValid(); ++mfi ) {
	//   const Box& bx = mfi.validbox();

	//   init_consvar(BL_TO_FORTRAN_BOX(bx),
	// 	       BL_TO_FORTRAN_FAB(cu[mfi]),
	// 	       dx, 
	// 	       // geom.ProbLo(), geom.ProbHi(),
	// 	       ZFILL(realDomain.lo()), ZFILL(realDomain.hi()));

	// }
	
	if (step > n_steps_skip) {
	  // evaluateStats(cu, cuMeans, cuVars, prim, primMeans, primVars, spatialCross, eta, etaMean, kappa, kappaMean, delHolder1, delHolder2, delHolder3, delHolder4, delHolder5, delHolder6, statsCount, dx);
	}

	///////////////////////////////////////////
	// Update structure factor
	///////////////////////////////////////////
	if (step > n_steps_skip && struct_fact_int > 0 && (step-n_steps_skip-1)%struct_fact_int == 0) {
	  MultiFab::Copy(struct_in_cc, cu, 0, 0, nvar_sf, 0);
	  structFact.FortStructure(struct_in_cc,geom);
        }
	///////////////////////////////////////////

        statsCount++;

        // if(step%100 == 0)
        // {    
	amrex::Print() << "Advanced step " << step << "\n";
        // }

        if (plot_int > 0 && step > 0 && step%plot_int == 0)
        {

           // yzAverage(cuMeans, cuVars, primMeans, primVars, spatialCross, etaMean, kappaMean, cuMeansAv, cuVarsAv, primMeansAv, primVarsAv, spatialCrossAv, etaMeanAv, kappaMeanAv);
           // WritePlotFile(step, time, geom, cu, cuMeansAv, cuVarsAv, prim, primMeansAv, primVarsAv, spatialCrossAv, etaMeanAv, kappaMeanAv);

           WritePlotFile(step, time, geom, cu, cuMeans, cuVars, prim, primMeans, primVars, spatialCross, eta, kappa);
	   if (struct_fact_int > 0) {
	       structFact.WritePlotFile(step,time,geom);
	   }
        }

        time = time + dt;
    }


    Real stop_time = ParallelDescriptor::second() - strt_time;
    ParallelDescriptor::ReduceRealMax(stop_time);
    amrex::Print() << "Run time = " << stop_time << std::endl;

    // amrex::Abort();
    // amrex::Finalize();

}

