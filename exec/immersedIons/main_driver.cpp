

#include "common_functions.H"
#include "gmres_functions.H"

#include "common_functions_F.H"
#include "common_namespace.H"
#include "common_namespace_declarations.H"

#include "gmres_functions_F.H"
#include "gmres_namespace.H"
#include "gmres_namespace_declarations.H"

#include "INS_functions.H"
#include "rng_functions_F.H"

#include "species.H"
#include "surfaces.H"

#include "analysis_functions_F.H"
#include "StochMFlux.H"
#include "StructFact.H"

#include "hydro_test_functions_F.H"

#include "hydro_functions.H"
#include "hydro_functions_F.H"

#include "electrostatic.H"

using namespace gmres;
using namespace common;
//using namespace amrex;

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

    const int n_rngs = 1;

    int fhdSeed = ParallelDescriptor::MyProc() + 1;
    int particleSeed = 2*ParallelDescriptor::MyProc() + 2;
    int selectorSeed = 3*ParallelDescriptor::MyProc() + 3;
    int thetaSeed = 4*ParallelDescriptor::MyProc() + 4;
    int phiSeed = 5*ParallelDescriptor::MyProc() + 5;
    int generalSeed = 6*ParallelDescriptor::MyProc() + 6;

//    int fhdSeed = 0;
//    int particleSeed = 0;
//    int selectorSeed = 0;
//    int thetaSeed = 0;
//    int phiSeed = 0;
//    int generalSeed = 0;

    //Initialise rngs
    rng_initialize(&fhdSeed,&particleSeed,&selectorSeed,&thetaSeed,&phiSeed,&generalSeed);

    // is the problem periodic?
    Vector<int> is_periodic(AMREX_SPACEDIM,0);  // set to 0 (not periodic) by default
    for (int i=0; i<AMREX_SPACEDIM; ++i) {
        if (bc_lo[i] == -1 && bc_hi[i] == -1) {
            is_periodic[i] = 1;
        }
    }

    //---------------------------------------

        //Terms prepended with a 'C' are related to the particle grid. Those with P are for the electostatic grid. Those without are for the hydro grid.
        //The particle grid created as a corsening or refinement of the hydro grid.

    //---------------------------------------

    // make BoxArray and Geometry
    BoxArray ba;
    BoxArray bc;
    BoxArray bp;
    Geometry geom;
    Geometry geomC;
    Geometry geomP;
    
    IntVect dom_lo(AMREX_D_DECL(           0,            0,            0));
    IntVect dom_hi(AMREX_D_DECL(n_cells[0]-1, n_cells[1]-1, n_cells[2]-1));
    Box domain(dom_lo, dom_hi);
    Box domainC = domain;
    Box domainP = domain;

    // Initialize the boxarray "ba" from the single box "bx"
    ba.define(domain);

    // Break up boxarray "ba" into chunks no larger than "max_grid_size" along a direction
    // note we are converting "Vector<int> max_grid_size" to an IntVect
    ba.maxSize(IntVect(max_grid_size));

    RealBox real_box({AMREX_D_DECL(prob_lo[0],prob_lo[1],prob_lo[2])},
                     {AMREX_D_DECL(prob_hi[0],prob_hi[1],prob_hi[2])});

    //This must be an even number for now?
    int sizeRatioC = 1;
    int sizeRatioP = 1;

    bc = ba;
    bp = ba;
    bc.coarsen(sizeRatioC);
    bp.coarsen(sizeRatioP);
    domainC.coarsen(sizeRatioC);
    domainP.coarsen(sizeRatioP);
    // This defines a Geometry object
    geom.define(domain,&real_box,CoordSys::cartesian,is_periodic.data());
    geomC.define(domainC,&real_box,CoordSys::cartesian,is_periodic.data());
    geomP.define(domainP,&real_box,CoordSys::cartesian,is_periodic.data());


    // how boxes are distrubuted among MPI processes
    DistributionMapping dmap(ba);

    const Real* dx = geom.CellSize();
    const Real* dxc = geomC.CellSize();
    const Real* dxp = geomP.CellSize();

    MultiFab cellVols(bc, dmap, 1, 0);

#if (AMREX_SPACEDIM == 2)
#elif (AMREX_SPACEDIM == 3)
    cellVols.setVal(dxc[0]*dxc[1]*cell_depth);
    cellVols.setVal(dxc[0]*dxc[1]*dxc[2]);
#endif

    const RealBox& realDomain = geom.ProbDomain();

    Real dt = fixed_dt;
    Real dtinv = 1.0/dt;

    const int* lims = domainC.hiVect();

#if (AMREX_SPACEDIM == 2)
    int totalCollisionCells = (lims[0]+1)*(lims[1]+1);
#elif (AMREX_SPACEDIM == 3)
    int totalCollisionCells = (lims[0]+1)*(lims[1]+1)*(lims[2]+1);
#endif

#if (AMREX_SPACEDIM == 2)
    double domainVol = (prob_hi[0] - prob_lo[0])*(prob_hi[1] - prob_lo[1])*cell_depth;
#endif
#if (AMREX_SPACEDIM == 3)
    double domainVol = (prob_hi[0] - prob_lo[0])*(prob_hi[1] - prob_lo[1])*(prob_hi[2] - prob_lo[2]);
#endif

    species ionParticle[nspecies];

    double realParticles = 0;
    double simParticles = 0;

    for(int i=0;i<nspecies;i++)
    {       
        ionParticle[i].m = mass[i];
        ionParticle[i].d = diameter[i];

        ionParticle[i].q = qval;

        if(particle_count[i] >= 0)
        {
            ionParticle[i].ppb = (int)ceil((double)particle_count[i]/(double)ba.size());

            ionParticle[i].total = ionParticle[i].ppb*ba.size();

            ionParticle[i].n0 = ionParticle[i].total/domainVol;

            Print() << "Species " << i << " count adjusted to " << ionParticle[i].total << "\n";

        }else
        {
            ionParticle[i].n0 = particle_n0[i];
            ionParticle[i].total = (int)ceil(ionParticle[i].n0*domainVol/particle_neff);

            ionParticle[i].ppb = (int)ceil((double)ionParticle[i].total/(double)ba.size());

            ionParticle[i].total = ionParticle[i].ppb*ba.size();

            ionParticle[i].n0 = ionParticle[i].total/domainVol;

            Print() << "Species " << i << " n0 adjusted to " << ionParticle[i].n0 << "\n";
        }

        Print() << "Species " << i << " particles per box: " <<  ionParticle[i].ppb << "\n";

        realParticles = realParticles + ionParticle[i].total;
        simParticles = simParticles + ionParticle[i].total*particle_neff;

        ionParticle[i].Neff = particle_neff;

        ionParticle[i].propulsion = 0;

        ionParticle[i].gamma1 = 1.27;
        ionParticle[i].R = k_B/ionParticle[i].m;
        ionParticle[i].T = 273;
    }


    //Members
    //Density
    //velx
    //vely
    //velz
    //Temperature
    //jx
    //jy
    //jz
    //energyDensity
    //pressure


    MultiFab particleInstant(bc, dmap, 11, 0);

    //Members
    //Density
    //velx
    //vely
    //velz
    //Temperature
    //jx
    //jy
    //jz
    //energyDensity
    //pressure
    //speed
    
    MultiFab particleMeans(bc, dmap, 12, 0);

    //Members
    //Density
    //velx
    //vely
    //velz
    //Temperature
    //jx
    //jy
    //jz
    //energyDensity
    //pressure
    //GVar
    //KGCross
    //KRhoCross
    //RhoGCross
    //SpatialCross1
    //SpatialCross2
    //SpatialCross3
    //SpatialCross4
    //SpatialCross5
    //SpatialCross6

    
    MultiFab particleVars(bc, dmap, 21, 0);
    

    Print() << "Total real particles: " << realParticles << "\n";
    Print() << "Total sim particles: " << simParticles << "\n";

    Print() << "Sim particles per box: " << simParticles/(double)ba.size() << "\n";

    Print() << "Collision cells: " << totalCollisionCells << "\n";
    Print() << "Sim particles per cell: " << simParticles/totalCollisionCells << "\n";

    
    //-----------------------------
    //  Hydro setup

    ///////////////////////////////////////////
    // rho, alpha, beta, gamma:
    ///////////////////////////////////////////


    int ng = 1;

    if(pkernel_fluid == 4)
    {
        ng = 3;
    }else if(pkernel_fluid == 6)
    {

        ng = 4;
    }

    MultiFab rho(ba, dmap, 1, 1);
    rho.setVal(1.);

    // alpha_fc arrays
    std::array< MultiFab, AMREX_SPACEDIM > alpha_fc;
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        alpha_fc[d].define(convert(ba,nodal_flag_dir[d]), dmap, 1, ng);
        alpha_fc[d].setVal(dtinv);
    }

    // beta cell centred
    MultiFab beta(ba, dmap, 1, 1);
    beta.setVal(visc_coef);

    // beta on nodes in 2d
    // beta on edges in 3d
    std::array< MultiFab, NUM_EDGE > beta_ed;
#if (AMREX_SPACEDIM == 2)
    beta_ed[0].define(convert(ba,nodal_flag), dmap, 1, ng);
#elif (AMREX_SPACEDIM == 3)
    beta_ed[0].define(convert(ba,nodal_flag_xy), dmap, 1, ng);
    beta_ed[1].define(convert(ba,nodal_flag_xz), dmap, 1, ng);
    beta_ed[2].define(convert(ba,nodal_flag_yz), dmap, 1, ng);
#endif
    for (int d=0; d<NUM_EDGE; ++d) {
        beta_ed[d].setVal(visc_coef);
    }

    // cell-centered gamma
    MultiFab gamma(ba, dmap, 1, ng);
    gamma.setVal(0.);

    ///////////////////////////////////////////

    ///////////////////////////////////////////
    // Define & initalize eta & temperature multifabs
    ///////////////////////////////////////////
    // eta & temperature
    const Real eta_const = visc_coef;
    const Real temp_const = T_init[0];      // [units: K]

    // eta & temperature cell centered
    MultiFab  eta_cc;
    MultiFab temp_cc;
    // eta & temperature nodal
    std::array< MultiFab, NUM_EDGE >  eta_ed;
    std::array< MultiFab, NUM_EDGE > temp_ed;
    // eta and temperature; cell-centered
    eta_cc.define(ba, dmap, 1, 1);
    temp_cc.define(ba, dmap, 1, 1);
    // eta and temperature; nodal
#if (AMREX_SPACEDIM == 2)
    eta_ed[0].define(convert(ba,nodal_flag), dmap, 1, 0);
    temp_ed[0].define(convert(ba,nodal_flag), dmap, 1, 0);
#elif (AMREX_SPACEDIM == 3)
    eta_ed[0].define(convert(ba,nodal_flag_xy), dmap, 1, 0);
    eta_ed[1].define(convert(ba,nodal_flag_xz), dmap, 1, 0);
    eta_ed[2].define(convert(ba,nodal_flag_yz), dmap, 1, 0);
    temp_ed[0].define(convert(ba,nodal_flag_xy), dmap, 1, 0);
    temp_ed[1].define(convert(ba,nodal_flag_xz), dmap, 1, 0);
    temp_ed[2].define(convert(ba,nodal_flag_yz), dmap, 1, 0);
#endif

    // Initalize eta & temperature multifabs
    eta_cc.setVal(eta_const);
    temp_cc.setVal(temp_const);
    for (int d=0; d<NUM_EDGE; ++d) {
        eta_ed[d].setVal(eta_const);
        temp_ed[d].setVal(temp_const);
    }
    ///////////////////////////////////////////

    ///////////////////////////////////////////
    // random fluxes:
    ///////////////////////////////////////////


    // mflux divergence, staggered in x,y,z

    std::array< MultiFab, AMREX_SPACEDIM >  stochMfluxdiv;
    // Define mfluxdiv predictor multifabs
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        stochMfluxdiv[d].define(convert(ba,nodal_flag_dir[d]), dmap, 1, ng);
        stochMfluxdiv[d].setVal(0.0);
    }

    Vector< amrex::Real > weights;
    weights = {1.0};

    // Declare object of StochMFlux class
    StochMFlux sMflux (ba,dmap,geom,n_rngs);

    ///////////////////////////////////////////

    // pressure for GMRES solve
    MultiFab pres(ba,dmap,1,1);
    pres.setVal(0.);  // initial guess

    // staggered velocities
    std::array< MultiFab, AMREX_SPACEDIM > umac;
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        umac[d].define(convert(ba,nodal_flag_dir[d]), dmap, 1, ng);
    }


    ///////////////////////////////////////////
    // structure factor:
    ///////////////////////////////////////////

    
    Vector< std::string > var_names;
    var_names.resize(AMREX_SPACEDIM);
    int cnt = 0;
    std::string x;
    for (int d=0; d<var_names.size(); d++) {
      x = "vel";
      x += (120+d);
      var_names[cnt++] = x;
    }

    MultiFab struct_in_cc;
    struct_in_cc.define(ba, dmap, AMREX_SPACEDIM, 0);
    
    amrex::Vector< int > s_pairA(AMREX_SPACEDIM);
    amrex::Vector< int > s_pairB(AMREX_SPACEDIM);

    // Select which variable pairs to include in structure factor:
    s_pairA[0] = 0;
    s_pairB[0] = 0;
    //
    s_pairA[1] = 1;
    s_pairB[1] = 1;
    //
#if (AMREX_SPACEDIM == 3)
    s_pairA[2] = 2;
    s_pairB[2] = 2;
#endif
    
    StructFact structFact(ba,dmap,var_names);
    // StructFact structFact(ba,dmap,var_names,s_pairA,s_pairB);


    int dm = 0;
    for ( MFIter mfi(beta); mfi.isValid(); ++mfi ) {
        const Box& bx = mfi.validbox();

        AMREX_D_TERM(dm=0; init_vel(BL_TO_FORTRAN_BOX(bx),
                                    BL_TO_FORTRAN_ANYD(umac[0][mfi]), geom.CellSize(),
                                    geom.ProbLo(), geom.ProbHi() ,&dm,
                                    ZFILL(realDomain.lo()), ZFILL(realDomain.hi()));,
                     dm=1; init_vel(BL_TO_FORTRAN_BOX(bx),
                                    BL_TO_FORTRAN_ANYD(umac[1][mfi]), geom.CellSize(),
                                    geom.ProbLo(), geom.ProbHi() ,&dm,
                                    ZFILL(realDomain.lo()), ZFILL(realDomain.hi()));,
                     dm=2; init_vel(BL_TO_FORTRAN_BOX(bx),
                                    BL_TO_FORTRAN_ANYD(umac[2][mfi]), geom.CellSize(),
                                    geom.ProbLo(), geom.ProbHi() ,&dm,
                                    ZFILL(realDomain.lo()), ZFILL(realDomain.hi())););

    }

    AMREX_D_TERM(umac[0].setVal(0);,
                 umac[1].setVal(0);,
                 umac[2].setVal(0););

    // fill periodic ghost cells
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        umac[d].FillBoundary(geom.periodicity());
    }

    // Add initial equilibrium fluctuations
    if(initial_variance_mom != 0.0) {
        sMflux.addMfluctuations(umac, rho, temp_cc, initial_variance_mom, geom);
    }


    // staggered real coordinates - fluid grid
    std::array< MultiFab, AMREX_SPACEDIM > RealFaceCoords;
    AMREX_D_TERM(RealFaceCoords[0].define(convert(ba,nodal_flag_x), dmap, AMREX_SPACEDIM, ng);,
                 RealFaceCoords[1].define(convert(ba,nodal_flag_y), dmap, AMREX_SPACEDIM, ng);,
                 RealFaceCoords[2].define(convert(ba,nodal_flag_z), dmap, AMREX_SPACEDIM, ng););

    // staggered source terms - fluid grid
    std::array< MultiFab, AMREX_SPACEDIM > source;
    AMREX_D_TERM(source[0].define(convert(ba,nodal_flag_x), dmap, 1, ng);,
                 source[1].define(convert(ba,nodal_flag_y), dmap, 1, ng);,
                 source[2].define(convert(ba,nodal_flag_z), dmap, 1, ng););

    // staggered temporary holder for calculating source terms - This may not be necesssary, review later.
    std::array< MultiFab, AMREX_SPACEDIM > sourceTemp;
    AMREX_D_TERM(sourceTemp[0].define(convert(ba,nodal_flag_x), dmap, 1, ng);,
                 sourceTemp[1].define(convert(ba,nodal_flag_y), dmap, 1, ng);,
                 sourceTemp[2].define(convert(ba,nodal_flag_z), dmap, 1, ng););

    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        source[d].setVal(0.0);
    }

    int step = 0;
    Real time = 0.;
    int statsCount = 1;

    //Define parametric surfaces for particle interaction - declare array for surfaces and then define properties in BuildSurfaces


#if (BL_SPACEDIM == 3)
    int surfaceCount = 7;
    surface surfaceList[surfaceCount];
    BuildSurfaces(surfaceList,surfaceCount,realDomain.lo(),realDomain.hi());
#endif
#if (BL_SPACEDIM == 2)
    int surfaceCount = 5;
    surface surfaceList[surfaceCount];
    BuildSurfaces(surfaceList,surfaceCount,realDomain.lo(),realDomain.hi());
#endif


    //Particles! Build on geom & box array for collision cells/ poisson grid?
    FhdParticleContainer particles(geomC, dmap, bc);

    //Find coordinates of cell faces. May be used for interpolating fields to particle locations
    FindFaceCoords(RealFaceCoords, geom); //May not be necessary to pass Geometry?

    //create particles
    particles.InitParticles(ionParticle[0]);

    //particles.InitializeFields(particleInstant, cellVols, ionParticle[0]);

    //setup initial DSMC collision parameters
    //particles.InitCollisionCells(collisionPairs, collisionFactor, cellVols, ionParticle[0], dt);


    //----------------------    
    // Electrostatic setup
    //----------------------

    int ngp = 1;

    if(pkernel_es == 4)
    {
        ngp = 3;

    }else if(pkernel_es == 6)
    {

        ngp = 4;
    }

    // cell centered real coordinates - es grid
    MultiFab RealCenteredCoords;
    RealCenteredCoords.define(bp, dmap, AMREX_SPACEDIM, ngp);

    FindCenterCoords(RealCenteredCoords, geomP);
    
    //Cell centred es potential
    MultiFab potential(ba, dmap, 1, ngp);
    MultiFab potentialTemp(ba, dmap, 1, ngp);
    potential.setVal(0);
    potentialTemp.setVal(0);

    MultiFab charge(ba, dmap, 1, ngp);
    MultiFab chargeTemp(ba, dmap, 1, ngp);
    charge.setVal(0);
    chargeTemp.setVal(0);

    //Staggered electric field - probably wont use this?
    std::array< MultiFab, AMREX_SPACEDIM > efield;
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        efield[d].define(convert(bp,nodal_flag_dir[d]), dmap, 1, ngp);
    }

    MLPoisson EsSolver;

    // write out initial state
    //WritePlotFile(step,time,geom,geomC,rhotot,umac,div,particleMembers,particleDensity,particleVelocity, particleTemperature, particlePressure, particleSpatialCross1, particleMembraneFlux, particles);

    particles.MoveIons(dt, dx, geom.ProbLo(), umac, RealFaceCoords, source, sourceTemp, surfaceList, surfaceCount, 2 /*1: interpolate only. 2: spread only. 3: both*/ );
    particles.Redistribute();
    particles.ReBin();

    //Time stepping loop

//    dt = dt*10e4;

    for(step=1;step<=max_step;++step)
    {


//        if(step==8)
//        {
//            dt = dt*10e-4;
//        }

        //HYDRO
        //--------------------------------------

	    // Fill stochastic terms
	    if(variance_coef_mom != 0.0) {

	      // compute the random numbers needed for the stochastic momentum forcing
	      sMflux.fillMStochastic();


	      // compute stochastic momentum force
	      sMflux.stochMforce(stochMfluxdiv,eta_cc,eta_ed,temp_cc,temp_ed,weights,dt);
	    }

	    // Advance umac
        //Print() << "STOKES SOLVE\n";
	    advance(umac,pres,stochMfluxdiv,source,alpha_fc,beta,gamma,beta_ed,geom,dt);

        if (plot_int > 0 && step%plot_int == 0)
        {
           
            //This write particle data and associated fields
            WritePlotFile(step,time,geom,geomC, particleInstant, particleMeans, particleVars, particles);

            //Writes instantaneous flow field and some other stuff? Check with Guy.
            WritePlotFileHydro(step,time,geom,umac,pres);
        }


        particles.MoveIons(dt, dx, geom.ProbLo(), umac, RealFaceCoords, source, sourceTemp, surfaceList, surfaceCount, 1 /*1: interpolate only. 2: spread only. 3: both*/ );

        particles.Redistribute();

        particles.ReBin();

       //Particles
        //--------------------------------------


        //Probably don't need to pass ProbLo(), check later.


        //particles.CollideParticles(collisionPairs, collisionFactor, cellVols, ionParticle[0], dt);

//        if(step == n_steps_skip)
//        {
//            particleMeans.setVal(0.0);
//            particleVars.setVal(0);
//            statsCount = 1;
//        }
       
        //particles.EvaluateStats(particleInstant, particleMeans, particleVars, delHolder1, delHolder2, delHolder3, delHolder4, delHolder5, delHolder6, particleMembraneFlux, cellVols, ionParticle[0], dt,statsCount);

        statsCount++;

        if(step%1 == 0)
        {    
                amrex::Print() << "Advanced step " << step << "\n";
        }
        
        time = time + dt;



    }

    // Call the timer again and compute the maximum difference between the start time 
    // and stop time over all processors
    Real stop_time = ParallelDescriptor::second() - strt_time;
    ParallelDescriptor::ReduceRealMax(stop_time);
    amrex::Print() << "Run time = " << stop_time << std::endl;

    amrex::Finalize();
}
