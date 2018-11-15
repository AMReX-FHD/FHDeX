

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

#include "hydro_test_functions.H"
#include "hydro_test_functions_F.H"

#include "hydro_functions.H"
#include "hydro_functions_F.H"

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

    //Initialise rngs
    rng_initialize(&fhdSeed,&particleSeed,&selectorSeed,&thetaSeed,&phiSeed,&generalSeed);

    // is the problem periodic?
    Vector<int> is_periodic(AMREX_SPACEDIM,0);  // set to 0 (not periodic) by default
   /* for (int i=0; i<AMREX_SPACEDIM; ++i) {
        if (bc_lo[i] == -1 && bc_hi[i] == -1) {
            is_periodic[i] = 1;
        }
    }*/

    //---------------------------------------

        //Terms prepended with a 'C' are related to the particle grid. Those without are for the hydro grid.
        //The particle grid created as a corsening or refinement of the hydro grid.

    //---------------------------------------

    // make BoxArray and Geometry
    BoxArray ba;
    BoxArray bc;
    Geometry geom;
    Geometry geomC;
    
    IntVect dom_lo(AMREX_D_DECL(           0,            0,            0));
    IntVect dom_hi(AMREX_D_DECL(n_cells[0]-1, n_cells[1]-1, n_cells[2]-1));
    Box domain(dom_lo, dom_hi);
    Box domainC = domain;

    // Initialize the boxarray "ba" from the single box "bx"
    ba.define(domain);

    // Break up boxarray "ba" into chunks no larger than "max_grid_size" along a direction
    // note we are converting "Vector<int> max_grid_size" to an IntVect
    ba.maxSize(IntVect(max_grid_size));

    RealBox real_box({AMREX_D_DECL(prob_lo[0],prob_lo[1],prob_lo[2])},
                     {AMREX_D_DECL(prob_hi[0],prob_hi[1],prob_hi[2])});

    //This must be an even number for now?
    int sizeRatio = 1;

    bc = ba;
    bc.coarsen(sizeRatio);
    domainC.coarsen(sizeRatio);
    // This defines a Geometry object
    geom.define(domain,&real_box,CoordSys::cartesian,is_periodic.data());
    geomC.define(domainC,&real_box,CoordSys::cartesian,is_periodic.data());

    // how boxes are distrubuted among MPI processes
    DistributionMapping dmap(ba);

    const Real* dx = geom.CellSize();
    const Real* dxc = geomC.CellSize();

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

    species activeParticle;

#if (AMREX_SPACEDIM == 2)
    double domainVol = (prob_hi[0] - prob_lo[0])*(prob_hi[1] - prob_lo[1])*cell_depth;
#endif
#if (AMREX_SPACEDIM == 3)
    double domainVol = (prob_hi[0] - prob_lo[0])*(prob_hi[1] - prob_lo[1])*(prob_hi[2] - prob_lo[2]);
#endif

    //1 micron particles with density of gold
       
    activeParticle.m = 1.5158e-8;
    activeParticle.d = 1e-6;        
    activeParticle.n0 = 200/domainVol;
    activeParticle.propulsion = 0;

    activeParticle.gamma1 = 1.27;
    activeParticle.R = 1.3806E-23/activeParticle.m;
    activeParticle.T = 273;
    //nitrogen.mu = 2.1E-5;
    //nitrogen.d = sqrt((nitrogen.gamma1*nitrogen.m*sqrt(nitrogen.R*nitrogen.T))/(4*sqrt(3.14159265359)*nitrogen.mu));
    //nitrogen.d = 3.66e-10;        

    //nitrogen.mfp = (6.26e-8);

    //nitrogen.n0 = 1.0/(sqrt(2)*3.14159265359*nitrogen.d*nitrogen.d*nitrogen.mfp);

    //nitrogen.n0 = 1.78/nitrogen.m;

    //Print() << nitrogen.n0*nitrogen.m << "\n";

    //nitrogen.P = nitrogen.R*nitrogen.m*nitrogen.n0*nitrogen.T;

    //nitrogen.cp = sqrt(2.0*nitrogen.R*nitrogen.T);

    double realParticles = domainVol*activeParticle.n0;

    Print() << "Real particles: " << realParticles << "\n";

    activeParticle.Neff = 1;
    const int totalParticles = realParticles;

    const int ppc  = (int)ceil(((double)totalParticles)/((double)totalCollisionCells));
    const int ppb = (int)ceil(((double)totalParticles)/((double)ba.size()));

    Print() << "Total particles: " << totalParticles << "\n";
    Print() << "Particles per box: " << ppb << "\n";

    Print() << "Collision cells: " << totalCollisionCells << "\n";
    Print() << "Particles per cell: " << ppc << "\n";

    Print() << "rho0: " << activeParticle.m*activeParticle.n0 << "\n";
    Print() << "Adjusted rho0: " << activeParticle.m*activeParticle.Neff*totalParticles/domainVol << "\n";


    MultiFab collisionPairs(bc, dmap, 1, 0);    
    MultiFab collisionFactor(bc, dmap, 1, 0);
    MultiFab particleMembraneFlux(bc, dmap, 1, 0);

    collisionFactor.setVal(0);
    collisionPairs.setVal(0);

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
    
    MultiFab particleVars(bc, dmap, 18, 0);

    particleInstant.setVal(0.0);
    particleMeans.setVal(0.0);
    particleVars.setVal(0.0);
    
    //-----------------------------
    //  Hydro setup

    ///////////////////////////////////////////
    // rho, alpha, beta, gamma:
    ///////////////////////////////////////////
    
    MultiFab rho(ba, dmap, 1, 1);
    rho.setVal(1.);

    MultiFab rhoNodal(convert(ba,nodal_flag), dmap, 1, 1);
    rhoNodal.setVal(1.);

    // alpha_fc arrays
    Real theta_alpha = 1.;
    std::array< MultiFab, AMREX_SPACEDIM > alpha_fc;
    AMREX_D_TERM(alpha_fc[0].define(convert(ba,nodal_flag_x), dmap, 1, 1);,
                 alpha_fc[1].define(convert(ba,nodal_flag_y), dmap, 1, 1);,
                 alpha_fc[2].define(convert(ba,nodal_flag_z), dmap, 1, 1););
    AMREX_D_TERM(alpha_fc[0].setVal(dtinv);,
                 alpha_fc[1].setVal(dtinv);,
                 alpha_fc[2].setVal(dtinv););

    // beta cell centred
    MultiFab beta(ba, dmap, 1, 1);
    beta.setVal(visc_coef);

    // beta on nodes in 2d
    // beta on edges in 3d
    std::array< MultiFab, NUM_EDGE > beta_ed;
#if (AMREX_SPACEDIM == 2)
    beta_ed[0].define(convert(ba,nodal_flag), dmap, 1, 1);
    beta_ed[0].setVal(visc_coef);
#elif (AMREX_SPACEDIM == 3)
    beta_ed[0].define(convert(ba,nodal_flag_xy), dmap, 1, 1);
    beta_ed[1].define(convert(ba,nodal_flag_xz), dmap, 1, 1);
    beta_ed[2].define(convert(ba,nodal_flag_yz), dmap, 1, 1);
    beta_ed[0].setVal(visc_coef);
    beta_ed[1].setVal(visc_coef);
    beta_ed[2].setVal(visc_coef);
#endif

    //Nodal beta for interpolations
    MultiFab betaNodal(convert(ba,nodal_flag), dmap, 1, 1);
    betaNodal.setVal(visc_coef);

    // cell-centered gamma
    MultiFab gamma(ba, dmap, 1, 1);
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
    std::array< MultiFab, NUM_EDGE >   eta_ed;
    std::array< MultiFab, NUM_EDGE >  temp_ed;
    // eta cell-centered
    eta_cc.define(ba, dmap, 1, 1);
    // temperature cell-centered
    temp_cc.define(ba, dmap, 1, 1);
#if (AMREX_SPACEDIM == 2)
    // eta nodal
    eta_ed[0].define(convert(ba,nodal_flag), dmap, 1, 0);
    // temperature nodal
    temp_ed[0].define(convert(ba,nodal_flag), dmap, 1, 0);
#elif (AMREX_SPACEDIM == 3)
    // eta nodal
    eta_ed[0].define(convert(ba,nodal_flag_xy), dmap, 1, 0);
    eta_ed[1].define(convert(ba,nodal_flag_xz), dmap, 1, 0);
    eta_ed[2].define(convert(ba,nodal_flag_yz), dmap, 1, 0);
    // temperature nodal
    temp_ed[0].define(convert(ba,nodal_flag_xy), dmap, 1, 0);
    temp_ed[1].define(convert(ba,nodal_flag_xz), dmap, 1, 0);
    temp_ed[2].define(convert(ba,nodal_flag_yz), dmap, 1, 0);
#endif

    // Initalize eta & temperature multifabs
    // eta cell-centered
    eta_cc.setVal(eta_const);
    // temperature cell-centered
    temp_cc.setVal(temp_const);
#if (AMREX_SPACEDIM == 2)
    // eta nodal
    eta_ed[0].setVal(eta_const);
    // temperature nodal
    temp_ed[0].setVal(temp_const);
#elif (AMREX_SPACEDIM == 3)
    // eta nodal
    eta_ed[0].setVal(eta_const);
    eta_ed[1].setVal(eta_const);
    eta_ed[2].setVal(eta_const);
    // temperature nodal
    temp_ed[0].setVal(temp_const);
    temp_ed[1].setVal(temp_const);
    temp_ed[2].setVal(temp_const);
#endif
    ///////////////////////////////////////////

    ///////////////////////////////////////////
    // random fluxes:
    ///////////////////////////////////////////

    // mflux divergence, staggered in x,y,z

    std::array< MultiFab, AMREX_SPACEDIM >  mfluxdiv_predict;
    // Define mfluxdiv predictor multifabs
    mfluxdiv_predict[0].define(convert(ba,nodal_flag_x), dmap, 1, 1);
    mfluxdiv_predict[1].define(convert(ba,nodal_flag_y), dmap, 1, 1);
#if (AMREX_SPACEDIM == 3)
    mfluxdiv_predict[2].define(convert(ba,nodal_flag_z), dmap, 1, 1);
#endif

    std::array< MultiFab, AMREX_SPACEDIM >  mfluxdiv_correct;
    // Define mfluxdiv corrector multifabs
    mfluxdiv_correct[0].define(convert(ba,nodal_flag_x), dmap, 1, 1);
    mfluxdiv_correct[1].define(convert(ba,nodal_flag_y), dmap, 1, 1);
#if (AMREX_SPACEDIM == 3)
    mfluxdiv_correct[2].define(convert(ba,nodal_flag_z), dmap, 1, 1);
#endif

    Vector< amrex::Real > weights;
    // weights = {std::sqrt(0.5), std::sqrt(0.5)};
    weights = {1.0};
    
    // Declare object of StochMFlux class 
    StochMFlux sMflux (ba,dmap,geom,n_rngs);

    ///////////////////////////////////////////

    // tracer
    MultiFab tracer(ba,dmap,1,1);

    // pressure for GMRES solve
    MultiFab pres(ba,dmap,1,1);
    pres.setVal(0.);  // initial guess

    // staggered velocities
    std::array< MultiFab, AMREX_SPACEDIM > umac;
    AMREX_D_TERM(umac[0].define(convert(ba,nodal_flag_x), dmap, 1, 1);,
                 umac[1].define(convert(ba,nodal_flag_y), dmap, 1, 1);,
                 umac[2].define(convert(ba,nodal_flag_z), dmap, 1, 1););

    std::array< MultiFab, AMREX_SPACEDIM > umacNew;
    AMREX_D_TERM(umacNew[0].define(convert(ba,nodal_flag_x), dmap, 1, 1);,
                 umacNew[1].define(convert(ba,nodal_flag_y), dmap, 1, 1);,
                 umacNew[2].define(convert(ba,nodal_flag_z), dmap, 1, 1););


    //Nodal velocity for interpolations
    std::array< MultiFab, AMREX_SPACEDIM > umacNodal;
    AMREX_D_TERM(umacNodal[0].define(convert(ba,IntVect{AMREX_D_DECL(1, 1, 1)}), dmap, 1, 1);,
                 umacNodal[1].define(convert(ba,IntVect{AMREX_D_DECL(1, 1, 1)}), dmap, 1, 1);,
                 umacNodal[2].define(convert(ba,IntVect{AMREX_D_DECL(1, 1, 1)}), dmap, 1, 1););

    // temporary placeholder for potential gradients on cell faces
    std::array< MultiFab, AMREX_SPACEDIM > umacT;
    AMREX_D_TERM(umacT[0].define(convert(ba,nodal_flag_x), dmap, 1, 1);,
                 umacT[1].define(convert(ba,nodal_flag_y), dmap, 1, 1);,
                 umacT[2].define(convert(ba,nodal_flag_z), dmap, 1, 1););

    // staggered real coordinates
    std::array< MultiFab, AMREX_SPACEDIM > RealFaceCoords;
    AMREX_D_TERM(RealFaceCoords[0].define(convert(ba,nodal_flag_x), dmap, AMREX_SPACEDIM, 1);,
                 RealFaceCoords[1].define(convert(ba,nodal_flag_y), dmap, AMREX_SPACEDIM, 1);,
                 RealFaceCoords[2].define(convert(ba,nodal_flag_z), dmap, AMREX_SPACEDIM, 1););

    // staggered source terms
    std::array< MultiFab, AMREX_SPACEDIM > source;
    AMREX_D_TERM(source[0].define(convert(ba,nodal_flag_x), dmap, 1, 1);,
                 source[1].define(convert(ba,nodal_flag_y), dmap, 1, 1);,
                 source[2].define(convert(ba,nodal_flag_z), dmap, 1, 1););

    // staggered temporary holder for calculating source terms - This may not be necesssary, review later.
    std::array< MultiFab, AMREX_SPACEDIM > sourceTemp;
    AMREX_D_TERM(sourceTemp[0].define(convert(ba,nodal_flag_x), dmap, 1, 1);,
                 sourceTemp[1].define(convert(ba,nodal_flag_y), dmap, 1, 1);,
                 sourceTemp[2].define(convert(ba,nodal_flag_z), dmap, 1, 1););

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
	for ( MFIter mfi(beta); mfi.isValid(); ++mfi )
    {
        const Box& bx = mfi.validbox();

        AMREX_D_TERM(dm=0; init_vel(BL_TO_FORTRAN_BOX(bx),
                                    BL_TO_FORTRAN_ANYD(umac[0][mfi]), geom.CellSize(),
                                    geom.ProbLo(), geom.ProbHi() ,&dm, ZFILL(realDomain.lo()), ZFILL(realDomain.hi()));,
                     dm=1; init_vel(BL_TO_FORTRAN_BOX(bx),
            			            BL_TO_FORTRAN_ANYD(umac[1][mfi]), geom.CellSize(),
            			            geom.ProbLo(), geom.ProbHi() ,&dm, ZFILL(realDomain.lo()), ZFILL(realDomain.hi()));,
                     dm=2; init_vel(BL_TO_FORTRAN_BOX(bx),
                                    BL_TO_FORTRAN_ANYD(umac[2][mfi]), geom.CellSize(),
                                    geom.ProbLo(), geom.ProbHi() ,&dm, ZFILL(realDomain.lo()), ZFILL(realDomain.hi())););
    }

    AMREX_D_TERM(umac[0].setVal(0);,
                 umac[1].setVal(0);,
                 umac[2].setVal(0););

    AMREX_D_TERM(
    MultiFab::Copy(umacNew[0], umac[0], 0, 0, 1, 0);,
    MultiFab::Copy(umacNew[1], umac[1], 0, 0, 1, 0);,
    MultiFab::Copy(umacNew[2], umac[2], 0, 0, 1, 0););

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


    //Particles! Build on geom & box array for collision cells
    FhdParticleContainer particles(geomC, dmap, bc);

    //Find coordinates of cell faces. May be used for interpolating fields to particle locations
    FindFaceCoords(RealFaceCoords, geom); //May not be necessary to pass Geometry?

    //create particles
    particles.InitParticles(ppc, activeParticle);

    particles.InitializeFields(particleInstant, cellVols, activeParticle);

    //setup initial DSMC collision parameters
    particles.InitCollisionCells(collisionPairs, collisionFactor, cellVols, activeParticle, dt);


    // write out initial state
    //WritePlotFile(step,time,geom,geomC,rhotot,umac,div,particleMembers,particleDensity,particleVelocity, particleTemperature, particlePressure, particleSpatialCross1, particleMembraneFlux, particles);

    //Time stepping loop
    for(step=1;step<=max_step;++step)
    {

        //HYDRO
        //--------------------------------------

	    // Fill stochastic terms
//	    sMflux.fillMStochastic();

	    // compute stochastic force terms
///	    sMflux.stochMforce(mfluxdiv_predict,eta_cc,eta_ed,temp_cc,temp_ed,weights,dt);
//	    sMflux.stochMforce(mfluxdiv_correct,eta_cc,eta_ed,temp_cc,temp_ed,weights,dt);
	
	    // Advance umac
 //           advance(umac,umacNew,pres,tracer,mfluxdiv_predict,mfluxdiv_correct,
//		    alpha_fc,beta,gamma,beta_ed,geom,dt);
	
	    //////////////////////////////////////////////////
	
	    ///////////////////////////////////////////
	    // Update structure factor
	    ///////////////////////////////////////////
//	    if (step > n_steps_skip && struct_fact_int > 0 && (step-n_steps_skip-1)%struct_fact_int == 0) {
//	      for(int d=0; d<AMREX_SPACEDIM; d++) {
//	        ShiftFaceToCC(umac[d], 0, struct_in_cc, d, 1);
//	      }
//	      structFact.FortStructure(struct_in_cc,geom);
//     }
	    ///////////////////////////////////////////


       //Particles
        //--------------------------------------


        //Probably don't need to pass ProbLo(), check later.
        
        particles.MoveParticles(dt, dx, geom.ProbLo(), umac, umacNodal, RealFaceCoords, beta, betaNodal, rho, rhoNodal, source, sourceTemp, surfaceList, surfaceCount);

        particles.Redistribute();

        particles.ReBin();

        particles.CollideParticles(collisionPairs, collisionFactor, cellVols, activeParticle, dt);

        if(step == n_steps_skip)
        {
            particleMeans.setVal(0.0);
            particleVars.setVal(0);
            statsCount = 1;
        }
       
        particles.EvaluateStats(particleInstant, particleMeans, particleVars, particleMembraneFlux, cellVols, activeParticle, dt,statsCount);

        statsCount++;

        if(step%5000 == 0)
        {    
                amrex::Print() << "Advanced step " << step << "\n";
        }
        
        time = time + dt;

        if (plot_int > 0 && step%plot_int == 0)
        {
           
          //  WritePlotFile(step,time,geom,geomConvert, particleInstant, particleMeans, particleVars, particleMembraneFlux, particles);
        }

    }

    // Call the timer again and compute the maximum difference between the start time 
    // and stop time over all processors
    Real stop_time = ParallelDescriptor::second() - strt_time;
    ParallelDescriptor::ReduceRealMax(stop_time);
    amrex::Print() << "Run time = " << stop_time << std::endl;

    amrex::Finalize();
}
