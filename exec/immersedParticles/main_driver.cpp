

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

using namespace gmres;
using namespace common;

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

    //RealBox real_box({AMREX_D_DECL(0,0,0)},
    //                 {AMREX_D_DECL( 1.0, 1.0, 1.0)});


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

    //DistributionMapping dmapc(bc);

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
    activeParticle.n0 = 16/domainVol;
    activeParticle.propulsion = 0;

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
    

    // how boxes are distrubuted among MPI processes
    DistributionMapping dmap(ba);

    // total density
    MultiFab rhotot(ba, dmap, 1, 1);

    // divergence
    MultiFab div(ba, dmap, 1, 1);

    // potential
    MultiFab phi(ba, dmap, 1, 1);

    // beta cell centred
    MultiFab betaCC(ba, dmap, 1, 1);

    // Nodal velocity for interpolations
    std::array< MultiFab, AMREX_SPACEDIM > umacNodal;
    AMREX_D_TERM(umacNodal[0].define(convert(ba,IntVect{AMREX_D_DECL(1, 1, 1)}), dmap, 1, 1);,
                 umacNodal[1].define(convert(ba,IntVect{AMREX_D_DECL(1, 1, 1)}), dmap, 1, 1);,
                 umacNodal[2].define(convert(ba,IntVect{AMREX_D_DECL(1, 1, 1)}), dmap, 1, 1););

    // gamma cell centred
    MultiFab gammaCC(ba, dmap, 1, 1);

    //Print() << nodal_flag_xy << "\n";
    //while(true);

    // beta on nodes in 2d
    // beta on edges in 3d
    std::array< MultiFab, NUM_EDGE > betaEdge;

    double tempVisc = 9e-4;

#if (AMREX_SPACEDIM == 2)
    betaEdge[0].define(convert(ba,nodal_flag), dmap, 1, 1);
    betaEdge[0].setVal(tempVisc);
#elif (AMREX_SPACEDIM == 3)
    betaEdge[0].define(convert(ba,nodal_flag_xy), dmap, 1, 1);
    betaEdge[1].define(convert(ba,nodal_flag_xz), dmap, 1, 1);
    betaEdge[2].define(convert(ba,nodal_flag_yz), dmap, 1, 1);
    betaEdge[0].setVal(tempVisc);  
    betaEdge[1].setVal(tempVisc);
    betaEdge[2].setVal(tempVisc);
#endif

    //Nodal beta. If running in 2D, betaEdge is already nodal.

#if (AMREX_SPACEDIM == 3)
    MultiFab betaNodal(convert(ba,IntVect{AMREX_D_DECL(1, 1, 1)}), dmap, 1, 1);
#endif

    //Replace with proper initialiser
    phi.setVal(100.);

    betaCC.setVal(1.);
    gammaCC.setVal(0);

    // temporary placeholder for potential gradients on cell faces
    std::array< MultiFab, AMREX_SPACEDIM > umacT;
    AMREX_D_TERM(umacT[0].define(convert(ba,nodal_flag_x), dmap, 1, 1);,
                 umacT[1].define(convert(ba,nodal_flag_y), dmap, 1, 1);,
                 umacT[2].define(convert(ba,nodal_flag_z), dmap, 1, 1););

    // set density to 1
    rhotot.setVal(1.);

    // staggered real coordinates
    std::array< MultiFab, AMREX_SPACEDIM > RealFaceCoords;
    AMREX_D_TERM(RealFaceCoords[0].define(convert(ba,nodal_flag_x), dmap, AMREX_SPACEDIM, 1);,
                 RealFaceCoords[1].define(convert(ba,nodal_flag_y), dmap, AMREX_SPACEDIM, 1);,
                 RealFaceCoords[2].define(convert(ba,nodal_flag_z), dmap, AMREX_SPACEDIM, 1););

    // staggered velocities
    std::array< MultiFab, AMREX_SPACEDIM > umac;
    AMREX_D_TERM(umac[0].define(convert(ba,nodal_flag_x), dmap, 1, 1);,
                 umac[1].define(convert(ba,nodal_flag_y), dmap, 1, 1);,
                 umac[2].define(convert(ba,nodal_flag_z), dmap, 1, 1););

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

    // alpha arrays
    std::array< MultiFab, AMREX_SPACEDIM > alpha;
    AMREX_D_TERM(alpha[0].define(convert(ba,nodal_flag_x), dmap, 1, 1);
                 alpha[0].setVal(0);,
                 alpha[1].define(convert(ba,nodal_flag_y), dmap, 1, 1);
                 alpha[1].setVal(0);,
                 alpha[2].define(convert(ba,nodal_flag_z), dmap, 1, 1);
                 alpha[2].setVal(0););

    // For testing timestepping
    std::array< MultiFab, AMREX_SPACEDIM > umacNew;
    AMREX_D_TERM(umacNew[0].define(convert(ba,nodal_flag_x), dmap, 1, 1);,
                 umacNew[1].define(convert(ba,nodal_flag_y), dmap, 1, 1);,
                 umacNew[2].define(convert(ba,nodal_flag_z), dmap, 1, 1););

    // For testing StagApplyOp
    std::array< MultiFab, AMREX_SPACEDIM > umacOut;
    AMREX_D_TERM(umacOut[0].define(convert(ba,nodal_flag_x), dmap, 1, 1);,
                 umacOut[1].define(convert(ba,nodal_flag_y), dmap, 1, 1);,
                 umacOut[2].define(convert(ba,nodal_flag_z), dmap, 1, 1););


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


	int dm = 0;
	for ( MFIter mfi(rhotot); mfi.isValid(); ++mfi )
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

    //Find coordinates of cell faces. Used for interpolating fields to particle locations
    FindFaceCoords(RealFaceCoords, geom); //May not be necessary to pass Geometry?

    particles.InitParticles(ppc, activeParticle);

    //particles.InitializeFields(particleInstant, cellVols, nitrogen);

    particles.InitCollisionCells(collisionPairs, collisionFactor, cellVols, activeParticle, dt);

    // write out initial state
    //WritePlotFile(step,time,geom,geomC,rhotot,umac,div,particleMembers,particleDensity,particleVelocity, particleTemperature, particlePressure, particleSpatialCross1, particleMembraneFlux, particles);

    //Time stepping loop
    for(step=1;step<=max_step;++step)
    {
        AMREX_D_TERM(
        umac[0].FillBoundary(geom.periodicity());,
        umac[1].FillBoundary(geom.periodicity());,
        umac[2].FillBoundary(geom.periodicity()););

        //eulerStep(betaCC, gammaCC, 
         //        betaEdge,
         //        umac, umacOut, umacNew, alpha, geom, &dt);

        AMREX_D_TERM(
        MultiFab::Copy(umac[0], umacNew[0], 0, 0, 1, 0);,
        MultiFab::Copy(umac[1], umacNew[1], 0, 0, 1, 0);,
        MultiFab::Copy(umac[2], umacNew[2], 0, 0, 1, 0););
        
        particles.MoveParticles(dt, surfaceList, surfaceCount);

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


    //Compute divergence, last arguement is flag for increment
    
    //ComputeDiv(div,umac,geom,0);

    //Compute gradient
    //ComputeGrad(phi,umacT,geom);


    // Call the timer again and compute the maximum difference between the start time 
    // and stop time over all processors
    Real stop_time = ParallelDescriptor::second() - strt_time;
    ParallelDescriptor::ReduceRealMax(stop_time);
    amrex::Print() << "Run time = " << stop_time << std::endl;

    amrex::Finalize();
}
