

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

using namespace common;
using namespace gmres;
using namespace amrex;

// argv contains the name of the inputs file entered at the command line
void main_driver(const char* argv)
{

    // store the current time so we can later compute total run time.
  //  Real strt_time = ParallelDescriptor::second();

    int argc = 0;

    char* ptr1 = (char*)argv;

    char** ptr2 = &ptr1;
        
    amrex::Initialize(argc,ptr2);

    std::string inputs_file = argv;

    // read in parameters from inputs file into F90 modules
    // we use "+1" because of amrex_string_c_to_f expects a null char termination
    read_common_namelist(inputs_file.c_str(),inputs_file.size()+1);
    read_gmres_namelist(inputs_file.c_str(),inputs_file.size()+1);

    // copy contents of F90 modules to C++ namespaces
    InitializeCommonNamespace();
    InitializeGmresNamespace();

//    int fhdSeed = ParallelDescriptor::MyProc()*20 + 1 + 10000;
//    int particleSeed = 2*ParallelDescriptor::MyProc()*30 + 2 + 10000;
//    int selectorSeed = 3*ParallelDescriptor::MyProc()*40 + 3 + 10000;
//    int thetaSeed = 4*ParallelDescriptor::MyProc()*50 + 4 + 10000;
//    int phiSeed = 5*ParallelDescriptor::MyProc()*60 + 5 + 10000;
//    int generalSeed = 6*ParallelDescriptor::MyProc()*70 + 6 + 10000;

    int fhdSeed = 0;
    int particleSeed = 0;
    int selectorSeed = 0;
    int thetaSeed = 0;
    int phiSeed = 0;
    int generalSeed = 0;



    //Initialise rngs
    rng_initialize(&fhdSeed,&particleSeed,&selectorSeed,&thetaSeed,&phiSeed,&generalSeed);

    // is the problem periodic?
    Vector<int> is_periodic(AMREX_SPACEDIM,1);  // set to 0 (not periodic) by default
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
    Geometry geomConvert;
    
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

    int cgsConv = 100;

    RealBox real_box_convert({AMREX_D_DECL(prob_lo[0]*cgsConv,prob_lo[1]*cgsConv,prob_lo[2]*cgsConv)},
                     {AMREX_D_DECL(prob_hi[0]*cgsConv,prob_hi[1]*cgsConv,prob_hi[2]*cgsConv)});


    //This must be an even number for now?
    int sizeRatio = 1;

    bc = ba;
    bc.coarsen(sizeRatio);
    domainC.coarsen(sizeRatio);
    // This defines a Geometry object
    geom.define(domain,&real_box,CoordSys::cartesian,is_periodic.data());
    geomC.define(domainC,&real_box,CoordSys::cartesian,is_periodic.data());

    geomConvert.define(domainC,&real_box_convert,CoordSys::cartesian,is_periodic.data());

    //DistributionMapping dmapc(bc);

    const int* lims = domainC.hiVect();

#if (AMREX_SPACEDIM == 2)
    int totalCollisionCells = (lims[0]+1)*(lims[1]+1);
#elif (AMREX_SPACEDIM == 3)
    int totalCollisionCells = (lims[0]+1)*(lims[1]+1)*(lims[2]+1);
#endif

    species nitrogen;

    //Hard sphere nitrogen gas
    /*nitrogen.gamma1 = 1.27;
    nitrogen.m = 4.65E-26;
    nitrogen.R = 1.3806E-23/nitrogen.m;
    nitrogen.T = 293;
    nitrogen.mu = 1.757E-5;
    nitrogen.d = sqrt((nitrogen.gamma1*nitrogen.m*sqrt(nitrogen.R*nitrogen.T))/(4*sqrt(3.14159265359)*nitrogen.mu));    

    nitrogen.mfp = 5e-7;

    nitrogen.n0 = 1.0/(sqrt(2)*3.14159265359*nitrogen.d*nitrogen.d*nitrogen.mfp);

    nitrogen.P = nitrogen.R*nitrogen.m*nitrogen.n0*nitrogen.T;

    nitrogen.cp = sqrt(2.0*nitrogen.R*nitrogen.T);*/

    nitrogen.gamma1 = 1.27;
    nitrogen.m = 6.63E-26;
    nitrogen.R = 1.3806E-23/nitrogen.m;
    nitrogen.T = 273;
    nitrogen.mu = 2.1E-5;
    //nitrogen.d = sqrt((nitrogen.gamma1*nitrogen.m*sqrt(nitrogen.R*nitrogen.T))/(4*sqrt(3.14159265359)*nitrogen.mu));
    nitrogen.d = 3.66e-10;        

    nitrogen.mfp = (6.26e-8);

    //nitrogen.n0 = 1.0/(sqrt(2)*3.14159265359*nitrogen.d*nitrogen.d*nitrogen.mfp);

    nitrogen.n0 = 1.78/nitrogen.m;

    //Print() << nitrogen.n0*nitrogen.m << "\n";

    nitrogen.P = nitrogen.R*nitrogen.m*nitrogen.n0*nitrogen.T;

    nitrogen.cp = sqrt(2.0*nitrogen.R*nitrogen.T);
    

#if (AMREX_SPACEDIM == 2)
    double domainVol = (prob_hi[0] - prob_lo[0])*(prob_hi[1] - prob_lo[1])*cell_depth;
#endif
#if (AMREX_SPACEDIM == 3)
    double domainVol = (prob_hi[0] - prob_lo[0])*(prob_hi[1] - prob_lo[1])*(prob_hi[2] - prob_lo[2]);
#endif

    double realParticles = domainVol*nitrogen.n0;

    Print() << "Real particles: " << realParticles << "\n";

    nitrogen.Neff = 1;
    const int totalParticles = realParticles;

    Print() << "Neff: " << nitrogen.Neff << "\n";

    const int ppc  = (int)ceil(((double)totalParticles)/((double)totalCollisionCells));
    const int ppb = (int)ceil(((double)totalParticles)/((double)ba.size()));

    nitrogen.ppb = ppb;    

    Print() << "Total particles: " << totalParticles << "\n";
    Print() << "Particles per box: " << ppb << "\n";

    Print() << "Collision cells: " << totalCollisionCells << "\n";
    Print() << "Particles per cell: " << ppc << "\n";

    Print() << "rho0: " << nitrogen.m*nitrogen.n0 << "\n";
    Print() << "Adjusted rho0: " << nitrogen.m*nitrogen.Neff*totalParticles/domainVol << "\n";
    

    // how boxes are distrubuted among MPI processes
    DistributionMapping dmap(ba);

    MultiFab collisionPairs(bc, dmap, 1, 0);    
    MultiFab collisionFactor(bc, dmap, 1, 0);

    collisionFactor.setVal(0);


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

    Real delHolder1[n_cells[1]*n_cells[2]];
    Real delHolder2[n_cells[1]*n_cells[2]];
    Real delHolder3[n_cells[1]*n_cells[2]];
    Real delHolder4[n_cells[1]*n_cells[2]];
    Real delHolder5[n_cells[1]*n_cells[2]];
    Real delHolder6[n_cells[1]*n_cells[2]];


    MultiFab particleMembraneFlux(bc, dmap, 1, 0);


    particleInstant.setVal(0.0);
    particleMeans.setVal(0.0);
    particleVars.setVal(0.0);
    particleMembraneFlux.setVal(0.0);

    const Real* dx = geom.CellSize();
    const Real* dxc = geomC.CellSize();

    MultiFab cellVols(bc, dmap, 1, 0);

#if (AMREX_SPACEDIM == 2)
    cellVols.setVal(dxc[0]*dxc[1]*cell_depth);
#elif (AMREX_SPACEDIM == 3)
    cellVols.setVal(dxc[0]*dxc[1]*dxc[2]);
#endif

    const RealBox& realDomain = geom.ProbDomain();


    // kinetic time step
    //Real dt = 0.1*nitrogen.mfp/sqrt(2.0*nitrogen.R*nitrogen.T);
    Real dt = fixed_dt;

    //dt = nitrogen.mfp/sqrt(2.0*nitrogen.R*nitrogen.T);
    Print() << "Step size: " << dt << "\n";

    //while(true);

    int step = 0;
    Real time = 0.;
    int statsCount = 1;

    //Define parametric surfaces for particle interaction - declare array for surfaces and then define properties in BuildSurfaces

#if (BL_SPACEDIM == 3)
    int surfaceCount = 6;
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

    particles.InitParticles(nitrogen);

    //particles.UpdateCellVectors();

    particles.InitializeFields(particleInstant, cellVols, nitrogen);

    particles.InitCollisionCells(collisionPairs, collisionFactor, cellVols, nitrogen, dt);

    // write out initial state
    WritePlotFile(step,time,geom,geomConvert, particleInstant, particleMeans, particleVars, particleMembraneFlux, particles);

    //Time stepping loop
    for(step=1;step<=max_step;++step)
    {


        particles.MoveParticles(dt, surfaceList, surfaceCount);


   //     Print() << "Here2\n";
        particles.Redistribute();
        particles.ReBin();


        particles.CollideParticles(collisionPairs, collisionFactor, cellVols, nitrogen, dt);

        //if((step-10)%20000 == 0)
        if(step == n_steps_skip)
        {
            particleMeans.setVal(0.0);
            particleVars.setVal(0);
            statsCount = 1;
        }

       
        particles.EvaluateStats(particleInstant, particleMeans, particleVars, delHolder1, delHolder2, delHolder3, delHolder4, delHolder5, delHolder6, particleMembraneFlux, cellVols, nitrogen, dt,statsCount);

        statsCount++;


        if(step%500 == 0)
        {    
                amrex::Print() << "Advanced step " << step << "\n";
        }
        
        time = time + dt;

        if (plot_int > 0 && step%plot_int == 0)
        {
           
            WritePlotFile(step,time,geom,geomConvert, particleInstant, particleMeans, particleVars, particleMembraneFlux, particles);
        }


    }

    // Call the timer again and compute the maximum difference between the start time 
    // and stop time over all processors
    //Real stop_time = ParallelDescriptor::second() - strt_time;
   // ParallelDescriptor::ReduceRealMax(stop_time);
   // amrex::Print() << "Run time = " << stop_time << std::endl;

    amrex::Finalize();
}
