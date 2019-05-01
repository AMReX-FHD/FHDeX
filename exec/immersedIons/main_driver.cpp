#include "INS_functions.H"

#include "common_functions.H"
#include "gmres_functions.H"

#include "common_functions_F.H"
#include "common_namespace.H"
#include "common_namespace_declarations.H"

#include "gmres_functions_F.H"
#include "gmres_namespace.H"
#include "gmres_namespace_declarations.H"


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

//#include "electrostatic.H"

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

//    int fhdSeed = ParallelDescriptor::MyProc() + 1;
//    int particleSeed = 2*ParallelDescriptor::MyProc() + 2;
//    int selectorSeed = 3*ParallelDescriptor::MyProc() + 3;
//    int thetaSeed = 4*ParallelDescriptor::MyProc() + 4;
//    int phiSeed = 5*ParallelDescriptor::MyProc() + 5;
//    int generalSeed = 6*ParallelDescriptor::MyProc() + 6;

    int fhdSeed = 0;
    int particleSeed = 0;
    int selectorSeed = 0;
    int thetaSeed = 0;
    int phiSeed = 0;
    int generalSeed = 0;

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

    // A - fluid, C - particle, P/E - electrostatic

    // AJN what is the particle grid used for?  Is it leftover from DSCM collision model?
    
    // AJN - change these to ba_fluid, ba_particle, etc.    
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
    bc = ba;
    bp = ba;

    int sizeRatio;

    if(particle_grid_refine < 1)
    {
        sizeRatio = (int)(1.0/particle_grid_refine);
        bc.refine(sizeRatio);
        domainC.refine(sizeRatio);
    }else
    {
        sizeRatio = (int)(particle_grid_refine);
        bc.coarsen(sizeRatio);
        domainC.coarsen(sizeRatio);
    }

    if(es_grid_refine < 1)
    {
        sizeRatio = (int)(1.0/es_grid_refine);
        bp.refine(sizeRatio);
        domainP.refine(sizeRatio);
    }else
    {
        sizeRatio = (int)(es_grid_refine);
        bp.coarsen(sizeRatio);
        domainP.coarsen(sizeRatio);
    }

    // This defines a Geometry object
    geom.define(domain,&real_box,CoordSys::cartesian,is_periodic.data());
    geomC.define(domainC,&real_box,CoordSys::cartesian,is_periodic.data());
    geomP.define(domainP,&real_box,CoordSys::cartesian,is_periodic.data());


    // how boxes are distrubuted among MPI processes
    // AJN needs to be fi
    DistributionMapping dmap(ba);

    Print() << geom << "\n";
    Print() << domain << "\n";

    const Real* dx = geom.CellSize();
    const Real* dxc = geomC.CellSize();
    const Real* dxp = geomP.CellSize();

    // AJN - does cellVols have to be a MultiFab (could it just a Real?)
    MultiFab cellVols(bc, dmap, 1, 0);

#if (AMREX_SPACEDIM == 2)
    cellVols.setVal(dxc[0]*dxc[1]*cell_depth);
#elif (AMREX_SPACEDIM == 3)
    cellVols.setVal(dxc[0]*dxc[1]*dxc[2]);
#endif

    // AJN - isn't this the same as real_box?
    const RealBox& realDomain = geom.ProbDomain();

    Real dt = fixed_dt;
    Real dtinv = 1.0/dt;

    const int* lims = domainC.hiVect();

    // AJN - get rid of collision stuff?
#if (AMREX_SPACEDIM == 2)
    int totalCollisionCells = (lims[0]+1)*(lims[1]+1);
#elif (AMREX_SPACEDIM == 3)
    int totalCollisionCells = (lims[0]+1)*(lims[1]+1)*(lims[2]+1);
#endif

#if (AMREX_SPACEDIM == 2)
    double domainVol = (prob_hi[0] - prob_lo[0])*(prob_hi[1] - prob_lo[1])*cell_depth;
#elif (AMREX_SPACEDIM == 3)
    double domainVol = (prob_hi[0] - prob_lo[0])*(prob_hi[1] - prob_lo[1])*(prob_hi[2] - prob_lo[2]);
#endif

    //Species type defined in species.H
    //array of length nspecies

    species ionParticle[nspecies];

    double realParticles = 0;
    double simParticles = 0;
    double dryRad, wetRad;
    double dxAv = (dx[0] + dx[1] + dx[2])/3.0; //This is probably the wrong way to do this.

    if(pkernel_fluid == 3)
    {
        wetRad = 0.9405*dxAv;                
    }
    else if(pkernel_fluid == 4)
    {
        wetRad = 1.3273*dxAv;
    }
    else if(pkernel_fluid == 6)
    {
        wetRad = 1.5705*dxAv;
    }

    for(int i=0;i<nspecies;i++) {
        
        ionParticle[i].m = mass[i];
        ionParticle[i].q = qval[i];

        if(diameter[i] > 0)
        {
            ionParticle[i].d = diameter[i];

            ionParticle[i].wetDiff = (k_B*T_init[0])/(6*3.14159265359*wetRad*visc_coef);

            ionParticle[i].totalDiff = (k_B*T_init[0])/(6*3.14159265359*(diameter[i]/2.0)*visc_coef);

            ionParticle[i].dryDiff = ionParticle[i].totalDiff - ionParticle[i].wetDiff; //This is probably wrong. Need to test.

        }else
        {
            ionParticle[i].totalDiff = diff[i];            

            ionParticle[i].d = 2.0*(k_B*T_init[0])/(6*3.14159265359*(ionParticle[i].totalDiff)*visc_coef);

            ionParticle[i].wetDiff = (k_B*T_init[0])/(6*3.14159265359*wetRad*visc_coef);

            ionParticle[i].dryDiff = ionParticle[i].totalDiff - ionParticle[i].wetDiff; //Test this
        }

        Print() << "Species " << i << " wet diffusion: " << ionParticle[i].wetDiff << ", dry diffusion: " << ionParticle[i].dryDiff << ", total: " << ionParticle[i].totalDiff << ", hydro radius: " << ionParticle[i].d/2.0 << "\n";

        ionParticle[i].Neff = particle_neff; // From DSMC, this will be set to 1 for electolyte calcs
        ionParticle[i].R = k_B/ionParticle[i].m; //used a lot in kinetic stats cals, bu not otherwise necessary for electrolytes

        ionParticle[i].sigma = sigma[i];
        ionParticle[i].eepsilon = eepsilon[i];

        // AJN - why round up particles so there are the same number in each box?
        if(particle_count[i] >= 0) {
            // adjust number of particles up so there is the same number per box            
            ionParticle[i].ppb = (int)ceil((double)particle_count[i]/(double)ba.size());
            ionParticle[i].total = ionParticle[i].ppb*ba.size();
            ionParticle[i].n0 = ionParticle[i].total/domainVol;

            Print() << "Species " << i << " count adjusted to " << ionParticle[i].total << "\n";
        }
        else {
            // if particle count is negative, we instead compute the number of particles based on particle density and particle_neff
            ionParticle[i].total = (int)ceil(particle_n0[i]*domainVol/particle_neff);
            // adjust number of particles up so there is the same number per box  
            ionParticle[i].ppb = (int)ceil((double)ionParticle[i].total/(double)ba.size());
            ionParticle[i].total = ionParticle[i].ppb*ba.size();
            ionParticle[i].n0 = ionParticle[i].total/domainVol;

            Print() << "Species " << i << " n0 adjusted to " << ionParticle[i].n0 << "\n";
        }

        Print() << "Species " << i << " particles per box: " <<  ionParticle[i].ppb << "\n";

        realParticles = realParticles + ionParticle[i].total;
        simParticles = simParticles + ionParticle[i].total*particle_neff;
    }
    
    Print() << "Total real particles: " << realParticles << "\n";
    Print() << "Total sim particles: " << simParticles << "\n";

    Print() << "Sim particles per box: " << simParticles/(double)ba.size() << "\n";

    Print() << "Collision cells: " << totalCollisionCells << "\n";
    Print() << "Sim particles per cell: " << simParticles/totalCollisionCells << "\n";


    // MFs for storing particle statistics

    // A lot of these relate to gas kinetics, but many are still useful so leave in for now.
    
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
    
    //-----------------------------
    //  Hydro setup

    ///////////////////////////////////////////
    // rho, alpha, beta, gamma:
    ///////////////////////////////////////////

    //set number of ghost cells to fit whole peskin kernel
    int ang = 1;
    int cng = 2;
    int png = 1;

    // AJN - for perdictor/corrector do we need one more ghost cell if the predictor pushes a particle into a ghost region?
    
    if(pkernel_fluid == 3) {
        ang = 2;
    }
    else if(pkernel_fluid == 4) {
        ang = 3;
    }
    else if(pkernel_fluid == 6) {
        ang = 4;
    }

    // AJN this only needs 1? ghost cell    
    MultiFab rho(ba, dmap, 1, 1);
    rho.setVal(1.);

    // AJN - the number of ghost cells needed for the GMRES solve for alpha, beta, etc., is a fixed number
    // and not dependent on the Peskin kernels.
    // alpha_fc -> 1 ghost cell
    // beta -> 1
    // beta_ed -> 1
    // gamma -> 1
    // I'm thinking only the velocities need the extra ghost cells.  Probably not density
    
    // alpha_fc arrays
    std::array< MultiFab, AMREX_SPACEDIM > alpha_fc;
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        alpha_fc[d].define(convert(ba,nodal_flag_dir[d]), dmap, 1, 1);
        alpha_fc[d].setVal(dtinv);
    }

    // beta cell centred
    MultiFab beta(ba, dmap, 1, 1);
    beta.setVal(visc_coef);

    // beta on nodes in 2d
    // beta on edges in 3d
    std::array< MultiFab, NUM_EDGE > beta_ed;
#if (AMREX_SPACEDIM == 2)
    beta_ed[0].define(convert(ba,nodal_flag), dmap, 1, 1);
#elif (AMREX_SPACEDIM == 3)
    beta_ed[0].define(convert(ba,nodal_flag_xy), dmap, 1, 1);
    beta_ed[1].define(convert(ba,nodal_flag_xz), dmap, 1, 1);
    beta_ed[2].define(convert(ba,nodal_flag_yz), dmap, 1, 1);
#endif
    for (int d=0; d<NUM_EDGE; ++d) {
        beta_ed[d].setVal(visc_coef);
    }

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

    // AJN - I think the number of ghost cells needed here is
    // eta_cc -> 1
    // temp_cc -> 1
    // eta_ed -> 0
    // temp_ed -> 0
    
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
    eta_ed[0].define(convert(ba,nodal_flag), dmap, 1, 1);
    temp_ed[0].define(convert(ba,nodal_flag), dmap, 1, 1);
#elif (AMREX_SPACEDIM == 3)
    eta_ed[0].define(convert(ba,nodal_flag_xy), dmap, 1, 1);
    eta_ed[1].define(convert(ba,nodal_flag_xz), dmap, 1, 1);
    eta_ed[2].define(convert(ba,nodal_flag_yz), dmap, 1, 1);
    temp_ed[0].define(convert(ba,nodal_flag_xy), dmap, 1, 1);
    temp_ed[1].define(convert(ba,nodal_flag_xz), dmap, 1, 1);
    temp_ed[2].define(convert(ba,nodal_flag_yz), dmap, 1, 1);
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

    // AJN - I think stochMfluxdiv needs 0 ghost cells //DRL setting it to 1 until I've checked.

    // mflux divergence, staggered in x,y,z

    std::array< MultiFab, AMREX_SPACEDIM >  stochMfluxdiv;
    // Define mfluxdiv predictor multifabs
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        stochMfluxdiv[d].define(convert(ba,nodal_flag_dir[d]), dmap, 1, 1);
        stochMfluxdiv[d].setVal(0.0);
    }

    Vector< amrex::Real > weights;
    weights = {1.0};

    // Declare object of StochMFlux class
    StochMFlux sMflux (ba,dmap,geom,n_rngs);

    ///////////////////////////////////////////

    // AJN - pres needs 1 ghost cell
    // but umac is the thing that needs extra ghost cells for Peskin kernels
    
    // pressure for GMRES solve
    MultiFab pres(ba,dmap,1,1);
    pres.setVal(0.);  // initial guess

    // staggered velocities
    std::array< MultiFab, AMREX_SPACEDIM > umac;
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        umac[d].define(convert(ba,nodal_flag_dir[d]), dmap, 1, ang);
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

    
    // AJN - don't need to initialize velocities in overdamped.  first gmres solve should get them as long as they start out with non-NaN values.
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


    // AJN - don't need this
    // Add initial equilibrium fluctuations
    if(initial_variance_mom != 0.0) {
        sMflux.addMfluctuations(umac, rho, temp_cc, initial_variance_mom, geom);
    }


    // staggered real coordinates - fluid grid
    std::array< MultiFab, AMREX_SPACEDIM > RealFaceCoords;
    AMREX_D_TERM(RealFaceCoords[0].define(convert(ba,nodal_flag_x), dmap, AMREX_SPACEDIM, ang);,
                 RealFaceCoords[1].define(convert(ba,nodal_flag_y), dmap, AMREX_SPACEDIM, ang);,
                 RealFaceCoords[2].define(convert(ba,nodal_flag_z), dmap, AMREX_SPACEDIM, ang););

    // staggered source terms - fluid grid
    std::array< MultiFab, AMREX_SPACEDIM > source;
    AMREX_D_TERM(source[0].define(convert(ba,nodal_flag_x), dmap, 1, ang);,
                 source[1].define(convert(ba,nodal_flag_y), dmap, 1, ang);,
                 source[2].define(convert(ba,nodal_flag_z), dmap, 1, ang););

    // staggered temporary holder for calculating source terms - This may not be necesssary, review later.
    std::array< MultiFab, AMREX_SPACEDIM > sourceTemp;
    AMREX_D_TERM(sourceTemp[0].define(convert(ba,nodal_flag_x), dmap, 1, ang);,
                 sourceTemp[1].define(convert(ba,nodal_flag_y), dmap, 1, ang);,
                 sourceTemp[2].define(convert(ba,nodal_flag_z), dmap, 1, ang););

    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        source[d].setVal(0.0);
    }

    int step = 0;
    Real time = 0.;
    int statsCount = 1;

    //Define parametric surfaces for particle interaction - declare array for surfaces and then define properties in BuildSurfaces


    // AJN - we don't understand why you need this for ions
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

    //int num_neighbor_cells = 4; replaced by input var
    //Particles! Build on geom & box array for collision cells/ poisson grid?
    FhdParticleContainer particles(geomC, dmap, bc, crange);

    //Find coordinates of cell faces (fluid grid). May be used for interpolating fields to particle locations
    FindFaceCoords(RealFaceCoords, geom); //May not be necessary to pass Geometry?

    //create particles
    particles.InitParticles(ionParticle);

    //----------------------    
    // Electrostatic setup
    //----------------------

    // AJN - should define 3 types of "ng" parameters.  fluid, repulsive force, peskin
    int ngp = 1;
    if(pkernel_es == 3)
    {
        ngp = 2;

    }else if(pkernel_es == 4)
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

    // AJN - what are all the Temp's for?
    
    //Cell centred es potential
    MultiFab potential(bp, dmap, 1, ngp);
    MultiFab potentialTemp(bp, dmap, 1, ngp);
    potential.setVal(0);
    potentialTemp.setVal(0);

    //charage density for RHS of Poisson Eq.
    MultiFab charge(bp, dmap, 1, ngp);
    MultiFab chargeTemp(bp, dmap, 1, ngp);
    charge.setVal(0);
    chargeTemp.setVal(0);

    //mass density on ES grid - not necessary
    MultiFab massFrac(bp, dmap, 1, 1);
    MultiFab massFracTemp(bp, dmap, 1, 1);
    massFrac.setVal(0);
    massFracTemp.setVal(0);

//    //?????
//    MultiFab trans(bp, dmap, 1, ngp);
//    MultiFab transTemp(bp, dmap, 1, ngp);
//    trans.setVal(0);
//    transTemp.setVal(0);

    //Staggered electric fields
    std::array< MultiFab, AMREX_SPACEDIM > efield;
    std::array< MultiFab, AMREX_SPACEDIM > external;

    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        efield[d].define(convert(bp,nodal_flag_dir[d]), dmap, 1, ngp);
       
    }

    //Centred electric fields
    std::array< MultiFab, AMREX_SPACEDIM > efieldCC;

    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        efieldCC[d].define(bp, dmap, 1, ngp);
        external[d].define(bp, dmap, 1, ngp);
    }

    AMREX_D_TERM(efieldCC[0].setVal(0);,
                 efieldCC[1].setVal(0);,
                 efieldCC[2].setVal(0););

    //Apply external field here.
    AMREX_D_TERM(external[0].setVal(0);,
                 external[1].setVal(0);,
                 external[2].setVal(0););

 
    //Time stepping loop
    for(step=1;step<=max_step;++step)
    {

        //Most of these functions are sensitive to the order of execution. We can fix this, but for now leave them in this order.
        particles.clearNeighbors();

        particles.DoRFD(dt, dx, dxp, geom, umac, efieldCC, RealFaceCoords, RealCenteredCoords, source, sourceTemp, surfaceList, surfaceCount, 3 /*this number currently does nothing, but we will use it later*/);
        std::cout << ParallelDescriptor::MyProc() << ", here2\n";
        particles.fillNeighbors();
        std::cout << ParallelDescriptor::MyProc() << ", here3\n";
        particles.computeForcesNL();
        std::cout << ParallelDescriptor::MyProc() << ", here4\n";


        if(es_tog==1)
        {
            //Spreads charge density from ions onto multifab 'charge'.
            particles.collectFields(dt, dxp, RealCenteredCoords, geomP, charge, chargeTemp, massFrac, massFracTemp);

            //Do Poisson solve using 'charge' for RHS, and put potential in 'potential'. Then calculate gradient and put in 'efield', then add 'external'.
            esSolve(potential, charge, efieldCC, external, geomP);
        }

        //compute other forces and spread to grid
        particles.SpreadIons(dt, dx, dxp, geom, umac, efieldCC, RealFaceCoords, RealCenteredCoords, source, sourceTemp, surfaceList, surfaceCount, 3 /*this number currently does nothing, but we will use it later*/);

        if((variance_coef_mom != 0.0) && fluid_tog == 1) {
          // compute the random numbers needed for the stochastic momentum forcing
          sMflux.fillMStochastic();
//          // compute stochastic momentum force
          sMflux.stochMforce(stochMfluxdiv,eta_cc,eta_ed,temp_cc,temp_ed,weights,dt);
        }
        if(fluid_tog ==1)
        {
    	    advance(umac,pres,stochMfluxdiv,source,alpha_fc,beta,gamma,beta_ed,geom,dt);
        }

        if(move_tog==1)
        {
            //Calls wet ion interpolation and movement.
            particles.MoveIons(dt, dx, dxp, geom, umac, efield, RealFaceCoords, source, sourceTemp, surfaceList, surfaceCount, 3 /*this number currently does nothing, but we will use it later*/);

            particles.Redistribute();
          //  particles.ReBin();            //We may not need to redist & rebin after seperately for wet & dry moves - check this later

            //particles.MoveParticlesDry(dt, dx, umac, RealFaceCoords, source, sourceTemp, surfaceList, surfaceCount);

            //These functions reorganise particles between cells and processes
    //Dout  particles.Redistribute();
    //Dout  particles.ReBin();
        }

        //Start collecting statistics after step n_steps_skip
        if(step == n_steps_skip)
        {
            particleMeans.setVal(0.0);
            particleVars.setVal(0);
            statsCount = 1;
        }
       
        particles.EvaluateStats(particleInstant, particleMeans, particleVars, cellVols, ionParticle[0], dt,statsCount);

        statsCount++;

        if (plot_int > 0 && step%plot_int == 0)
        {
           
            //This write particle data and associated fields
            WritePlotFile(step,time,geom,geomC,geomP,particleInstant, particleMeans, particleVars, particles, charge, potential, efieldCC);

            //Writes instantaneous flow field and some other stuff? Check with Guy.
            WritePlotFileHydro(step,time,geom,umac,pres);
        }

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
}
