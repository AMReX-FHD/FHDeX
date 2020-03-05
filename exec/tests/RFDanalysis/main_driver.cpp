#include "INS_functions.H"
#include <iostream>

#include "common_functions.H"
#include "gmres_functions.H"

#include "common_namespace_declarations.H"

#include "gmres_namespace_declarations.H"



#include "species.H"
#include "surfaces.H"

//#include "analysis_functions_F.H"
//#include "StructFact_F.H"
#include "StochMomFlux.H"
//#include "StructFact.H"

#include "hydro_test_functions_F.H"

#include "hydro_functions.H"

#include "electrostatic.H"

#include "particle_functions.H"

//#include "electrostatic.H"

#include "debug_functions_F.H"
#include "AMReX_ArrayLim.H"

//#include <IBMarkerContainer.H>

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

    remove("potential.dat");
    remove("kinetic.dat");

    const int n_rngs = 1;

    int fhdSeed = 0;
    int particleSeed = 0;
    int selectorSeed = 0;
    int thetaSeed = 0;
    int phiSeed = 0;
    int generalSeed = 0;

    if(seed > 0)
    {
        fhdSeed = ParallelDescriptor::MyProc() + 1 + seed;
        particleSeed = 2*ParallelDescriptor::MyProc() + 2 + seed;
        selectorSeed = 3*ParallelDescriptor::MyProc() + 3 + seed;
        thetaSeed = 4*ParallelDescriptor::MyProc() + 4 + seed;
        phiSeed = 5*ParallelDescriptor::MyProc() + 5 + seed;
        generalSeed = 6*ParallelDescriptor::MyProc() + 6 + seed;
    }

    //Initialise rngs
    rng_initialize(&fhdSeed,&particleSeed,&selectorSeed,&thetaSeed,&phiSeed,&generalSeed);

    // is the problem periodic?
    Vector<int> is_periodic(AMREX_SPACEDIM,0);  // set to 0 (not periodic) by default
    Vector<int> is_periodic_c(AMREX_SPACEDIM,0);  // set to 0 (not periodic) by default
    Vector<int> is_periodic_p(AMREX_SPACEDIM,0);  // set to 0 (not periodic) by default
    for (int i=0; i<AMREX_SPACEDIM; ++i) {
        if (bc_vel_lo[i] == -1 && bc_vel_hi[i] == -1) {
            is_periodic[i] = 1;
            is_periodic_c[i] = 1;
        }
        if (bc_es_lo[i] == -1 && bc_es_hi[i] == -1) {
            is_periodic_p[i] = 1;
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
    geomC.define(domainC,&real_box,CoordSys::cartesian,is_periodic_c.data());
    geomP.define(domainP,&real_box,CoordSys::cartesian,is_periodic_p.data());

    // how boxes are distrubuted among MPI processes
    // AJN needs to be fi
    DistributionMapping dmap(ba);

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

        Print() << "Species " << i << " wet diffusion: " << ionParticle[i].wetDiff << ", dry diffusion: " << ionParticle[i].dryDiff << ", total diffusion: " << ionParticle[i].totalDiff << "\n";
        Print() << "Species " << i << " wet radius: " << wetRad << ", dry radius: " << (k_B*T_init[0])/(6*3.14159265359*(ionParticle[i].dryDiff)*visc_coef) << ", total radius: " << ionParticle[i].d/2.0 << "\n";

        if(ionParticle[i].dryDiff < 0)
        {
            Print() << "Negative dry diffusion in species " << i << "\n";
            abort();
        }

        ionParticle[i].Neff = particle_neff; // From DSMC, this will be set to 1 for electolyte calcs
        ionParticle[i].R = k_B/ionParticle[i].m; //used a lot in kinetic stats cals, bu not otherwise necessary for electrolytes

        ionParticle[i].sigma = sigma[i];
        ionParticle[i].eepsilon = eepsilon[i];

        // AJN - why round up particles so there are the same number in each box? DRL - Have to divide them into whole numbers of particles somehow. 
        if(particle_count[i] >= 0) {


            ionParticle[i].ppb = (double)particle_count[i]/(double)ba.size();
            ionParticle[i].total = particle_count[i];
            ionParticle[i].n0 = ionParticle[i].total/domainVol;
            
            Print() << "Species " << i << " count adjusted to " << ionParticle[i].total << "\n";
        }
        else {
            // if particle count is negative, we instead compute the number of particles based on particle density and particle_neff
            ionParticle[i].total = (int)ceil(particle_n0[i]*domainVol/particle_neff);
            // adjust number of particles up so there is the same number per box  
            ionParticle[i].ppb = (int)ceil((double)ionParticle[i].total/(double)ba.size());
            //ionParticle[i].total = ionParticle[i].ppb*ba.size();
            ionParticle[i].n0 = ionParticle[i].total/domainVol;

            Print() << "Species " << i << " n0 adjusted to " << ionParticle[i].n0 << "\n";
        }

        Print() << "Species " << i << " particles per box: " <<  ionParticle[i].ppb << "\n";

        realParticles = realParticles + ionParticle[i].total*particle_neff;
        simParticles = simParticles + ionParticle[i].total;
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
    //Cx
    //Cy
    //Cz
    MultiFab particleInstant(bc, dmap, 14, 0);

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
    //Cx
    //Cy
    //Cz    
    MultiFab particleMeans(bc, dmap, 14, 0);

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
    //Cx
    //Cy
    //Cz 
   
    MultiFab particleVars(bc, dmap, 18, 0);
    
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
    std::array< MultiFab, AMREX_SPACEDIM >  stochMfluxdivC;
    // Define mfluxdiv predictor/corrector multifabs
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        stochMfluxdiv[d].define(convert(ba,nodal_flag_dir[d]), dmap, 1, 1);
        stochMfluxdiv[d].setVal(0.0);
        stochMfluxdivC[d].define(convert(ba,nodal_flag_dir[d]), dmap, 1, 1);
        stochMfluxdivC[d].setVal(0.0);
    }

    Vector< amrex::Real > weights;
    weights = {1.0};
//weights = {std::sqrt(0.5), std::sqrt(0.5)};

    // Declare object of StochMomFlux class
    StochMomFlux sMflux (ba,dmap,geom,n_rngs);

    ///////////////////////////////////////////

    // AJN - pres needs 1 ghost cell
    // but umac is the thing that needs extra ghost cells for Peskin kernels
    
    // pressure for GMRES solve
    MultiFab pres(ba,dmap,1,1);
   pres.setVal(0.);  // initial guess

    // staggered velocities
    std::array< MultiFab, AMREX_SPACEDIM > umac;
    std::array< MultiFab, AMREX_SPACEDIM > umacNew;
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        umac[d].define(convert(ba,nodal_flag_dir[d]), dmap, 1, ang);
        umacNew[d].define(convert(ba,nodal_flag_dir[d]), dmap, 1, ang);
    }


    // staggered mfabs for storing some basic fluid stats
    std::array< MultiFab, AMREX_SPACEDIM > umacM;
    std::array< MultiFab, AMREX_SPACEDIM > umacV;
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        umacM[d].define(convert(ba,nodal_flag_dir[d]), dmap, 1, ang);
        umacV[d].define(convert(ba,nodal_flag_dir[d]), dmap, 1, ang);
    }

    // tracer - get rid of this.
    MultiFab tracer(ba, dmap, 1,1);
    tracer.setVal(0.);


    ///////////////////////////////////////////
    // structure factor:
    ///////////////////////////////////////////

    Vector< std::string > var_names;
    int nvar_sf = 1;
    // int nvar_sf = AMREX_SPACEDIM;
    var_names.resize(nvar_sf);
    var_names[0] = "charge";

    MultiFab struct_in_cc;
    struct_in_cc.define(bp, dmap, nvar_sf, 0);

//    amrex::Vector< int > s_pairA(nvar_sf);
//    amrex::Vector< int > s_pairB(nvar_sf);

//    // Select which variable pairs to include in structure factor:
//    for (int d=0; d<nvar_sf; d++) {
//      s_pairA[d] = d;
//      s_pairB[d] = d;
//    }

 //   StructFact structFact(bp,dmap,var_names);






    // AJN - don't need to initialize velocities in overdamped.  first gmres solve should get them as long as they start out with non-NaN values.

    // DRL - This is actually useful for dubugging, to get a known velocity field.
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

            AMREX_D_TERM(umacM[0].setVal(0);,
                     umacM[1].setVal(0);,
                     umacM[2].setVal(0););

            AMREX_D_TERM(umacV[0].setVal(0);,
                     umacV[1].setVal(0);,
                     umacV[2].setVal(0););

    // fill periodic ghost cells
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        umac[d].FillBoundary(geom.periodicity());
    }


    // AJN - don't need this
    // Add initial equilibrium fluctuations
    if(initial_variance_mom != 0.0) {
        //sMflux.addMfluctuations(umac, rho, temp_cc, initial_variance_mom);
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

	// IBMarkerContainerBase default behaviour is to do tiling. Turn off here:

    Vector<int> ts(BL_SPACEDIM);

    //AMREX_D_TERM(ts[0] = max_grid_size[0];, ts[1] = max_grid_size[1];, ts[2] = max_grid_size[2];);

    AMREX_D_TERM(
        if(max_particle_tile_size[0] > 0)
        {
            ts[0] = max_particle_tile_size[0];
        }else
        {
            ts[0] = max_grid_size[0];
        },
        if(max_particle_tile_size[1] > 0)
        {
            ts[1] = max_particle_tile_size[1];
        }else
        {
            ts[1] = max_grid_size[1];
        },
        if(max_particle_tile_size[2] > 0)
        {
            ts[2] = max_particle_tile_size[2];
        }else
        {
            ts[2] = max_grid_size[2];
        }
    );

	ParmParse pp ("particles");
    pp.addarr("tile_size", ts);



	//pp.add("do_tiling", false);

    //int num_neighbor_cells = 4; replaced by input var
    //Particles! Build on geom & box array for collision cells/ poisson grid?
    FhdParticleContainer particles(geomC, dmap, bc, crange);

    //Find coordinates of cell faces (fluid grid). May be used for interpolating fields to particle locations
    FindFaceCoords(RealFaceCoords, geom); //May not be necessary to pass Geometry?

    //create particles

    particles.InitParticles(ionParticle, dxp);

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
    
    MultiFab dryMobility(ba, dmap, nspecies*AMREX_SPACEDIM, ang);

    ComputeDryMobility(dryMobility, ionParticle, geom);
 
    //Time stepping loop
    for(step=1;step<=max_step;++step)
    {

        //Most of these functions are sensitive to the order of execution. We can fix this, but for now leave them in this order.

        //Apply external field here.
        AMREX_D_TERM(external[0].setVal(eamp[0]*cos(efreq[0]*time + ephase[0]));,
                     external[1].setVal(eamp[1]*cos(efreq[1]*time + ephase[1]));,
                     external[2].setVal(eamp[2]*cos(efreq[2]*time + ephase[2])););

        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            source[d].setVal(0.0);
            sourceTemp[d].setVal(0.0);
        }

        if(rfd_tog==1) {
            // Apply RFD force to fluid
            particles.RFD(0, dx, sourceTemp, RealFaceCoords);
            particles.ResetMarkers(0);
            //particles.DoRFD(dt, dx, dxp, geom, umac, efieldCC, RealFaceCoords, RealCenteredCoords, source, sourceTemp, surfaceList, surfaceCount, 3 /*this number currently does nothing, but we will use it later*/);
        }
        else {
            // set velx/y/z and forcex/y/z for each particle to zero
            particles.ResetMarkers(0);
        }

        // sr_tog is short range forces
        // es_tog is electrostatic solve (0=off, 1=Poisson, 2=Pairwise, 3=P3M)
        if(sr_tog==1 || es_tog==3)
        {
            // each tile clears its neighbors
            particles.clearNeighbors();
            // fill the neighbor buffers for each tile with the proper data
            particles.fillNeighbors();

            // compute short range forces (if sr_tog=1)
            // compute P3M short range correction (if es_tog=3)
            particles.computeForcesNL(charge, RealCenteredCoords, dxp);

        }

        if(es_tog==1 || es_tog==3)
        {
            // spreads charge density from ions onto multifab 'charge'.
            particles.collectFields(dt, dxp, RealCenteredCoords, geomP, charge, chargeTemp, massFrac, massFracTemp);
        }
        
        // do Poisson solve using 'charge' for RHS, and put potential in 'potential'. Then calculate gradient and put in 'efieldCC', then add 'external'.
        esSolve(potential, charge, efieldCC, external, geomP);


        // compute other forces and spread to grid
        particles.SpreadIons(dt, dx, dxp, geom, umac, efieldCC, charge, RealFaceCoords, RealCenteredCoords, source, sourceTemp, surfaceList,
                             surfaceCount, 3 /*this number currently does nothing, but we will use it later*/);

        if((variance_coef_mom != 0.0) && fluid_tog != 0) {
          // compute the random numbers needed for the stochastic momentum forcing
          sMflux.fillMStochastic();

          // compute stochastic momentum force
          sMflux.StochMomFluxDiv(stochMfluxdiv,0,eta_cc,eta_ed,temp_cc,temp_ed,weights,dt);

          // integrator containing inertial terms and predictor/corrector requires 2 RNG stages
          if(fluid_tog ==2) {
              sMflux.StochMomFluxDiv(stochMfluxdivC,0,eta_cc,eta_ed,temp_cc,temp_ed,weights,dt);
          }
        }

        // AJN - should this be an if/else fluid_tog==2? DRL - No.
    	advanceStokes(umac,pres,stochMfluxdiv,source,alpha_fc,beta,gamma,beta_ed,geom,dt);
        if(fluid_tog ==2) {
            advanceLowMach(umac, umacNew, pres, tracer, stochMfluxdiv, stochMfluxdivC, alpha_fc, beta, gamma, beta_ed, geom,dt);
        }

        // total particle move (1=single step, 2=midpoint)
        if(move_tog != 0)
        {
            //Calls wet ion interpolation and movement.
            Print() << "Start move.\n";
            particles.MoveIons(dt, dx, dxp, geom, umac, efield, RealFaceCoords, source, sourceTemp, dryMobility, surfaceList,
                               surfaceCount, 3 /*this number currently does nothing, but we will use it later*/);
            particles.Redistribute();
            particles.ReBin();
            Print() << "Finish move.\n";
        }

        //Start collecting statistics after step n_steps_skip
        if(step == n_steps_skip)
        {
            particleMeans.setVal(0.0);
            particleVars.setVal(0);
            AMREX_D_TERM(umacM[0].setVal(0);,
                         umacM[1].setVal(0);,
                         umacM[2].setVal(0););
            AMREX_D_TERM(umacV[0].setVal(0);,
                         umacV[1].setVal(0);,
                         umacV[2].setVal(0););

            Print() << "Resetting stat collection.\n";

            statsCount = 1;
        }
       
        particles.EvaluateStats(particleInstant, particleMeans, particleVars, cellVols, ionParticle[0], dt,statsCount);

        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            ComputeBasicStats(umac[d], umacM[d], umacV[d], 1, 1, statsCount);
        }

        statsCount++;
	//_______________________________________________________________________
	// Update structure factor




       // if(step > n_steps_skip && struct_fact_int > 0 && (step-n_steps_skip-1)%struct_fact_int == 0) {
//	      MultiFab::Copy(struct_in_cc, charge, 0, 0, nvar_sf, 0);
//	      structFact.FortStructure(struct_in_cc,geomP);
  //      }

        if (plot_int > 0 && step%plot_int == 0)
        {

            //This write particle data and associated fields and electrostatic fields
            WritePlotFile(step,time,geom,geomC,geomP,particleInstant, particleMeans, particleVars, particles, charge, potential, efieldCC, dryMobility);

            //Writes instantaneous flow field and some other stuff? Check with Guy.
            WritePlotFileHydro(step,time,geom,umac,pres, umacM, umacV);
        }

        if(step%1 == 0)
        {    
                amrex::Print() << "Advanced step " << step << "\n";
        }
        
        time = time + dt;

    }
    ///////////////////////////////////////////
    if (struct_fact_int > 0) {

        Real dVol = dx[0]*dx[1];
        int tot_n_cells = n_cells[0]*n_cells[1];
      if (AMREX_SPACEDIM == 2) {
	    dVol *= cell_depth;
      } else if (AMREX_SPACEDIM == 3) {
	    dVol *= dx[2];
	    tot_n_cells = n_cells[2]*tot_n_cells;
      }

      // let rho = 1
      //Real SFscale = dVol/(rho0*k_B*T_init[0]);
      Real SFscale = 1;
      // SFscale = 1.0;
      // Print() << "Hack: structure factor scaling = " << SFscale << std::endl;
      
    //  structFact.Finalize(SFscale);
    //  structFact.WritePlotFile(step,time,geomP,"plt_SF");

    }

    // Call the timer again and compute the maximum difference between the start time 
    // and stop time over all processors
    Real stop_time = ParallelDescriptor::second() - strt_time;
    ParallelDescriptor::ReduceRealMax(stop_time);
    amrex::Print() << "Run time = " << stop_time << std::endl;

}
