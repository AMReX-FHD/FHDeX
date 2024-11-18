#include "common_functions.H"
#include "compressible_functions.H"
#include "compressible_functions_stag.H"

#include <AMReX_Vector.H>
#include <AMReX_MPMD.H>

#include "rng_functions.H"

#include "StructFact.H"

#if defined(TURB)
#include "TurbForcingComp.H"
#endif

#include "chemistry_functions.H"

#include "MFsurfchem_functions.H"

#include "chrono"

#if defined(MUI) || defined(USE_AMREX_MPMD)
#include "surfchem_mui_functions.H"
using namespace surfchem_mui;
#endif

using namespace std::chrono;
using namespace amrex;

// argv contains the name of the inputs file entered at the command line
void main_driver(const char* argv)
{
    BL_PROFILE_VAR("main_driver()",main_driver);
    
    // store the current time so we can later compute total run time.
    Real strt_time = ParallelDescriptor::second();

    std::string inputs_file = argv;

    amrex::Print() << "Compiled with support for maximum species = " << MAX_SPECIES << "\n";
    
    // copy contents of F90 modules to C++ namespaces
    InitializeCommonNamespace();

    InitializeCompressibleNamespace();

    if (nvars != AMREX_SPACEDIM + 2 + nspecies) {
        Abort("nvars must be equal to AMREX_SPACEDIM + 2 + nspecies");
    }

    if (nprimvars != AMREX_SPACEDIM + 3 + 2*nspecies) {
        Abort("nprimvars must be equal to AMREX_SPACEDIM + 3 + 2*nspecies");
    }

    //if (advection_type != 2) {
    //    Abort("only interpolation of conserved quantities works for advective fluxes in the staggered code. this corresponds to advection_type = 2");
    //}

    // read the inputs file for chemistry
    InitializeChemistryNamespace();

    // read the inputs file for MFsurfchem
    InitializeMFSurfchemNamespace();

#if defined(MUI) || defined(USE_AMREX_MPMD)
    // read the inputs file for surfchem_mui
    InitializeSurfChemMUINamespace();

    if (n_ads_spec>0) {
        Abort("MFsurfchem cannot be used in compressible_stag_mui");
    }

    if (nspec_mui<1) {
        Abort("nspec_mui must be at least one");
    }

    if (restart>0) {
        Abort("restart not supported in compressible_stag_mui");
    }
#endif

    int step_start, statsCount;
    amrex::Real time;

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

    if (((do_1D) or (do_2D)) and (amrex::Math::abs(visc_type) == 3)) Abort("1D and 2D version only work for zero bulk viscosity currently. Use visc_type 1 or 2");

    if ((do_1D) and (do_2D)) Abort("Can not have both 1D and 2D mode on at the same time");

    // for each direction, if bc_vel_lo/hi is periodic, then
    // set the corresponding bc_mass_lo/hi and bc_therm_lo/hi to periodic
    SetupBCStag();

    // if multispecies
    if (algorithm_type == 2) {
        // compute wall concentrations if BCs call for it
        SetupCWallStag();
    }

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

#ifdef MUI
    // init MUI
    mui::uniface2d uniface( "mpi://FHD-side/FHD-KMC-coupling" );
#endif

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

    // MFsurfchem
    MultiFab surfcov;       // also used in surfchem_mui for stats and plotfiles
    MultiFab dNadsdes;

#if defined(MUI) || defined(USE_AMREX_MPMD)
    MultiFab Ntot;          // saves total number of sites
#endif

#if defined(USE_AMREX_MPMD)
    std::unique_ptr<MPMD::Copier> mpmd_copier;
#endif

    //statistics    
    MultiFab cuMeans;
    MultiFab cuVars;
    
    MultiFab primMeans;
    MultiFab primVars;

    MultiFab coVars;

    MultiFab surfcovMeans;  // used in either MFsurfchem or surfchem_mui
    MultiFab surfcovVars;   // used in either MFsurfchem or surfchem_mui
    MultiFab surfcovcoVars; // used in either MFsurfchem or surfchem_mui

    std::array< MultiFab, AMREX_SPACEDIM > velMeans;
    std::array< MultiFab, AMREX_SPACEDIM > velVars;
    std::array< MultiFab, AMREX_SPACEDIM > cumomMeans;
    std::array< MultiFab, AMREX_SPACEDIM > cumomVars;

    if ((plot_cross) and ((cross_cell < 0) or (cross_cell > n_cells[0]-1))) {
        Abort("Cross cell needs to be within the domain: 0 <= cross_cell <= n_cells[0] - 1");
    }
    if ((do_slab_sf) and ((membrane_cell <= 0) or (membrane_cell >= n_cells[0]-1))) {
        Abort("Slab structure factor needs a membrane cell within the domain: 0 < cross_cell < n_cells[0] - 1");
    }
    if ((project_dir >= 0) and ((do_1D) or (do_2D))) {
        Abort("Projected structure factors (project_dir) works only for 3D case");
    }
    if ((all_correl > 1) or (all_correl < 0)) {
        Abort("all_correl can be 0 or 1");
    }
    if ((all_correl == 1) and (cross_cell > 0) and (cross_cell < n_cells[0]-1)) {
        amrex::Print() << "Correlations will be done at four equi-distant x* because all_correl = 1" << "\n";
    }

    // contains yz-averaged running & instantaneous averages of conserved variables (2*nvars) + primitive variables [vx, vy, vz, T, Yk]: 2*4 + 2*nspecies 
    Vector<Real> dataSliceMeans_xcross(2*nvars+8+2*nspecies, 0.0); 
    
    // see statsStag for the list
    // can add more -- change main_driver, statsStag, writeplotfilestag, and Checkpoint
    int ncross = 37+nspecies+3;
    MultiFab spatialCross1D;
    MultiFab spatialCross2D;
    Vector<Real> spatialCross3D(n_cells[0]*ncross, 0.0);
    
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
    Real sys_volume = 1.0;
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        sys_volume *= (realDomain.hi(d) - realDomain.lo(d));
    }

    std::string filename = "crossMeans";
    std::ofstream outfile;

#if defined(TURB)
    // data structure for turbulence diagnostics
    std::string turbfilename = "turbstats";
    std::ofstream turboutfile;
    std::string turbfilenamedecomp = "turbstatsdecomp";
    std::ofstream turboutfiledecomp;
#endif

    /////////////////////////////////////////////
    // Setup Structure factor variables & scaling
    /////////////////////////////////////////////

    // Standard 3D structure factors
    StructFact structFactPrim;
    StructFact structFactCons;
    MultiFab structFactPrimMF;
    MultiFab structFactConsMF;

    // Structure factor for 2D averaged data
    StructFact structFactPrimVerticalAverage;
    StructFact structFactConsVerticalAverage;

    // Structure factor for 2D averaged data (across a membrane)
    StructFact structFactPrimVerticalAverage0;
    StructFact structFactPrimVerticalAverage1;
    StructFact structFactConsVerticalAverage0;
    StructFact structFactConsVerticalAverage1;
    MultiFab master_project_rot_prim;
    MultiFab master_project_rot_cons;

    // Vector of structure factors for 2D simulation
    Vector < StructFact > structFactPrimArray;
    Vector < StructFact > structFactConsArray;
    MultiFab master_2D_rot_prim;
    MultiFab master_2D_rot_cons;

#if defined(TURB)
    // Structure factor for compressible turbulence
    StructFact turbStructFactVelTotal; // total velocity
    StructFact turbStructFactVelDecomp; // decomposed velocity
    StructFact turbStructFactScalar; // scalars 
#endif

    // surface coverage structure factor
    StructFact surfcovStructFact;
    MultiFab surfcovFlattenedRotMaster;
    Geometry surfcov_geom_flat;
    BoxArray surfcov_ba_flat;
    DistributionMapping surfcov_dmap_flat;
    
    Geometry geom_flat;
    Geometry geom_flat_2D;
    BoxArray ba_flat;
    BoxArray ba_flat_2D;
    DistributionMapping dmap_flat;
    DistributionMapping dmap_flat_2D;

    // "primitive" variable structure factor will contain
    // rho
    // vel (shifted)
    // T
    // Yk
    // vel (averaged)
    // rhoYk (copy from cons)
    int structVarsPrim = 2*AMREX_SPACEDIM+2*nspecies+2;

    Vector< std::string > prim_var_names;
    prim_var_names.resize(structVarsPrim);

    int cnt = 0;
    int numvars;
    std::string x;

    // rho
    prim_var_names[cnt] = "rho";
    ++cnt;

    // velx, vely, velz
    for (int d=0; d<AMREX_SPACEDIM; d++) {
        x = "velCC";
        x += (120+d);
        prim_var_names[cnt] = x;
        ++cnt;
    }

    // Temp
    prim_var_names[cnt] = "Temp";
    ++cnt;

    // Yk
    for (int d=0; d<nspecies; d++) {
        x = "Y";
        x += (49+d);
        prim_var_names[cnt] = x;
        ++cnt;
    }

    // velx, vely, velz
    for (int d=0; d<AMREX_SPACEDIM; d++) {
        x = "velFACE";
        x += (120+d);
        prim_var_names[cnt] = x;
        ++cnt;
    }

    // rho*Yk
    for (int d=0; d<nspecies; d++) {
        x = "rhoY";
        x += (49+d);
        prim_var_names[cnt] = x;
        ++cnt;
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
    cons_var_names[cnt] = "rho";
    ++cnt;

    // jx, jy, jz
    for (int d=0; d<AMREX_SPACEDIM; d++) {
        x = "jCC";
        x += (120+d);
        cons_var_names[cnt] = x;
        ++cnt;
    }

    // rho*E
    cons_var_names[cnt] = "rhoE";
    ++cnt;

    // rho*Yk
    for (int d=0; d<nspecies; d++) {
        x = "rhoY";
        x += (49+d);
        cons_var_names[cnt] = x;
        ++cnt;
    }

    // Temp
    cons_var_names[cnt] = "Temp";
    ++cnt;

    // jx, jy, jz
    for (int d=0; d<AMREX_SPACEDIM; d++) {
        x = "jFACE";
        x += (120+d);
        cons_var_names[cnt] = x;
        ++cnt;
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

#if defined(TURB)
    //////////////////////////////////////////////////////////////
    // structure factor variables names and scaling for turbulence
    // variables are velocities, density, pressure and temperature
    //////////////////////////////////////////////////////////////
    // need to use dVol for scaling
    Real dVol = (AMREX_SPACEDIM==2) ? dx[0]*dx[1]*cell_depth : dx[0]*dx[1]*dx[2];
    Real dVolinv = 1.0/dVol;
    
    MultiFab structFactMFTurbVel;
    MultiFab structFactMFTurbScalar;
    MultiFab vel_decomp;

    Vector< std::string > var_names_turbVelTotal{"ux","uy","uz"};
    Vector<Real> var_scaling_turbVelTotal(3, dVolinv);
    amrex::Vector< int > s_pairA_turbVelTotal(3);
    amrex::Vector< int > s_pairB_turbVelTotal(3);
    for (int d=0; d<3; ++d) {
        s_pairA_turbVelTotal[d] = d;
        s_pairB_turbVelTotal[d] = d;
    }
    
    Vector<Real> var_scaling_turbVelDecomp(6, dVolinv);
    
    Vector< std::string > var_names_turbScalar{"rho","tenp","press"};
    Vector<Real> var_scaling_turbScalar(3, dVolinv);
    amrex::Vector< int > s_pairA_turbScalar(3);
    amrex::Vector< int > s_pairB_turbScalar(3);
    for (int d=0; d<3; ++d) {
        s_pairA_turbScalar[d] = d;
        s_pairB_turbScalar[d] = d;
    }
#endif
    
    // object for turbulence forcing
    TurbForcingComp turbforce;

    /////////////////////////////////////////////
    // Initialize based on fresh start or restart
    /////////////////////////////////////////////

    if (restart > 0) {
        
        if (do_1D) {
            ReadCheckPoint1D(step_start, time, statsCount, geom, domain, cu, cuMeans, cuVars, prim,
                             primMeans, primVars, cumom, cumomMeans, cumomVars, 
                             vel, velMeans, velVars, coVars, spatialCross1D, ncross, ba, dmap);
        }
        else if (do_2D) {
            ReadCheckPoint2D(step_start, time, statsCount, geom, domain, cu, cuMeans, cuVars, prim,
                             primMeans, primVars, cumom, cumomMeans, cumomVars, 
                             vel, velMeans, velVars, coVars, spatialCross2D, ncross, ba, dmap);
        }
        else {
            ReadCheckPoint3D(step_start, time, statsCount, geom, domain, cu, cuMeans, cuVars, prim,
                             primMeans, primVars, cumom, cumomMeans, cumomVars, 
                             vel, velMeans, velVars, coVars, surfcov, surfcovMeans, surfcovVars, surfcovcoVars, spatialCross3D, ncross, turbforce, ba, dmap);
        }

        if (reset_stats == 1) statsCount = 1;

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
        
        if (n_ads_spec>0) {
            dNadsdes.define(ba,dmap,n_ads_spec,0);
            nspec_surfcov = n_ads_spec;
        }

        if ((plot_cross) and (do_1D==0) and (do_2D==0)) {
            if (ParallelDescriptor::IOProcessor()) outfile.open(filename, std::ios::app);
        }

#if defined(TURB)
        if (turbForcing >= 1) { // temporary fab for turbulent
            if (ParallelDescriptor::IOProcessor()) {
                turboutfile.open(turbfilename, std::ios::app);
                turboutfiledecomp.open(turbfilenamedecomp, std::ios::app);
            }
        }
#endif
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

        if (n_ads_spec>0) {
            surfcov.define(ba,dmap,n_ads_spec,0);
            dNadsdes.define(ba,dmap,n_ads_spec,0);
            nspec_surfcov = n_ads_spec;
        }

#if defined(MUI) || defined(USE_AMREX_MPMD)
        surfcov.define(ba,dmap,nspec_mui,0);
        Ntot.define(ba,dmap,1,0);
        nspec_surfcov = nspec_mui;
#endif

        cuMeans.define(ba,dmap,nvars,ngc);
        cuVars.define(ba,dmap,nvars,ngc);
        cuMeans.setVal(0.0);
        cuVars.setVal(0.0);
        
        primMeans.define(ba,dmap,nprimvars+3,ngc); // the last three have COM velocity
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
        // 10: <rhoYkL rhoYkH>
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
        // 21: <YkH YkL>
        // 22: <YkL velx>
        // 23: <YkH velx>
        // 24: <rhoYkL velx>
        // 25: <rhoYkH velx>
        coVars.define(ba,dmap,26,0);
        coVars.setVal(0.0);

        if (nspec_surfcov>0) {
            surfcovMeans.define(ba,dmap,nspec_surfcov,0);
            surfcovVars.define(ba,dmap,nspec_surfcov,0);
	    surfcovcoVars.define(ba,dmap,nspec_surfcov*6,0);
            surfcovMeans.setVal(0.0);
            surfcovVars.setVal(0.0);
	    surfcovcoVars.setVal(0.0);
        }

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

        if (do_1D) {
            if (all_correl) spatialCross1D.define(ba,dmap,ncross*5,0); // for five x*: [0, fl(n_cells[0]/4), fl(n_cells[0]/2), fl(n_cells[0]*3/4), n_cells[0]-1]
            else spatialCross1D.define(ba,dmap,ncross,0);
            spatialCross1D.setVal(0.0);
        }
        else if (do_2D) {
            spatialCross2D.define(ba,dmap,ncross,0);
            spatialCross2D.setVal(0.0);
        }

#if defined(TURB)
        if (turbForcing >= 1) { // temporary fab for turbulent
            if (ParallelDescriptor::IOProcessor()) {
                turboutfile.open(turbfilename);
                turboutfile << "step " << "time " << "turbKE " << "RMSu " 
                            << "<c> " << "TaylorLen " << "TaylorRe " << "TaylorMa "
                            << "skew " << "kurt "
                            << "eps_s " << "eps_d " << "eps_d/eps_s "
                            << "kolm_s " << "kolm_s" << "kolm_t"
                            << std::endl;

                turboutfiledecomp.open(turbfilenamedecomp);
                turboutfiledecomp << "step " << "time " 
                                  << "turbKE_s " << "turbKE_d " << "delta_turbKE "
                                  << "u_rms_s " << "u_rms_d " << "delta_u_rms " 
                                  << "TaylorMa_d "
                                  << "skew_s " << "kurt_s "
                                  << "skew_d " << "kurt_d "
                                  << std::endl;

            }
        }
#endif

        ///////////////////////////////////////////
        // Initialize everything
        ///////////////////////////////////////////

        // initialize conserved variables
        InitConsVarStag(cu,cumom,geom); // Need to add for staggered -- Ishan

        // initialize primitive variables
        //prim.setVal(0.0,0,nprimvars,ngc);
        //for (int d=0; d<AMREX_SPACEDIM; d++) { // staggered momentum & velocities
        //    vel[d].setVal(0.,ngc);
        //}
        conservedToPrimitiveStag(prim, vel, cu, cumom);

        if (n_ads_spec>0) init_surfcov(surfcov, geom);

#if defined(MUI)
        mui_fetch_Ntot(Ntot, dx, uniface, 0);

        mui_fetch_surfcov(Ntot, surfcov, dx, uniface, 0);

        mui_forget(uniface, 0);

#elif defined(USE_AMREX_MPMD)
        mpmd_copier = std::make_unique<MPMD::Copier>(Ntot.boxArray(),
                                                     Ntot.DistributionMap());
        amrex_fetch_Ntot(Ntot, *mpmd_copier);
        amrex_fetch_surfcov(Ntot, surfcov, *mpmd_copier);
#endif

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
                          prim, primMeans, primVars, vel, velMeans, velVars, coVars, surfcov, surfcovMeans, surfcovVars, surfcovcoVars, eta, kappa, zeta);
#if defined(TURB)
            if (turbForcing > 0) {
                EvaluateWritePlotFileVelGrad(0, 0.0, geom, vel, vel_decomp);
            }
#endif

            if (plot_cross) {
                if (do_1D) {
                    WriteSpatialCross1D(spatialCross1D, 0, geom, ncross);
                }
                else if (do_2D) {
                //    WriteSpatialCross2D(spatialCross2D, 0, geom, ncross); // (do later)
                }
                else {
                    WriteSpatialCross3D(spatialCross3D, 0, geom, ncross);
                }
            }
        }

        if ((plot_cross) and (do_1D==0) and (do_2D==0)) {
            if (ParallelDescriptor::IOProcessor()) outfile.open(filename);
        }

        step_start = 1;
        time = 0.;
        statsCount = 1;

#if defined(TURB)
        // define turbulence forcing object
        if (turbForcing > 1) {
            turbforce.define(ba,dmap,turb_a,turb_b,turb_c,turb_d,turb_alpha);
        }
#endif


    } // end t=0 setup

    ///////////////////////////////////////////
    // Setup Structure factor
    ///////////////////////////////////////////

    if (struct_fact_int > 0) {
        structFactPrimMF.define(ba, dmap, structVarsPrim, 0);
        structFactPrim.define(ba,dmap,prim_var_names,var_scaling_prim);
            
        structFactConsMF.define(ba, dmap, structVarsCons, 0);
        structFactCons.define(ba,dmap,cons_var_names,var_scaling_cons);
        
        // structure factor class for vertically-averaged dataset
        if (project_dir >= 0) {

            {
                MultiFab X, XRot;
                ComputeVerticalAverage(prim, X, geom, project_dir, 0, nprimvars);
                XRot = RotateFlattenedMF(X);
                ba_flat = XRot.boxArray();
                dmap_flat = XRot.DistributionMap();
                master_project_rot_prim.define(ba_flat,dmap_flat,structVarsPrim,0);
                master_project_rot_cons.define(ba_flat,dmap_flat,structVarsCons,0);

                IntVect dom_lo_flat(AMREX_D_DECL(0,0,0));
                IntVect dom_hi_flat;
#if (AMREX_SPACEDIM == 2)
                if (project_dir == 0) {
                    dom_hi_flat[0] = n_cells[1]-1;
                    dom_hi_flat[1] = 0;
                }
                else if (project_dir == 1) {
                    dom_hi_flat[0] = n_cells[0]-1;
                    dom_hi_flat[1] = 0;
                }
#elif (AMREX_SPACEDIM == 3)
                if (project_dir == 0) {
                    dom_hi_flat[0] = n_cells[1]-1;
                    dom_hi_flat[1] = n_cells[2]-1;
                    dom_hi_flat[2] = 0;
                } else if (project_dir == 1) {
                    dom_hi_flat[0] = n_cells[0]-1;
                    dom_hi_flat[1] = n_cells[2]-1;
                    dom_hi_flat[2] = 0;
                } else if (project_dir == 2) {
                    dom_hi_flat[0] = n_cells[0]-1;
                    dom_hi_flat[1] = n_cells[1]-1;
                    dom_hi_flat[2] = 0;
                }
#endif
                Box domain_flat(dom_lo_flat, dom_hi_flat);

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

                RealBox real_box_flat({AMREX_D_DECL(     prob_lo[0],     prob_lo[1],     prob_lo[2])},
                                    {AMREX_D_DECL(projected_hi[0],projected_hi[1],projected_hi[2])});
          
                // This defines a Geometry object
                geom_flat.define(domain_flat,&real_box_flat,CoordSys::cartesian,is_periodic.data());

            }

            if (do_slab_sf == 0) {
                structFactPrimVerticalAverage.define(ba_flat,dmap_flat,prim_var_names,var_scaling_prim,2);
                structFactConsVerticalAverage.define(ba_flat,dmap_flat,cons_var_names,var_scaling_cons,2);
            }
            else {
                structFactPrimVerticalAverage0.define(ba_flat,dmap_flat,prim_var_names,var_scaling_prim);
                structFactPrimVerticalAverage1.define(ba_flat,dmap_flat,prim_var_names,var_scaling_prim);
                structFactConsVerticalAverage0.define(ba_flat,dmap_flat,cons_var_names,var_scaling_cons);
                structFactConsVerticalAverage1.define(ba_flat,dmap_flat,cons_var_names,var_scaling_cons);
            }
    
        }

        if (do_2D) { // 2D is coded only for XY plane

            {
                MultiFab X, XRot;
                ExtractSlice(prim, X, geom, 2, 0, 0, nprimvars);
                XRot = RotateFlattenedMF(X);
                ba_flat_2D = XRot.boxArray();
                dmap_flat_2D = XRot.DistributionMap();
                master_2D_rot_prim.define(ba_flat_2D,dmap_flat_2D,structVarsPrim,0);
                master_2D_rot_cons.define(ba_flat_2D,dmap_flat_2D,structVarsCons,0);

                IntVect dom_lo_flat(AMREX_D_DECL(0,0,0));
                IntVect dom_hi_flat;
                dom_hi_flat[0] = n_cells[0]-1;
                dom_hi_flat[1] = n_cells[1]-1;
                dom_hi_flat[2] = 0;
                Box domain_flat(dom_lo_flat, dom_hi_flat);

                // This defines the physical box
                Vector<Real> projected_hi(AMREX_SPACEDIM);
                for (int d=0; d<AMREX_SPACEDIM; d++) {
                    projected_hi[d] = prob_hi[d];
                }
                projected_hi[AMREX_SPACEDIM-1] = prob_hi[2] / n_cells[2];

                RealBox real_box_flat({AMREX_D_DECL(     prob_lo[0],     prob_lo[1],     prob_lo[2])},
                                      {AMREX_D_DECL(projected_hi[0],projected_hi[1],projected_hi[2])});
          
                // This defines a Geometry object
                geom_flat_2D.define(domain_flat,&real_box_flat,CoordSys::cartesian,is_periodic.data());

            }

            structFactPrimArray.resize(n_cells[2]);
            structFactConsArray.resize(n_cells[2]);

            for (int i = 0; i < n_cells[2]; ++i) { 
                structFactPrimArray[i].define(ba_flat_2D,dmap_flat_2D,prim_var_names,var_scaling_prim,2);
                structFactConsArray[i].define(ba_flat_2D,dmap_flat_2D,cons_var_names,var_scaling_cons,2);
            }

        }
    }
    
#if defined(TURB)
    if (turbForcing >= 1) {
        
        structFactMFTurbVel.define(ba, dmap, 3, 0);
        structFactMFTurbScalar.define(ba, dmap, 6, 0);
        vel_decomp.define(ba, dmap, 6, 0);
        vel_decomp.setVal(0.0);

        turbStructFactVelTotal.define(ba,dmap,
                var_names_turbVelTotal,var_scaling_turbVelTotal,
                s_pairA_turbVelTotal,s_pairB_turbVelTotal);
        turbStructFactScalar.define(ba,dmap,
                var_names_turbScalar,var_scaling_turbScalar,
                s_pairA_turbScalar,s_pairB_turbScalar);
        turbStructFactVelDecomp.defineDecomp(ba,dmap,
                var_names_turbVelTotal,var_scaling_turbVelDecomp,
                s_pairA_turbVelTotal,s_pairB_turbVelTotal);
    }
#endif

    if (n_ads_spec>0) {

        MultiFab Flattened;  // flattened multifab defined below

        // we are only calling ExtractSlice here to obtain
        // a built version of Flattened so can obtain what we need to build the
        // structure factor and geometry objects for flattened data
        // assume surface covered is stored in the "k" direction in the k=0 coordinate.
        int surfcov_dir = 2;
        int surfcov_plane = 0;
        int surfcov_structVars = n_ads_spec;
        int surfcov_nPairs = surfcov_structVars*(surfcov_structVars+1)/2;

        Vector< std::string > surfcov_var_names;
        surfcov_var_names.resize(surfcov_structVars);
        for (int d=0; d<surfcov_structVars; d++) {
            x = "surfCov";
            x += (48+d);
            surfcov_var_names[d] = x;
        }

        Vector<Real> surfcov_var_scaling(surfcov_nPairs);
        for (int d=0; d<surfcov_var_scaling.size(); ++d) {
            surfcov_var_scaling[d] = 1.;
        }
      
        ExtractSlice(surfcov, Flattened, geom, surfcov_dir, surfcov_plane, 0, surfcov_structVars);
        // we rotate this flattened MultiFab to have normal in the z-direction since
        // our structure factor class assumes this for flattened
        MultiFab FlattenedRot = RotateFlattenedMF(Flattened);
        BoxArray surfcov_surfcov_ba_flat = FlattenedRot.boxArray();
        const DistributionMapping& dmap_flat = FlattenedRot.DistributionMap();
        surfcovFlattenedRotMaster.define(surfcov_surfcov_ba_flat,dmap_flat,surfcov_structVars,0);
        {
            IntVect dom_lo(AMREX_D_DECL(0,0,0));
            IntVect dom_hi;

            // yes you could simplify this code but for now
            // these are written out fully to better understand what is happening
            // we wanted dom_hi[AMREX_SPACEDIM-1] to be equal to 0
            // and need to transmute the other indices depending on surfcov_dir
#if (AMREX_SPACEDIM == 2)
            if (surfcov_dir == 0) {
                dom_hi[0] = n_cells[1]-1;
            }
            else if (surfcov_dir == 1) {
                dom_hi[0] = n_cells[0]-1;
            }
            dom_hi[1] = 0;
#elif (AMREX_SPACEDIM == 3)
            if (surfcov_dir == 0) {
                dom_hi[0] = n_cells[1]-1;
                dom_hi[1] = n_cells[2]-1;
            } else if (surfcov_dir == 1) {
                dom_hi[0] = n_cells[0]-1;
                dom_hi[1] = n_cells[2]-1;
            } else if (surfcov_dir == 2) {
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
            // and need to transmute the other indices depending on surfcov_dir
#if (AMREX_SPACEDIM == 2)
            if (surfcov_dir == 0) {
                projected_hi[0] = prob_hi[1];
            } else if (surfcov_dir == 1) {
                projected_hi[0] = prob_hi[0];
            }
            projected_hi[1] = prob_hi[surfcov_dir] / n_cells[surfcov_dir];
#elif (AMREX_SPACEDIM == 3)
            if (surfcov_dir == 0) {
                projected_hi[0] = prob_hi[1];
                projected_hi[1] = prob_hi[2];
            } else if (surfcov_dir == 1) {
                projected_hi[0] = prob_hi[0];
                projected_hi[1] = prob_hi[2];
            } else if (surfcov_dir == 2) {
                projected_hi[0] = prob_hi[0];
                projected_hi[1] = prob_hi[1];
            }
            projected_hi[2] = prob_hi[surfcov_dir] / n_cells[surfcov_dir];
#endif

            RealBox real_box({AMREX_D_DECL(     prob_lo[0],     prob_lo[1],     prob_lo[2])},
                             {AMREX_D_DECL(projected_hi[0],projected_hi[1],projected_hi[2])});
        
            // This defines a Geometry object
            surfcov_geom_flat.define(domain,&real_box,CoordSys::cartesian,is_periodic.data());
        }

        surfcovStructFact.define(surfcov_surfcov_ba_flat,dmap_flat,surfcov_var_names,surfcov_var_scaling);
    }

    /////////////////////////////////////////////////
    // Initialize Fluxes and Sources
    /////////////////////////////////////////////////
    
    // external source term - possibly for later
    MultiFab source(ba,dmap,nprimvars,ngc);
    source.setVal(0.0);

    MultiFab ranchem;
    if (nreaction>0) ranchem.define(ba,dmap,nreaction,ngc);

    //fluxes (except momentum) at faces
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
    std::array< MultiFab, AMREX_SPACEDIM > faceflux;
    AMREX_D_TERM(faceflux[0].define(convert(ba,nodal_flag_x), dmap, nvars+4, 0);,
                 faceflux[1].define(convert(ba,nodal_flag_y), dmap, nvars+4, 0);,
                 faceflux[2].define(convert(ba,nodal_flag_z), dmap, nvars+4, 0););

    //momentum flux (edge + center)
#if (AMREX_SPACEDIM == 3)
    std::array< MultiFab, 2 > edgeflux_x;
    std::array< MultiFab, 2 > edgeflux_y;
    std::array< MultiFab, 2 > edgeflux_z;

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
    

#if defined(TURB)
    // Initialize Turbulence Forcing Object
    if (turbForcing > 1) {
        turbforce.Initialize(geom);
    }
#endif
                
    /////////////////////////////////////////////////
    //Time stepping loop
    /////////////////////////////////////////////////
    
    for (int step=step_start;step<=max_step;++step) {

        // timer
        Real ts1 = ParallelDescriptor::second();

        // sample surface chemistry
#if defined(MUI)
        mui_push(cu, prim, dx, uniface, step);

        mui_commit(uniface, step);
#elif defined(USE_AMREX_MPMD)
        amrex_push(cu, prim, *mpmd_copier);
#endif
        if (n_ads_spec>0) {
	    if (splitting_MFsurfchem == 0) sample_MFsurfchem(cu, prim, surfcov, dNadsdes, geom, dt);
	    else if (splitting_MFsurfchem == 1) {
	        sample_MFsurfchem(cu, prim, surfcov, dNadsdes, geom, dt/2.0);
		update_MFsurfchem(cu, prim, surfcov, dNadsdes, geom);

		for (int d=0; d<AMREX_SPACEDIM; d++) {
		    cumom[d].FillBoundary(geom.periodicity());
		}
		cu.FillBoundary(geom.periodicity());

		conservedToPrimitiveStag(prim, vel, cu, cumom);
		            // Set BC: 1) fill boundary 2) physical
                for (int d=0; d<AMREX_SPACEDIM; d++) {
                    vel[d].FillBoundary(geom.periodicity());
                }
                prim.FillBoundary(geom.periodicity());
                cu.FillBoundary(geom.periodicity());

                setBCStag(prim, cu, cumom, vel, geom);
	    }
	    else Abort("splitting_MFsurfchem can be 0 or 1");
	}

        // FHD
        RK3stepStag(cu, cumom, prim, vel, source, eta, zeta, kappa, chi, D, 
            faceflux, edgeflux_x, edgeflux_y, edgeflux_z, cenflux, ranchem, geom, dt, step, turbforce);

	if (n_ads_spec>0 && splitting_MFsurfchem == 1) sample_MFsurfchem(cu, prim, surfcov, dNadsdes, geom, dt/2.0);

        // update surface chemistry (via either surfchem_mui or MFsurfchem)
#if defined(MUI) || defined(USE_AMREX_MPMD)
#if defined(MUI)
        mui_fetch(cu, prim, dx, uniface, step);

        mui_fetch_surfcov(Ntot, surfcov, dx, uniface, step);

        mui_forget(uniface, step);
#elif defined(USE_AMREX_MPMD)
        amrex_fetch(cu, prim, geom.CellSizeArray(), *mpmd_copier);
        amrex_fetch_surfcov(Ntot, surfcov, *mpmd_copier);
#endif

        for (int d=0; d<AMREX_SPACEDIM; d++) {
            cumom[d].FillBoundary(geom.periodicity());
        }
        cu.FillBoundary(geom.periodicity());

        conservedToPrimitiveStag(prim, vel, cu, cumom);

        // Set BC: 1) fill boundary 2) physical
        for (int d=0; d<AMREX_SPACEDIM; d++) {
            vel[d].FillBoundary(geom.periodicity());
        }
        prim.FillBoundary(geom.periodicity());
        cu.FillBoundary(geom.periodicity());

        setBCStag(prim, cu, cumom, vel, geom);
#endif

        if (n_ads_spec>0) {

            update_MFsurfchem(cu, prim, surfcov, dNadsdes, geom);

            for (int d=0; d<AMREX_SPACEDIM; d++) {
                cumom[d].FillBoundary(geom.periodicity());
            }
            cu.FillBoundary(geom.periodicity());

            conservedToPrimitiveStag(prim, vel, cu, cumom);

            // Set BC: 1) fill boundary 2) physical
            for (int d=0; d<AMREX_SPACEDIM; d++) {
                vel[d].FillBoundary(geom.periodicity());
            }
            prim.FillBoundary(geom.periodicity());
            cu.FillBoundary(geom.periodicity());

            setBCStag(prim, cu, cumom, vel, geom);
        }

        // timer
        Real ts2 = ParallelDescriptor::second() - ts1;
        ParallelDescriptor::ReduceRealMax(ts2, ParallelDescriptor::IOProcessorNumber());
        if (step%100 == 0) {
            amrex::Print() << "Advanced step " << step << " in " << ts2 << " seconds\n";
        }

        // timer
        Real aux1 = ParallelDescriptor::second();
        
        // reset statistics after n_steps_skip
        // if n_steps_skip is negative, we use it as an interval
        if ((n_steps_skip > 0 && step == n_steps_skip) ||
            (n_steps_skip < 0 && step%amrex::Math::abs(n_steps_skip) == 0) ) {

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

            if (nspec_surfcov>0) {
                surfcovMeans.setVal(0.0);
                surfcovVars.setVal(0.0);
		surfcovcoVars.setVal(0.0);
            }

            if (do_1D) {
                spatialCross1D.setVal(0.0);
            }
            else if (do_2D) {
                spatialCross2D.setVal(0.0);
            }
            else {
                spatialCross3D.assign(spatialCross3D.size(), 0.0);
            }

            std::printf("Resetting stat collection.\n");

            statsCount = 1;

        }

        // Evaluate Statistics
        if (do_1D) {
            evaluateStatsStag1D(cu, cuMeans, cuVars, prim, primMeans, primVars, vel, 
                                velMeans, velVars, cumom, cumomMeans, cumomVars, coVars,
                                surfcov, surfcovMeans, surfcovVars, surfcovcoVars,
                                spatialCross1D, ncross, statsCount, geom);
        }
        else if (do_2D) {
            evaluateStatsStag2D(cu, cuMeans, cuVars, prim, primMeans, primVars, vel, 
                                velMeans, velVars, cumom, cumomMeans, cumomVars, coVars,
                                surfcov, surfcovMeans, surfcovVars, surfcovcoVars,
                                spatialCross2D, ncross, statsCount, geom);
        }
        else {
            evaluateStatsStag3D(cu, cuMeans, cuVars, prim, primMeans, primVars, vel, 
                                velMeans, velVars, cumom, cumomMeans, cumomVars, coVars,
                                surfcov, surfcovMeans, surfcovVars, surfcovcoVars,
                                dataSliceMeans_xcross, spatialCross3D, ncross, domain,
                                statsCount, geom);
        }
        statsCount++;
        if (step%100 == 0) {
            amrex::Print() << "Mean Rho: "      << ComputeSpatialMean(cu, 0) 
                           << " Mean Temp.:"    << ComputeSpatialMean(prim, 4)
                           << " Mean Press.:"   << ComputeSpatialMean(prim, 5)
                           << " Mean Mom. (x):" << ComputeSpatialMean(cumom[0], 0) 
                           << " Mean Mom. (y):" << ComputeSpatialMean(cumom[1], 0) 
                           << " Mean Mom. (z):" << ComputeSpatialMean(cumom[2], 0) 
                           << " Mean En.:"      << ComputeSpatialMean(cu, 4) 
                           << "\n";
        }


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
                              prim, primMeans, primVars, vel, velMeans, velVars, coVars, surfcov, surfcovMeans, surfcovVars, surfcovcoVars, eta, kappa, zeta);
            
#if defined(TURB)
            if (turbForcing > 0) {
                EvaluateWritePlotFileVelGrad(step, time, geom, vel, vel_decomp);
            }
#endif

            if (plot_cross) {
                if (do_1D) {
                    WriteSpatialCross1D(spatialCross1D, step, geom, ncross);
                }
                else if (do_2D) {
                //    WriteSpatialCross2D(spatialCross2D, step, geom, ncross); // (do later)
                }
                else {
                    WriteSpatialCross3D(spatialCross3D, step, geom, ncross);
                    if (ParallelDescriptor::IOProcessor()) {
                        outfile << step << " ";
                        for (auto l=0; l<2*nvars+8+2*nspecies; ++l) {
                            outfile << dataSliceMeans_xcross[l] << " ";
                        }
                        outfile << std::endl;
                    }
                }
            }

#if defined(TURB)
            // compressible turbulence structure factor snapshot
            if (turbForcing >= 1) {

                // copy velocities into structFactMFTurb
                for(int d=0; d<AMREX_SPACEDIM; d++) {
                    ShiftFaceToCC(vel[d], 0, structFactMFTurbVel, d, 1);
                }
                MultiFab::Copy(structFactMFTurbScalar, prim, 0, 0, 1, 0);
                MultiFab::Copy(structFactMFTurbScalar, prim, 4, 1, 1, 0);
                MultiFab::Copy(structFactMFTurbScalar, prim, 5, 2, 1, 0);
                
                 // decomposed velocities
                    turbStructFactVelDecomp.FortStructureDecomp(structFactMFTurbVel,geom,1);
                    turbStructFactVelDecomp.GetDecompVel(vel_decomp,geom);
                    turbStructFactVelDecomp.CallFinalize(geom);
                    turbStructFactVelDecomp.IntegratekShellsDecomp(step,geom,"vel_solenoid","vel_dilation");
                
                 // total velocity
                    turbStructFactVelTotal.FortStructure(structFactMFTurbVel,geom,1);
                    turbStructFactVelTotal.CallFinalize(geom);
                    turbStructFactVelTotal.IntegratekShells(step,geom,"vel_total");
                
                 // scalars
                    turbStructFactScalar.FortStructure(structFactMFTurbScalar,geom,1);
                    turbStructFactScalar.CallFinalize(geom);
                    turbStructFactScalar.IntegratekShellsScalar(step,geom,var_names_turbScalar);
            }
#endif
        }


#if defined(TURB)
        // turbulence outputs
        if ((turbForcing >= 1) and (step%1000 == 0)) {

            Real turbKE, c_speed, u_rms, taylor_len, taylor_Re, taylor_Ma,
            skew, kurt, eps_s, eps_d, eps_ratio, kolm_s, kolm_d, kolm_t;
            for (int i=0; i<AMREX_SPACEDIM; ++i) {
                vel[i].FillBoundary(geom.periodicity());
                cumom[i].FillBoundary(geom.periodicity());
            }
            GetTurbQty(vel, cumom, prim, eta, geom,
                       turbKE, c_speed, u_rms,
                       taylor_len, taylor_Re, taylor_Ma,
                       skew, kurt,
                       eps_s, eps_d, eps_ratio,
                       kolm_s, kolm_d, kolm_t);
            
            turboutfile << step << " ";
            turboutfile << time << " ";
            turboutfile << turbKE << " ";
            turboutfile << u_rms << " ";
            turboutfile << c_speed << " ";
            turboutfile << taylor_len << " " ;
            turboutfile << taylor_Re << " "; 
            turboutfile << taylor_Ma << " ";
            turboutfile << skew << " ";
            turboutfile << kurt << " ";
            turboutfile << eps_s << " ";
            turboutfile << eps_d << " "; 
            turboutfile << eps_ratio << " ";
            turboutfile << kolm_s << " ";
            turboutfile << kolm_d << " ";
            turboutfile << kolm_t;
            turboutfile << std::endl;
        }
        
        if ((turbForcing >= 1) and (writePlt)) {
            Real turbKE_s, turbKE_d, delta_turbKE;
            Real u_rms_s, u_rms_d, delta_u_rms;
            Real taylor_Ma_d;
            Real skew_s, kurt_s;
            Real skew_d, kurt_d;
            GetTurbQtyDecomp(vel_decomp, prim, geom,
                             turbKE_s, turbKE_d, delta_turbKE,
                             u_rms_s, u_rms_d, delta_u_rms,
                             taylor_Ma_d,
                             skew_s, kurt_s,
                             skew_d, kurt_d);
            
            turboutfiledecomp << step << " ";
            turboutfiledecomp << time << " ";
            turboutfiledecomp << turbKE_s << " ";
            turboutfiledecomp << turbKE_d << " ";
            turboutfiledecomp << delta_turbKE << " ";
            turboutfiledecomp << u_rms_s << " ";
            turboutfiledecomp << u_rms_d << " ";
            turboutfiledecomp << delta_u_rms << " ";
            turboutfiledecomp << taylor_Ma_d << " ";
            turboutfiledecomp << skew_s << " ";
            turboutfiledecomp << kurt_s << " ";
            turboutfiledecomp << skew_d << " ";
            turboutfiledecomp << kurt_d;
            turboutfiledecomp << std::endl;
        }
#endif
        
        // collect a snapshot for structure factor
        if (step > amrex::Math::abs(n_steps_skip) && 
            struct_fact_int > 0 && 
            (step-amrex::Math::abs(n_steps_skip))%struct_fact_int == 0) {
            
            /////////// First structFactPrimMF ////////////////
            cnt = 0;
            
            // copy [rho, vx, vy, vz, T]
            numvars = 5;
            MultiFab::Copy(structFactPrimMF, prim, 0, cnt, numvars, 0);
            cnt+=numvars;

            // copy Yk
            numvars = nspecies;
            MultiFab::Copy(structFactPrimMF, prim, AMREX_SPACEDIM+3, cnt, numvars, 0);
            cnt+=numvars;

            // copy velFACE
            for (int d=0; d<AMREX_SPACEDIM; ++d) {
                ShiftFaceToCC(vel[d],0,structFactPrimMF,cnt,1);
                ++cnt;
            }

            // copy rhoYk
            numvars = nspecies;
            MultiFab::Copy(structFactPrimMF, cu, AMREX_SPACEDIM+2, cnt, numvars, 0);
            ////////////////////////////////////////////////////

            ////////////// Second structFactConsMF /////////////
            cnt = 0;

            // copy [rho, jx, jy, jz, rhoE, rhoYk]
            numvars = nvars;
            MultiFab::Copy(structFactConsMF, cu, 0, cnt, numvars, 0);
            cnt+=numvars;

            // T
            numvars = 1;
            MultiFab::Copy(structFactConsMF, prim, AMREX_SPACEDIM+1, cnt, numvars, 0);
            cnt+=numvars;

            // copy jxFACE
            for (int d=0; d<AMREX_SPACEDIM; ++d) {
                ShiftFaceToCC(cumom[d],0,structFactConsMF,cnt,1);
                ++cnt;
            }
            ////////////////////////////////////////////////////

            if ((do_1D==0) and (do_2D==0)) {
                structFactPrim.FortStructure(structFactPrimMF,geom);
                structFactCons.FortStructure(structFactConsMF,geom);
            }

            if (project_dir >= 0) {

                if (do_slab_sf == 0) {
                    
                    {
                        MultiFab X, XRot;

                        ComputeVerticalAverage(structFactPrimMF, X, geom, project_dir, 0, structVarsPrim);
                        XRot = RotateFlattenedMF(X);
                        master_project_rot_prim.ParallelCopy(XRot, 0, 0, structVarsPrim); 
                        structFactPrimVerticalAverage.FortStructure(master_project_rot_prim,geom_flat);
                    }

                    {
                        MultiFab X, XRot;

                        ComputeVerticalAverage(structFactConsMF, X, geom, project_dir, 0, structVarsCons);
                        XRot = RotateFlattenedMF(X);
                        master_project_rot_cons.ParallelCopy(XRot, 0, 0, structVarsCons);
                        structFactConsVerticalAverage.FortStructure(master_project_rot_cons,geom_flat);
                    }

                }
                else {
                    
                    {
                        MultiFab X, XRot;

                        ComputeVerticalAverage(structFactPrimMF, X, geom, project_dir, 0, structVarsPrim, 0, membrane_cell-1);
                        XRot = RotateFlattenedMF(X);
                        master_project_rot_prim.ParallelCopy(XRot, 0, 0, structVarsPrim);
                        structFactPrimVerticalAverage0.FortStructure(master_project_rot_prim,geom_flat);
                    }

                    {
                        MultiFab X, XRot;

                        ComputeVerticalAverage(structFactPrimMF, X, geom, project_dir, 0, structVarsPrim, membrane_cell, n_cells[project_dir]-1);
                        XRot = RotateFlattenedMF(X);
                        master_project_rot_prim.ParallelCopy(XRot, 0, 0, structVarsPrim); 
                        structFactPrimVerticalAverage1.FortStructure(master_project_rot_prim,geom_flat);
                    }

                    {
                        MultiFab X, XRot;

                        ComputeVerticalAverage(structFactConsMF, X, geom, project_dir, 0, structVarsCons, 0, membrane_cell-1);
                        XRot = RotateFlattenedMF(X);
                        master_project_rot_cons.ParallelCopy(XRot, 0, 0, structVarsCons); 
                        structFactConsVerticalAverage0.FortStructure(master_project_rot_cons,geom_flat);
                    }

                    {
                        MultiFab X, XRot;

                        ComputeVerticalAverage(structFactConsMF, X, geom, project_dir, 0, structVarsCons, membrane_cell, n_cells[project_dir]-1);
                        XRot = RotateFlattenedMF(X);
                        master_project_rot_cons.ParallelCopy(XRot, 0, 0, structVarsCons); 
                        structFactConsVerticalAverage1.FortStructure(master_project_rot_cons,geom_flat);
                    }
                }
            }

            if (do_2D) {

                for (int i=0; i<n_cells[2]; ++i) {

                    {
                        MultiFab X, XRot;

                        ExtractSlice(structFactPrimMF, X, geom, 2, i, 0, structVarsPrim);
                        XRot = RotateFlattenedMF(X);
                        master_2D_rot_prim.ParallelCopy(XRot, 0, 0, structVarsPrim); 
                        structFactPrimArray[i].FortStructure(master_2D_rot_prim,geom_flat_2D);
                    }

                    {
                        MultiFab X, XRot;

                        ExtractSlice(structFactConsMF, X, geom, 2, i, 0, structVarsCons);
                        XRot = RotateFlattenedMF(X);
                        master_2D_rot_cons.ParallelCopy(XRot, 0, 0, structVarsCons); 
                        structFactConsArray[i].FortStructure(master_2D_rot_cons,geom_flat_2D);
                    }

                }
            }

            if (n_ads_spec > 0) {
                int surfcov_dir = 2;
                int surfcov_plane = 0;
                int surfcov_structVars = n_ads_spec;
                MultiFab Flattened;  // flattened multifab defined below
                ExtractSlice(surfcov, Flattened, geom, surfcov_dir, surfcov_plane, 0, surfcov_structVars);
                // we rotate this flattened MultiFab to have normal in the z-direction since
                // our structure factor class assumes this for flattened
                MultiFab FlattenedRot = RotateFlattenedMF(Flattened);
                surfcovFlattenedRotMaster.ParallelCopy(FlattenedRot,0,0,surfcov_structVars);
                surfcovStructFact.FortStructure(surfcovFlattenedRotMaster,surfcov_geom_flat);
            }

        }

        // write out structure factor
        if (step > amrex::Math::abs(n_steps_skip) && 
            struct_fact_int > 0 && plot_int > 0 && 
            step%plot_int == 0) {

            if ((do_1D==0) and (do_2D==0)) {
                structFactPrim.WritePlotFile(step,time,geom,"plt_SF_prim");
                structFactCons.WritePlotFile(step,time,geom,"plt_SF_cons");
            }

            if (project_dir >= 0) {
                if (do_slab_sf == 0) {
                    structFactPrimVerticalAverage.WritePlotFile(step,time,geom_flat,"plt_SF_prim_VerticalAverage");
                    structFactConsVerticalAverage.WritePlotFile(step,time,geom_flat,"plt_SF_cons_VerticalAverage");
                }
                else {
                    structFactPrimVerticalAverage0.WritePlotFile(step,time,geom_flat,"plt_SF_prim_VerticalAverageSlab0");
                    structFactPrimVerticalAverage1.WritePlotFile(step,time,geom_flat,"plt_SF_prim_VerticalAverageSlab1");
                    structFactConsVerticalAverage0.WritePlotFile(step,time,geom_flat,"plt_SF_cons_VerticalAverageSlab0");
                    structFactConsVerticalAverage1.WritePlotFile(step,time,geom_flat,"plt_SF_cons_VerticalAverageSlab1");
                }
            }

            if (do_2D) {
                    
                MultiFab prim_mag, prim_realimag, cons_mag, cons_realimag;

                prim_mag.define(ba_flat_2D,dmap_flat_2D,structFactPrimArray[0].get_ncov(),0);
                prim_realimag.define(ba_flat_2D,dmap_flat_2D,2*structFactPrimArray[0].get_ncov(),0);
                cons_mag.define(ba_flat_2D,dmap_flat_2D,structFactConsArray[0].get_ncov(),0);
                cons_realimag.define(ba_flat_2D,dmap_flat_2D,2*structFactConsArray[0].get_ncov(),0);

                prim_mag.setVal(0.0);
                cons_mag.setVal(0.0);
                prim_realimag.setVal(0.0);
                cons_realimag.setVal(0.0);

                for (int i=0; i<n_cells[2]; ++i) {
                    structFactPrimArray[i].AddToExternal(prim_mag,prim_realimag,geom_flat_2D);
                    structFactConsArray[i].AddToExternal(cons_mag,cons_realimag,geom_flat_2D);
                }
                    
                Real ncellsinv = 1.0/n_cells[2];
                prim_mag.mult(ncellsinv);
                cons_mag.mult(ncellsinv);
                prim_realimag.mult(ncellsinv);
                cons_realimag.mult(ncellsinv);

                WritePlotFilesSF_2D(prim_mag,prim_realimag,geom_flat_2D,step,time,
                                    structFactPrimArray[0].get_names(),"plt_SF_prim_2D");
                WritePlotFilesSF_2D(cons_mag,cons_realimag,geom_flat_2D,step,time,
                                    structFactConsArray[0].get_names(),"plt_SF_cons_2D");

            }

            if (n_ads_spec > 0) {
                surfcovStructFact.WritePlotFile(step,time,surfcov_geom_flat,"plt_SF_surfcov");
            }
        }

        // write checkpoint file
        if (chk_int > 0 && step > 0 && step%chk_int == 0)
        {
            if (do_1D) {
                WriteCheckPoint1D(step, time, statsCount, geom, cu, cuMeans, cuVars, prim,
                                  primMeans, primVars, cumom, cumomMeans, cumomVars, 
                                  vel, velMeans, velVars, coVars, spatialCross1D, ncross);
            }
            else if (do_2D) {
                WriteCheckPoint2D(step, time, statsCount, geom, cu, cuMeans, cuVars, prim,
                                  primMeans, primVars, cumom, cumomMeans, cumomVars, 
                                  vel, velMeans, velVars, coVars, spatialCross2D, ncross);
            }
            else {
                WriteCheckPoint3D(step, time, statsCount, geom, cu, cuMeans, cuVars, prim,
                                  primMeans, primVars, cumom, cumomMeans, cumomVars, 
                                  vel, velMeans, velVars, coVars, surfcov, surfcovMeans, surfcovVars, surfcovcoVars, spatialCross3D, ncross, turbforce);
            }
        }

        // timer
        Real aux2 = ParallelDescriptor::second() - aux1;
        ParallelDescriptor::ReduceRealMax(aux2,  ParallelDescriptor::IOProcessorNumber());
        if (step%100 == 0) {
            amrex::Print() << "Aux time (stats, struct fac, plotfiles) " << aux2 << " seconds\n";
        }
        
        time = time + dt;

        // MultiFab memory usage
        const int IOProc = ParallelDescriptor::IOProcessorNumber();

        amrex::Long min_fab_megabytes  = amrex::TotalBytesAllocatedInFabsHWM()/1048576;
        amrex::Long max_fab_megabytes  = min_fab_megabytes;

        ParallelDescriptor::ReduceLongMin(min_fab_megabytes, IOProc);
        ParallelDescriptor::ReduceLongMax(max_fab_megabytes, IOProc);

        if (step%100 == 0) {
            amrex::Print() << "High-water FAB megabyte spread across MPI nodes: ["
                           << min_fab_megabytes << " ... " << max_fab_megabytes << "]\n";
        }

        min_fab_megabytes  = amrex::TotalBytesAllocatedInFabs()/1048576;
        max_fab_megabytes  = min_fab_megabytes;

        ParallelDescriptor::ReduceLongMin(min_fab_megabytes, IOProc);
        ParallelDescriptor::ReduceLongMax(max_fab_megabytes, IOProc);

        if (step%100 == 0) {
            amrex::Print() << "Curent     FAB megabyte spread across MPI nodes: ["
                           << min_fab_megabytes << " ... " << max_fab_megabytes << "]\n";
        }

        if (step%100 == 0) {
            amrex::Real cfl_max = GetMaxAcousticCFL(prim, vel, dt, geom);
            amrex::Print() << "Max convective-acoustic CFL is: " << cfl_max << "\n";
        }
    }

    if (ParallelDescriptor::IOProcessor()) outfile.close();
#if defined(TURB)
    if (turbForcing >= 1) {
        if (ParallelDescriptor::IOProcessor()) turboutfile.close();
        if (ParallelDescriptor::IOProcessor()) turboutfiledecomp.close();
    }
#endif

    // timer
    Real stop_time = ParallelDescriptor::second() - strt_time;
    ParallelDescriptor::ReduceRealMax(stop_time, ParallelDescriptor::IOProcessorNumber());
    amrex::Print() << "Run time = " << stop_time << std::endl;
}
