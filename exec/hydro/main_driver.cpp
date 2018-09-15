
#include "hydro_test_functions.H"
#include "hydro_test_functions_F.H"

#include "hydro_functions.H"
#include "hydro_functions_F.H"
#include "StochMFlux.H"
#include "StructFact.H"

#include "rng_functions_F.H"

#include "common_functions.H"
#include "common_functions_F.H"

#include "gmres_functions.H"
#include "gmres_functions_F.H"

#include "common_namespace.H"
#include "common_namespace_declarations.H"

#include "gmres_namespace.H"
#include "gmres_namespace_declarations.H"

#include <AMReX_VisMF.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_MultiFabUtil.H>

using namespace amrex;
using namespace common;
using namespace gmres;

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
    read_gmres_namelist(inputs_file.c_str(),inputs_file.size()+1);

    // copy contents of F90 modules to C++ namespaces
    InitializeCommonNamespace();
    InitializeGmresNamespace();

    // is the problem periodic?
    Vector<int> is_periodic(AMREX_SPACEDIM,0);  // set to 0 (not periodic) by default
    for (int i=0; i<AMREX_SPACEDIM; ++i) {
        if (bc_lo[i] == -1 && bc_hi[i] == -1) {
            is_periodic[i] = 1;
        }
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
    Real dtinv = 1.0/dt;
    const Real* dx = geom.CellSize();

    // how boxes are distrubuted among MPI processes
    DistributionMapping dmap(ba);

    /////////////////////////////////////////
    //Initialise rngs
    /////////////////////////////////////////
    const int n_rngs = 1;

    int fhdSeed = 1;
    int particleSeed = 2;
    int selectorSeed = 3;
    int thetaSeed = 4;
    int phiSeed = 5;
    int generalSeed = 6;

    // const int proc = ParallelDescriptor::MyProc();
    // fhdSeed += proc;
    // particleSeed += proc;
    // selectorSeed += proc;
    // thetaSeed += proc;
    // phiSeed += proc;
    // generalSeed += proc;

    //Initialise rngs
    rng_initialize(&fhdSeed,&particleSeed,&selectorSeed,&thetaSeed,&phiSeed,&generalSeed);
    /////////////////////////////////////////

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

    // cell-centered gamma
    MultiFab gamma(ba, dmap, 1, 1);
    gamma.setVal(0.);

    ///////////////////////////////////////////
    // Scaled alpha, beta, gamma:
    ///////////////////////////////////////////
    // alpha_fc_0 arrays
    std::array< MultiFab, AMREX_SPACEDIM > alpha_fc_0;
    AMREX_D_TERM(alpha_fc_0[0].define(convert(ba,nodal_flag_x), dmap, 1, 1);,
                 alpha_fc_0[1].define(convert(ba,nodal_flag_y), dmap, 1, 1);,
                 alpha_fc_0[2].define(convert(ba,nodal_flag_z), dmap, 1, 1););
    AMREX_D_TERM(alpha_fc_0[0].setVal(0.);,
                 alpha_fc_0[1].setVal(0.);,
                 alpha_fc_0[2].setVal(0.););

    // Scaled by 1/2:
    // beta_hlf cell centered
    MultiFab beta_hlf(ba, dmap, 1, 1);
    MultiFab::Copy(beta_hlf, beta, 0, 0, 1, 1);
    beta_hlf.mult(0.5, 1);

    // beta_hlf on nodes in 2d
    // beta_hlf on edges in 3d
    std::array< MultiFab, NUM_EDGE > beta_ed_hlf;
#if (AMREX_SPACEDIM == 2)
    beta_ed_hlf[0].define(convert(ba,nodal_flag), dmap, 1, 1);
    MultiFab::Copy(beta_ed_hlf[0], beta_ed[0], 0, 0, 1, 1);
    beta_ed_hlf[0].mult(0.5, 1);
#elif (AMREX_SPACEDIM == 3)
    beta_ed_hlf[0].define(convert(ba,nodal_flag_xy), dmap, 1, 1);
    beta_ed_hlf[1].define(convert(ba,nodal_flag_xz), dmap, 1, 1);
    beta_ed_hlf[2].define(convert(ba,nodal_flag_yz), dmap, 1, 1);
    MultiFab::Copy(beta_ed_hlf[0], beta_ed[0], 0, 0, 1, 1);
    MultiFab::Copy(beta_ed_hlf[1], beta_ed[1], 0, 0, 1, 1);
    MultiFab::Copy(beta_ed_hlf[2], beta_ed[2], 0, 0, 1, 1);
    beta_ed_hlf[0].mult(0.5, 1);
    beta_ed_hlf[1].mult(0.5, 1);
    beta_ed_hlf[2].mult(0.5, 1);
#endif

    // cell-centered gamma_hlf
    MultiFab gamma_hlf(ba, dmap, 1, 1);
    MultiFab::Copy(gamma_hlf, gamma, 0, 0, 1, 1);
    gamma_hlf.mult(-0.5, 1);

    // Scaled by -1/2:
    // beta_neghlf cell centered
    MultiFab beta_neghlf(ba, dmap, 1, 1);
    MultiFab::Copy(beta_neghlf, beta, 0, 0, 1, 1);
    beta_neghlf.mult(-0.5, 1);

    // beta_neghlf on nodes in 2d
    // beta_neghlf on edges in 3d
    std::array< MultiFab, NUM_EDGE > beta_ed_neghlf;
#if (AMREX_SPACEDIM == 2)
    beta_ed_neghlf[0].define(convert(ba,nodal_flag), dmap, 1, 1);
    MultiFab::Copy(beta_ed_neghlf[0], beta_ed[0], 0, 0, 1, 1);
    beta_ed_neghlf[0].mult(-0.5, 1);
#elif (AMREX_SPACEDIM == 3)
    beta_ed_neghlf[0].define(convert(ba,nodal_flag_xy), dmap, 1, 1);
    beta_ed_neghlf[1].define(convert(ba,nodal_flag_xz), dmap, 1, 1);
    beta_ed_neghlf[2].define(convert(ba,nodal_flag_yz), dmap, 1, 1);
    MultiFab::Copy(beta_ed_neghlf[0], beta_ed[0], 0, 0, 1, 1);
    MultiFab::Copy(beta_ed_neghlf[1], beta_ed[1], 0, 0, 1, 1);
    MultiFab::Copy(beta_ed_neghlf[2], beta_ed[2], 0, 0, 1, 1);
    beta_ed_neghlf[0].mult(-0.5, 1);
    beta_ed_neghlf[1].mult(-0.5, 1);
    beta_ed_neghlf[2].mult(-0.5, 1);
#endif

    // cell-centered gamma
    MultiFab gamma_neghlf(ba, dmap, 1, 1);
    MultiFab::Copy(gamma_neghlf, gamma, 0, 0, 1, 1);
    gamma_neghlf.mult(-0.5, 1);
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

    // FIXME: eta & temperature = 1.0 (uniform) here for convienience

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
    std::array< MultiFab, AMREX_SPACEDIM >  mfluxdiv;
    // Define mfluxdiv multifabs
    mfluxdiv[0].define(convert(ba,nodal_flag_x), dmap, 1, 1);
    mfluxdiv[1].define(convert(ba,nodal_flag_y), dmap, 1, 1);
#if (AMREX_SPACEDIM == 3)
    mfluxdiv[2].define(convert(ba,nodal_flag_z), dmap, 1, 1);
#endif

    Vector< amrex::Real > weights;
    // weights = {std::sqrt(0.5), std::sqrt(0.5)};
    weights = {1.0};
    
    // Declare object of StochMFlux class 
    // StochMFlux sMflux (ba,dmap,geom);
    StochMFlux sMflux (ba,dmap,geom,n_rngs);

    ///////////////////////////////////////////

    /////////////// Test/Hack /////////////////////////
    // sMflux.fillMStochastic();
    // sMflux.stochMforce(mfluxdiv,eta_cc,eta_ed,temp_cc,temp_ed,weights,dt);
    // sMflux.writeMFs(mfluxdiv);

    // Abort("Done with hack");
    // exit(0);
    //////////////////////////////////////////////////

    // tracer
    MultiFab tracer(ba,dmap,1,1);
    MultiFab tracerPred(ba,dmap,1,1);
    MultiFab advFluxdivS(ba,dmap,1,1);

    // rhs_p GMRES solve
    MultiFab gmres_rhs_p(ba, dmap, 1, 0);
    gmres_rhs_p.setVal(0.);

    // pressure for GMRES solve
    MultiFab pres(ba,dmap,1,1);
    pres.setVal(0.);  // initial guess

    // rhs_u GMRES solve
    std::array< MultiFab, AMREX_SPACEDIM > gmres_rhs_u;
    AMREX_D_TERM(gmres_rhs_u[0].define(convert(ba,nodal_flag_x), dmap, 1, 1);,
                 gmres_rhs_u[1].define(convert(ba,nodal_flag_y), dmap, 1, 1);,
                 gmres_rhs_u[2].define(convert(ba,nodal_flag_z), dmap, 1, 1););

    // staggered velocities
    std::array< MultiFab, AMREX_SPACEDIM > umac;
    AMREX_D_TERM(umac[0].define(convert(ba,nodal_flag_x), dmap, 1, 1);,
                 umac[1].define(convert(ba,nodal_flag_y), dmap, 1, 1);,
                 umac[2].define(convert(ba,nodal_flag_z), dmap, 1, 1););

    std::array< MultiFab, AMREX_SPACEDIM > umacNew;
    AMREX_D_TERM(umacNew[0].define(convert(ba,nodal_flag_x), dmap, 1, 1);,
                 umacNew[1].define(convert(ba,nodal_flag_y), dmap, 1, 1);,
                 umacNew[2].define(convert(ba,nodal_flag_z), dmap, 1, 1););

    std::array< MultiFab, AMREX_SPACEDIM > Lumac;
    AMREX_D_TERM(Lumac[0].define(convert(ba,nodal_flag_x), dmap, 1, 1);,
                 Lumac[1].define(convert(ba,nodal_flag_y), dmap, 1, 1);,
                 Lumac[2].define(convert(ba,nodal_flag_z), dmap, 1, 1););

    // advective terms
    std::array< MultiFab, AMREX_SPACEDIM > advFluxdiv;
    AMREX_D_TERM(advFluxdiv[0].define(convert(ba,nodal_flag_x), dmap, 1, 1);,
                 advFluxdiv[1].define(convert(ba,nodal_flag_y), dmap, 1, 1);,
                 advFluxdiv[2].define(convert(ba,nodal_flag_z), dmap, 1, 1););

    std::array< MultiFab, AMREX_SPACEDIM > advFluxdivPred;
    AMREX_D_TERM(advFluxdivPred[0].define(convert(ba,nodal_flag_x), dmap, 1, 1);,
                 advFluxdivPred[1].define(convert(ba,nodal_flag_y), dmap, 1, 1);,
                 advFluxdivPred[2].define(convert(ba,nodal_flag_z), dmap, 1, 1););

    // staggered momentum
    std::array< MultiFab, AMREX_SPACEDIM > uMom;
    AMREX_D_TERM(uMom[0].define(convert(ba,nodal_flag_x), dmap, 1, 1);,
                 uMom[1].define(convert(ba,nodal_flag_y), dmap, 1, 1);,
                 uMom[2].define(convert(ba,nodal_flag_z), dmap, 1, 1););

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

    Vector< MultiFab > struct_in_cc;
    struct_in_cc.resize(AMREX_SPACEDIM);
    for (int d=0; d<struct_in_cc.size(); d++) {
      struct_in_cc[d].define(ba, dmap, 1, 0);
    }
    
    StructFact structFact(ba,dmap,var_names);

    ///////////////////////////////////////////

    const RealBox& realDomain = geom.ProbDomain();
    int dm;

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

    	// initialize tracer
        init_s_vel(BL_TO_FORTRAN_BOX(bx),
    		   BL_TO_FORTRAN_ANYD(tracer[mfi]),
    		   dx, ZFILL(realDomain.lo()), ZFILL(realDomain.hi()));

    }

    // initial guess for new solution
    AMREX_D_TERM(MultiFab::Copy(umacNew[0], umac[0], 0, 0, 1, 0);,
                 MultiFab::Copy(umacNew[1], umac[1], 0, 0, 1, 0);,
                 MultiFab::Copy(umacNew[2], umac[2], 0, 0, 1, 0););

    int step = 0;
    Real time = 0.;

    Real norm_pre_rhs;

    // write out initial state
    if (plot_int > 0) 
      {
	WritePlotFile(step,time,geom,umac,tracer,pres);
      }
    
    //////////////////////////
    //// FFT test
    if (struct_fact_int > 0) {
      // // std::array <MultiFab, AMREX_SPACEDIM> mf_cc;
      // // mf_cc[0].define(ba, dmap, 1, 0);
      // // mf_cc[1].define(ba, dmap, 1, 0);
      // // mf_cc[2].define(ba, dmap, 1, 0);
      // // for ( MFIter mfi(beta); mfi.isValid(); ++mfi ) {
      // //   const Box& bx = mfi.validbox();
      // //   init_s_vel(BL_TO_FORTRAN_BOX(bx),
      // // 		   BL_TO_FORTRAN_ANYD(mf_cc[0][mfi]),
      // // 		   dx, ZFILL(realDomain.lo()), ZFILL(realDomain.hi()));
      // // }

      // // structFact.ComputeFFT(mf_cc,geom);
      // structFact.FortStructure(umac,geom);
      // structFact.WritePlotFile(step,time,geom,1.0);
      // exit(0);
    }
    //////////////////////////

    //Time stepping loop
    for(step=1;step<=max_step;++step) {

        AMREX_D_TERM(umac[0].FillBoundary(geom.periodicity());,
                     umac[1].FillBoundary(geom.periodicity());,
                     umac[2].FillBoundary(geom.periodicity()););

	// Fill stochastic terms
	sMflux.fillMStochastic();
	
	//////////////////////////
	// Advance tracer
	//////////////////////////

    	// Compute tracer:
	tracer.FillBoundary(geom.periodicity());
	MkAdvSFluxdiv(umac,tracer,advFluxdivS,dx,geom,0);
	advFluxdivS.mult(dt, 1);

	// compute predictor
	MultiFab::Copy(tracerPred, tracer, 0, 0, 1, 0);
	MultiFab::Add(tracerPred, advFluxdivS, 0, 0, 1, 0);
	tracerPred.FillBoundary(geom.periodicity());
	MkAdvSFluxdiv(umac,tracerPred,advFluxdivS,dx,geom,0);
	advFluxdivS.mult(dt, 1);

	// advance in time
	MultiFab::Add(tracer, tracerPred, 0, 0, 1, 0);
	MultiFab::Add(tracer, advFluxdivS, 0, 0, 1, 0);
	tracer.mult(0.5, 1);

	// amrex::Print() << "tracer L0 norm = " << tracer.norm0() << "\n";
	//////////////////////////

	//////////////////////////////////////////////////
	// ADVANCE
	//////////////////////////////////////////////////

    	// PREDICTOR STEP (heun's method: part 1)
    	// compute advective term
        AMREX_D_TERM(MultiFab::Copy(uMom[0], umac[0], 0, 0, 1, 0);,
                     MultiFab::Copy(uMom[1], umac[1], 0, 0, 1, 0);,
                     MultiFab::Copy(uMom[2], umac[2], 0, 0, 1, 0););

    	// let rho = 1
    	for (int d=0; d<AMREX_SPACEDIM; d++) {
    	  uMom[d].mult(1.0, 1);
    	}

        AMREX_D_TERM(uMom[0].FillBoundary(geom.periodicity());,
                     uMom[1].FillBoundary(geom.periodicity());,
                     uMom[2].FillBoundary(geom.periodicity()););

    	MkAdvMFluxdiv(umac,uMom,advFluxdiv,dx,0);

        // crank-nicolson terms
        StagApplyOp(beta_neghlf,gamma_neghlf,beta_ed_neghlf,umac,Lumac,alpha_fc_0,dx,theta_alpha);

	// compute stochastic force terms
	sMflux.stochMforce(mfluxdiv,eta_cc,eta_ed,temp_cc,temp_ed,weights,dt);

    	AMREX_D_TERM(MultiFab::Copy(gmres_rhs_u[0], umac[0], 0, 0, 1, 0);,
                     MultiFab::Copy(gmres_rhs_u[1], umac[1], 0, 0, 1, 0);,
                     MultiFab::Copy(gmres_rhs_u[2], umac[2], 0, 0, 1, 0););
    	for (int d=0; d<AMREX_SPACEDIM; d++) {
    	  gmres_rhs_u[d].mult(dtinv, 1);
    	}
    	AMREX_D_TERM(MultiFab::Add(gmres_rhs_u[0], mfluxdiv[0], 0, 0, 1, 0);,
                     MultiFab::Add(gmres_rhs_u[1], mfluxdiv[1], 0, 0, 1, 0);,
                     MultiFab::Add(gmres_rhs_u[2], mfluxdiv[2], 0, 0, 1, 0););
	AMREX_D_TERM(MultiFab::Add(gmres_rhs_u[0], Lumac[0], 0, 0, 1, 0);,
                     MultiFab::Add(gmres_rhs_u[1], Lumac[1], 0, 0, 1, 0);,
                     MultiFab::Add(gmres_rhs_u[2], Lumac[2], 0, 0, 1, 0););
    	AMREX_D_TERM(MultiFab::Add(gmres_rhs_u[0], advFluxdiv[0], 0, 0, 1, 0);,
                     MultiFab::Add(gmres_rhs_u[1], advFluxdiv[1], 0, 0, 1, 0);,
                     MultiFab::Add(gmres_rhs_u[2], advFluxdiv[2], 0, 0, 1, 0););

        AMREX_D_TERM(gmres_rhs_u[0].FillBoundary(geom.periodicity());,
                     gmres_rhs_u[1].FillBoundary(geom.periodicity());,
                     gmres_rhs_u[2].FillBoundary(geom.periodicity()););

    	// initial guess for new solution
    	AMREX_D_TERM(MultiFab::Copy(umacNew[0], umac[0], 0, 0, 1, 0);,
    		     MultiFab::Copy(umacNew[1], umac[1], 0, 0, 1, 0);,
    		     MultiFab::Copy(umacNew[2], umac[2], 0, 0, 1, 0););
    	pres.setVal(0.);  // initial guess

        // call GMRES to compute predictor
    	GMRES(gmres_rhs_u,gmres_rhs_p,umacNew,pres,alpha_fc,beta_hlf,beta_ed_hlf,gamma_hlf,theta_alpha,geom,norm_pre_rhs);

    	// Compute predictor advective term
        AMREX_D_TERM(umacNew[0].FillBoundary(geom.periodicity());,
                     umacNew[1].FillBoundary(geom.periodicity());,
                     umacNew[2].FillBoundary(geom.periodicity()););

        AMREX_D_TERM(MultiFab::Copy(uMom[0], umacNew[0], 0, 0, 1, 0);,
                     MultiFab::Copy(uMom[1], umacNew[1], 0, 0, 1, 0);,
                     MultiFab::Copy(uMom[2], umacNew[2], 0, 0, 1, 0););

    	// let rho = 1
    	for (int d=0; d<AMREX_SPACEDIM; d++) {
    	  uMom[d].mult(1.0, 1);
    	}

        AMREX_D_TERM(uMom[0].FillBoundary(geom.periodicity());,
                     uMom[1].FillBoundary(geom.periodicity());,
                     uMom[2].FillBoundary(geom.periodicity()););

    	MkAdvMFluxdiv(umacNew,uMom,advFluxdivPred,dx,0);

    	// ADVANCE STEP (crank-nicolson + heun's method)

    	// Compute gmres_rhs

        // trapezoidal advective terms
    	for (int d=0; d<AMREX_SPACEDIM; d++) {
    	  advFluxdiv[d].mult(0.5, 1);
    	  advFluxdivPred[d].mult(0.5, 1);
    	}

        // crank-nicolson terms
        StagApplyOp(beta_neghlf,gamma_neghlf,beta_ed_neghlf,umac,Lumac,alpha_fc_0,dx,theta_alpha);

	// compute stochastic force terms
	sMflux.stochMforce(mfluxdiv,eta_cc,eta_ed,temp_cc,temp_ed,weights,dt);

    	AMREX_D_TERM(MultiFab::Copy(gmres_rhs_u[0], umac[0], 0, 0, 1, 0);,
                     MultiFab::Copy(gmres_rhs_u[1], umac[1], 0, 0, 1, 0);,
                     MultiFab::Copy(gmres_rhs_u[2], umac[2], 0, 0, 1, 0););
    	for (int d=0; d<AMREX_SPACEDIM; d++) {
    	  gmres_rhs_u[d].mult(dtinv, 1);
    	}
    	AMREX_D_TERM(MultiFab::Add(gmres_rhs_u[0], Lumac[0], 0, 0, 1, 0);,
                     MultiFab::Add(gmres_rhs_u[1], Lumac[1], 0, 0, 1, 0);,
                     MultiFab::Add(gmres_rhs_u[2], Lumac[2], 0, 0, 1, 0););
	AMREX_D_TERM(MultiFab::Add(gmres_rhs_u[0], mfluxdiv[0], 0, 0, 1, 0);,
                     MultiFab::Add(gmres_rhs_u[1], mfluxdiv[1], 0, 0, 1, 0);,
                     MultiFab::Add(gmres_rhs_u[2], mfluxdiv[2], 0, 0, 1, 0););
    	AMREX_D_TERM(MultiFab::Add(gmres_rhs_u[0], advFluxdiv[0], 0, 0, 1, 0);,
                     MultiFab::Add(gmres_rhs_u[1], advFluxdiv[1], 0, 0, 1, 0);,
                     MultiFab::Add(gmres_rhs_u[2], advFluxdiv[2], 0, 0, 1, 0););
    	AMREX_D_TERM(MultiFab::Add(gmres_rhs_u[0], advFluxdivPred[0], 0, 0, 1, 0);,
                     MultiFab::Add(gmres_rhs_u[1], advFluxdivPred[1], 0, 0, 1, 0);,
                     MultiFab::Add(gmres_rhs_u[2], advFluxdivPred[2], 0, 0, 1, 0););

        AMREX_D_TERM(gmres_rhs_u[0].FillBoundary(geom.periodicity());,
                     gmres_rhs_u[1].FillBoundary(geom.periodicity());,
                     gmres_rhs_u[2].FillBoundary(geom.periodicity()););

    	// initial guess for new solution
    	AMREX_D_TERM(MultiFab::Copy(umacNew[0], umac[0], 0, 0, 1, 0);,
    		     MultiFab::Copy(umacNew[1], umac[1], 0, 0, 1, 0);,
    		     MultiFab::Copy(umacNew[2], umac[2], 0, 0, 1, 0););
    	pres.setVal(0.);  // initial guess

        // call GMRES here
    	GMRES(gmres_rhs_u,gmres_rhs_p,umacNew,pres,alpha_fc,beta_hlf,beta_ed_hlf,gamma_hlf,theta_alpha,geom,norm_pre_rhs);

        AMREX_D_TERM(MultiFab::Copy(umac[0], umacNew[0], 0, 0, 1, 0);,
                     MultiFab::Copy(umac[1], umacNew[1], 0, 0, 1, 0);,
                     MultiFab::Copy(umac[2], umacNew[2], 0, 0, 1, 0););
	//////////////////////////////////////////////////
	
	///////////////////////////////////////////
	// Update structure factor
	///////////////////////////////////////////
	if (step > n_steps_skip && struct_fact_int > 0 && (step-n_steps_skip-1)%struct_fact_int == 0) {
	  for(int d=0; d<AMREX_SPACEDIM; d++) {
	    ShiftFaceToCC(umac[d], 0, struct_in_cc[d], 0, 1);
	  }
	  structFact.FortStructure(struct_in_cc,geom);
	  // Print() << "Executed FortStructure() at step = " << step << std::endl;
        }
	///////////////////////////////////////////

        amrex::Print() << "Advanced step " << step << "\n";

        time = time + dt;

        if (plot_int > 0 && step%plot_int == 0) {
          // write out umac & pres to a plotfile
    	  WritePlotFile(step,time,geom,umac,tracer,pres);
        }
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
      Real SFscale = dVol/(k_B*temp_const);
      // Print() << "Hack: structure factor scaling = " << SFscale << std::endl;

      structFact.WritePlotFile(step,time,geom,SFscale);
      // amrex::Vector< MultiFab > struct_out;
      // structFact.StructOut(struct_out);
    }

    // Call the timer again and compute the maximum difference between the start time 
    // and stop time over all processors
    Real stop_time = ParallelDescriptor::second() - strt_time;
    ParallelDescriptor::ReduceRealMax(stop_time);
    amrex::Print() << "Run time = " << stop_time << std::endl;
    amrex::Print() << "dt = " << dt << "\n";

}
