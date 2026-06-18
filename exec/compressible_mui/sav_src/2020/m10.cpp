#include "common_functions.H"
#include "compressible_functions.H"


#include "rng_functions.H"

#include "StructFact.H"

using namespace amrex;

// mui
#include <mui.h>
using namespace mui;

// this routine pushes the following information to MUI
// - species number densities and temperature of FHD cells contacting the interface
void mui_push(MultiFab& cu, MultiFab& prim, const amrex::Real* dx, mui::uniface2d &uniface, const int step)
{
    // assuming the interface is perpendicular to the z-axis
    // and includes cells with the smallest value of z (i.e. k=0)

    for (MFIter mfi(cu,false); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        Dim3 lo = lbound(bx);
        Dim3 hi = ubound(bx);
        const Array4<Real> & cu_fab = cu.array(mfi);
        const Array4<Real> & prim_fab = prim.array(mfi);

        // unless bx contains cells at the interface, skip
        int k = 0;
        if (k<lo.z || k>hi.z) continue;

        for (int j = lo.y; j<= hi.y; ++j) {
            for (int i = lo.x; i<=hi.x; ++i) {

                double x = prob_lo[0]+(i+0.5)*dx[0];
                double y = prob_lo[1]+(j+0.5)*dx[1];

                std::string channel;

                for (int n = 0; n < nspecies; ++n) {

                    channel = "CH_density";
                    channel += '0'+(n+1);   // assuming nspecies<10

                    double dens = cu_fab(i,j,k,5+n);    // mass density
                    dens *= 6.02e23/molmass[n];         // number density

                    uniface.push(channel,{x,y},dens);
                }

                channel = "CH_temp";

                uniface.push(channel,{x,y},prim_fab(i,j,k,4));
            }
        }
    }

    uniface.commit(step);

    return;
}

// this routine fetches the following information from MUI:
// - adsoprtion and desoprtion counts of each species between time points
void mui_fetch(MultiFab& cu, MultiFab& prim, const amrex::Real* dx, mui::uniface2d &uniface, const int step)
{
    // assuming the interface is perpendicular to the z-axis
    // and includes cells with the smallest value of z (i.e. k=0)

    mui::sampler_kmc_fhd2d<int> s({dx[0],dx[1]});
    mui::chrono_sampler_exact2d t;

    for (MFIter mfi(cu,false); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        Dim3 lo = lbound(bx);
        Dim3 hi = ubound(bx);
        const Array4<Real> & cu_fab = cu.array(mfi);
        const Array4<Real> & prim_fab = prim.array(mfi);

        // unless bx contains cells at the interface, skip
        int k = 0;
        if (k<lo.z || k>hi.z) continue;

        for (int j = lo.y; j<= hi.y; ++j) {
            for (int i = lo.x; i<=hi.x; ++i) {

                double x = prob_lo[0]+(i+0.5)*dx[0];
                double y = prob_lo[1]+(j+0.5)*dx[1];
                double dV = dx[0]*dx[1]*dx[2];
                double temp = prim_fab(i,j,k,4);

                for (int n = 0; n < nspecies; ++n) {

                    std::string channel;
                    int ac,dc;

                    channel = "CH_ac";
                    channel += '0'+(n+1);   // assuming nspecies<10
                    ac = uniface.fetch(channel,{x,y},step,s,t);

                    channel = "CH_dc";
                    channel += '0'+(n+1);   // assuming nspecies<10
                    dc = uniface.fetch(channel,{x,y},step,s,t);

                    double mass = molmass[n]/6.02e23;
                    double kBTm = k_B*temp/mass;
                    double sqrtkBTm = sqrt(kBTm);
                    double vx,vy,vz;
                    double dmomx,dmomy,dmomz,derg;

                    dmomx = dmomy = dmomz = derg = 0.;

                    for (int l=0;l<ac;l++)
                    {
                        // colliding velocity
                        vx = RandomNormal(0.,sqrtkBTm);
                        vy = RandomNormal(0.,sqrtkBTm);
                        vz = -sqrt(-2.*kBTm*log(1.-Random()));

                        dmomx -= mass*vx;
                        dmomy -= mass*vy;
                        dmomz -= mass*vz;
                        derg  -= 0.5*mass*(vx*vx+vy*vy+vz*vz);
                    }

                    for (int l=0;l<dc;l++)
                    {
                        // new velocity
                        vx = RandomNormal(0.,sqrtkBTm);
                        vy = RandomNormal(0.,sqrtkBTm);
                        vz = sqrt(-2.*kBTm*log(1.-Random()));

                        dmomx += mass*vx;
                        dmomy += mass*vy;
                        dmomz += mass*vz;
                        derg  += 0.5*mass*(vx*vx+vy*vy+vz*vz);
                    }

                    cu_fab(i,j,k,0) += (dc-ac)*mass/dV;
                    cu_fab(i,j,k,5+n) += (dc-ac)*mass/dV;

                    cu_fab(i,j,k,1) += dmomx/dV;
                    cu_fab(i,j,k,2) += dmomy/dV;
                    cu_fab(i,j,k,3) += dmomz/dV;
                    cu_fab(i,j,k,4) += derg/dV;
                }
            }
        }
    }

    uniface.forget(step);

    return;
}

// argv contains the name of the inputs file entered at the command line
void main_driver(const char* argv)
{
    BL_PROFILE_VAR("main_driver()",main_driver);

    // store the current time so we can later compute total run time.
    Real strt_time = ParallelDescriptor::second();

    std::string inputs_file = argv;

    read_compressible_namelist(inputs_file.c_str(),inputs_file.size()+1);

    // copy contents of F90 modules to C++ namespaces
    InitializeCommonNamespace();
    InitializeCompressibleNamespace();

    // if gas heat capacities in the namelist are negative, calculate them using using dofs.
    // This will only update the Fortran values.
    get_hc_gas();
    // now update C++ values
    for (int i=0; i<nspecies; ++i) {
        if (hcv[i] < 0.) {
            hcv[i] = 0.5*dof[i]*Runiv/molmass[i];
        }
        if (hcp[i] < 0.) {
            hcp[i] = 0.5*(2.+dof[i])*Runiv/molmass[i];
        }
    }

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
    int fhdSeed      = seed;
    int particleSeed = seed;
    int selectorSeed = seed;
    int thetaSeed    = seed;
    int phiSeed      = seed;
    int generalSeed  = seed;

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

    // initializes the seed for C++ random number calls
    InitRandom(seed+ParallelDescriptor::MyProc(),
               ParallelDescriptor::NProcs(),
               seed+ParallelDescriptor::MyProc());

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

    MultiFab cu_temp(ba,dmap,nvars,ngc);

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
    MultiFab cuMeansAv(ba,dmap,nvars,ngc);
    MultiFab cuVarsAv (ba,dmap,nvars,ngc);

    cuMeans.setVal(0.0);
    cuVars.setVal(0.0);

    MultiFab primVertAvg;  // flattened multifab defined below

    MultiFab primMeans  (ba,dmap,nprimvars  ,ngc);
    MultiFab primVars   (ba,dmap,nprimvars+5,ngc);
    MultiFab primMeansAv(ba,dmap,nprimvars  ,ngc);
    MultiFab primVarsAv (ba,dmap,nprimvars+5,ngc);
    primMeans.setVal(0.0);
    primVars.setVal(0.0);

    //Miscstats
    // 0        time averaged kinetic energy density

    MultiFab miscStats(ba,dmap,10,ngc);
    Real miscVals[20];
    MultiFab spatialCross(ba,dmap,6,ngc);
    MultiFab spatialCrossAv(ba,dmap,6,ngc);

    miscStats.setVal(0.0);
    spatialCross.setVal(0.0);
    spatialCrossAv.setVal(0.0);

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

    // "primitive" variable structure factor will contain
    // rho
    // vel
    // T
    // Yk
    int structVarsPrim = AMREX_SPACEDIM+nspecies+2;

    Vector< std::string > prim_var_names;
    prim_var_names.resize(structVarsPrim);

    int cnt = 0;
    std::string x;

    // rho
    prim_var_names[cnt++] = "rho";

    // velx, vely, velz
    for (int d=0; d<AMREX_SPACEDIM; d++) {
      x = "vel";
      x += (120+d);
      prim_var_names[cnt++] = x;
    }

    // Temp
    prim_var_names[cnt++] = "Temp";

    // Yk
    for (int d=0; d<nspecies; d++) {
      x = "Y";
      x += (49+d);
      prim_var_names[cnt++] = x;
    }

    MultiFab structFactPrimMF;
    structFactPrimMF.define(ba, dmap, structVarsPrim, 0);

    // scale SF results by inverse cell volume
    Vector<Real> var_scaling;
    var_scaling.resize(structVarsPrim*(structVarsPrim+1)/2);
    for (int d=0; d<var_scaling.size(); ++d) {
        var_scaling[d] = 1./(dx[0]*dx[1]*dx[2]);
    }

    // compute all pairs
    // note: StructFactPrim option to compute only speicified pairs not written yet
    StructFact structFactPrim(ba,dmap,prim_var_names,var_scaling);

    //////////////////////////////////////////////

    // "conserved" variable structure factor will contain
    // rho
    // j
    // rho*E
    // rho*Yk
    // Temperature (not in the conserved array; will have to copy it in)
    int structVarsCons = AMREX_SPACEDIM+nspecies+3;

    Vector< std::string > cons_var_names;
    cons_var_names.resize(structVarsCons);

    cnt = 0;

    // rho
    cons_var_names[cnt++] = "rho";

    // velx, vely, velz
    for (int d=0; d<AMREX_SPACEDIM; d++) {
      x = "j";
      x += (120+d);
      cons_var_names[cnt++] = x;
    }

    // rho*E
    cons_var_names[cnt++] = "rhoE";

    // rho*Yk
    for (int d=0; d<nspecies; d++) {
      x = "rhoY";
      x += (49+d);
      cons_var_names[cnt++] = x;
    }

    // Temp
    cons_var_names[cnt++] = "Temp";

    MultiFab structFactConsMF;
    structFactConsMF.define(ba, dmap, structVarsCons, 0);

    // scale SF results by inverse cell volume
    var_scaling.resize(structVarsCons*(structVarsCons+1)/2);
    for (int d=0; d<var_scaling.size(); ++d) {
        var_scaling[d] = 1./(dx[0]*dx[1]*dx[2]);
    }

    // compute all pairs
    // note: StructFactCons option to compute only speicified pairs not written yet
    StructFact structFactCons(ba,dmap,cons_var_names,var_scaling);

    //////////////////////////////////////////////

    // structure factor class for vertically-averaged dataset
    StructFact structFactPrimVerticalAverage;

    Geometry geom_flat;

    if(project_dir >= 0){
      prim.setVal(0.0);
      ComputeVerticalAverage(prim, primVertAvg, project_dir, 0, structVarsPrim);
      BoxArray ba_flat = primVertAvg.boxArray();
      const DistributionMapping& dmap_flat = primVertAvg.DistributionMap();
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

      structFactPrimVerticalAverage.~StructFact(); // destruct
      new(&structFactPrimVerticalAverage) StructFact(ba_flat,dmap_flat,prim_var_names,var_scaling); // reconstruct

    }

    ///////////////////////////////////////////

    // Initialize everything

    prim.setVal(0.0,0,nprimvars,ngc);
    prim.setVal(rho0,0,1,ngc);      // density
    prim.setVal(0.,1,3,ngc);        // x/y/z velocity
    prim.setVal(T_init[0],4,1,ngc); // temperature
                                    // pressure computed later in conservedToPrimitive
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
        InitConsVar(cu,prim,geom);
    }

    statsCount = 1;

    // Write initial plotfile
    conservedToPrimitive(prim, cu);

    // Set BC: 1) fill boundary 2) physical
    cu.FillBoundary(geom.periodicity());
    prim.FillBoundary(geom.periodicity());
    setBC(prim, cu);

    if (plot_int > 0) {
        WritePlotFile(0, 0.0, geom, cu, cuMeans, cuVars,
                      prim, primMeans, primVars, spatialCross, eta, kappa);
    }

    mui::uniface2d uniface( "mpi://FHD-side/FHD-KMC-coupling" );

    //Time stepping loop
    for(step=1;step<=max_step;++step) {

        if (restart > 0 && step==1) {
            ReadCheckPoint(step, time, statsCount, geom, cu, cuMeans, cuVars, prim,
                           primMeans, primVars, spatialCross, miscStats, eta, kappa);
        }

        // timer
        Real ts1 = ParallelDescriptor::second();

        mui_push(cu, prim, dx, uniface, step);

        RK3step(cu, cup, cup2, cup3, prim, source, eta, zeta, kappa, chi, D, flux,
                stochFlux, cornx, corny, cornz, visccorn, rancorn, geom, dx, dt);

        // set cu_temp to zero since mui_fetch *increments* cu_temp
        cu_temp.setVal(0.);

        // store the *increment* in cu_temp
        mui_fetch(cu_temp, prim, dx, uniface, step);

        // add cu_temp to cu
        MultiFab::Add(cu,cu_temp,0,0,nvars,ngc);

        conservedToPrimitive(prim, cu);

        // Set BC: 1) fill boundary 2) physical
        cu.FillBoundary(geom.periodicity());
        prim.FillBoundary(geom.periodicity());
        setBC(prim, cu);

        // timer
        Real ts2 = ParallelDescriptor::second() - ts1;
        ParallelDescriptor::ReduceRealMax(ts2);
        amrex::Print() << "Advanced step " << step << " in " << ts2 << " seconds\n";

        // timer
        Real aux1 = ParallelDescriptor::second();

        // compute mean and variances
        if (step > n_steps_skip) {
            evaluateStats(cu, cuMeans, cuVars, prim, primMeans, primVars,
                          spatialCross, miscStats, miscVals, statsCount, dx);
            statsCount++;
        }

        // write a plotfile
        if (plot_int > 0 && step > 0 && step%plot_int == 0) {
/*
           yzAverage(cuMeans, cuVars, primMeans, primVars, spatialCross,
                     cuMeansAv, cuVarsAv, primMeansAv, primVarsAv, spatialCrossAv);
           WritePlotFile(step, time, geom, cu, cuMeansAv, cuVarsAv,
                         prim, primMeansAv, primVarsAv, spatialCrossAv, eta, kappa);
*/

           // also horizontal average
           WriteHorizontalAverage(cu,2,5,nspecies,step,geom);
        }

        if (chk_int > 0 && step > 0 && step%chk_int == 0)
        {
           WriteCheckPoint(step, time, statsCount, geom, cu, cuMeans,
                           cuVars, prim, primMeans, primVars, spatialCross, miscStats, eta, kappa);
        }

        // collect a snapshot for structure factor
        if (step > n_steps_skip && struct_fact_int > 0 && (step-n_steps_skip)%struct_fact_int == 0) {
           MultiFab::Copy(structFactPrimMF, prim, 0,                0,                structVarsPrim,   0);
           MultiFab::Copy(structFactConsMF, cu,   0,                0,                structVarsCons-1, 0);
           MultiFab::Copy(structFactConsMF, prim, AMREX_SPACEDIM+1, structVarsCons-1, 1,                0); // temperature too
           structFactPrim.FortStructure(structFactPrimMF);
           structFactCons.FortStructure(structFactConsMF);
           if(project_dir >= 0) {
                ComputeVerticalAverage(prim, primVertAvg, project_dir, 0, structVarsPrim);
                structFactPrimVerticalAverage.FortStructure(primVertAvg);
           }
        }

        // write out structure factor
        if (step > n_steps_skip && struct_fact_int > 0 && plot_int > 0 && step%plot_int == 0) {
            structFactPrim.WritePlotFile(step,time,geom,"plt_SF_prim");
            structFactCons.WritePlotFile(step,time,geom,"plt_SF_cons");
            if(project_dir >= 0) {
                structFactPrimVerticalAverage.WritePlotFile(step,time,geom_flat,"plt_SF_prim_VerticalAverage");
            }
        }

        // timer
        Real aux2 = ParallelDescriptor::second() - aux1;
        ParallelDescriptor::ReduceRealMax(aux2);
        amrex::Print() << "Aux time (stats, struct fac, plotfiles) " << aux2 << " seconds\n";

        time = time + dt;

        // MultiFab memory usage
        const int IOProc = ParallelDescriptor::IOProcessorNumber();

        amrex::Long min_fab_megabytes  = amrex::TotalBytesAllocatedInFabsHWM()/1048576;
        amrex::Long max_fab_megabytes  = min_fab_megabytes;

        ParallelDescriptor::ReduceLongMin(min_fab_megabytes, IOProc);
        ParallelDescriptor::ReduceLongMax(max_fab_megabytes, IOProc);

        amrex::Print() << "High-water FAB megabyte spread across MPI nodes: ["
                       << min_fab_megabytes << " ... " << max_fab_megabytes << "]\n";

        min_fab_megabytes  = amrex::TotalBytesAllocatedInFabs()/1048576;
        max_fab_megabytes  = min_fab_megabytes;

        ParallelDescriptor::ReduceLongMin(min_fab_megabytes, IOProc);
        ParallelDescriptor::ReduceLongMax(max_fab_megabytes, IOProc);

        amrex::Print() << "Curent     FAB megabyte spread across MPI nodes: ["
                       << min_fab_megabytes << " ... " << max_fab_megabytes << "]\n";
    }

    // timer
    Real stop_time = ParallelDescriptor::second() - strt_time;
    ParallelDescriptor::ReduceRealMax(stop_time);
    amrex::Print() << "Run time = " << stop_time << std::endl;
}