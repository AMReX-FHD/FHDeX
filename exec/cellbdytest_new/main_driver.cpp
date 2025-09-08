#include "INS_functions.H"



#include "species.H"
#include "paramPlane.H"

#include "StructFact.H"

#include "StochMomFlux.H"

#include "hydro_functions.H"

#include "electrostatic.H"

#include "particle_functions.H"

#include "chrono"

using namespace std::chrono;

// argv contains the name of the inputs file entered at the command line
void main_driver(const char* argv)
{
    // timer for total simulation time
    Real strt_time = ParallelDescriptor::second();

    std::string inputs_file = argv;


    // copy contents of F90 modules to C++ namespaces
    InitializeCommonNamespace();
    InitializeGmresNamespace();

    int step = 1;
    Real time = 0.;
    int statsCount = 1;

    /*
      Terms prepended with a 'C' are related to the particle grid; only used for finding neighbor lists
      Those with 'P' are for the electostatic grid.
      Those without are for the fluid grid.
      The particle grid and es grids are created as a corsening or refinement of the fluid grid.
    */

    // BoxArray for the fluid
    BoxArray ba;

    // BoxArray for the particles
    BoxArray bc;

    // BoxArray electrostatic grid
    BoxArray bp;

    // Box for the fluid
    IntVect dom_lo(AMREX_D_DECL(           0,            0,            0));
    IntVect dom_hi(AMREX_D_DECL(n_cells[0]-1, n_cells[1]-1, n_cells[2]-1));
    Box domain(dom_lo, dom_hi);

    // how boxes are distrubuted among MPI processes
    DistributionMapping dmap;

    // set number of ghost cells to fit whole peskin kernel
    // AJN - for perdictor/corrector do we need one more ghost cell if the predictor pushes
    //       a particle into a ghost region?
    int ang = 1;
    for(int i=0;i<nspecies;i++)
    {
        int tempang = 1;
        if (pkernel_fluid[i] == 3) {
            tempang = 2;
        }
        else if (pkernel_fluid[i] == 4) {
            tempang = 3;
        }
        else if (pkernel_fluid[i] == 6) {
            tempang = 4;
        }
        else if (eskernel_fluid[i] == 4) {
            tempang = 3;
        }
        else if (eskernel_fluid[i] == 5) {
            tempang = 3;
        }
        else if (eskernel_fluid[i] == 6) {
            tempang = 4;
        }
        else {
            tempang = floor(eskernel_fluid[i]/2)+1;
        }

        if(tempang > ang)
        {
            ang = tempang;
        }
    }
    //if (bond_tog != 0) {
    //    ang = std::max(ang, 8);
    //}

    int ngp = 1;
    // using maximum number of peskin kernel points to determine the ghost cells for the whole grid.
    //     not sure if it will cause problem for BCs.
    if (*(std::max_element(pkernel_es.begin(),pkernel_es.begin()+nspecies)) == 3) {
        ngp = 2;
    }
    else if (*(std::max_element(pkernel_es.begin(),pkernel_es.begin()+nspecies)) == 4) {
        ngp = 3;
    }
    else if (*(std::max_element(pkernel_es.begin(),pkernel_es.begin()+nspecies)) == 6) {
        ngp = 4;
    }

    //// TODO: need a better way to determine ghost cells for bonds
    //if (bond_tog != 0) {
    //    ngp = std::max(ngp, 6);
    //}

    // staggered velocities
    // umac needs extra ghost cells for Peskin kernels
    // note if we are restarting, these are defined and initialized to the checkpoint data
    std::array< MultiFab, AMREX_SPACEDIM > umac;
    std::array< MultiFab, AMREX_SPACEDIM > umacM;    // mean

    std::array< MultiFab, AMREX_SPACEDIM > touched;

    // pressure for GMRES solve; 1 ghost cell
    MultiFab pres;

    // MFs for storing particle statistics
    // A lot of these relate to gas kinetics, but many are still useful so leave in for now.
    MultiFab particleMeans;
    MultiFab particleVars;

    // MF for electric potential
    MultiFab potential;
    MultiFab potentialM;

    // MF for charge mean and variance
    MultiFab chargeM;

    // variables that save FFT (real and imag parts) at t0
    MultiFab struct_cc_numdens0_real;
    MultiFab struct_cc_numdens0_imag;

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

        // Initialize the boxarray "ba" from the single box "bx"
        ba.define(domain);

        // Break up boxarray "ba" into chunks no larger than "max_grid_size" along a direction
        // note we are converting "Vector<int> max_grid_size" to an IntVect
        ba.maxSize(IntVect(max_grid_size));

        // how boxes are distrubuted among MPI processes
        dmap.define(ba);

        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            umac [d].define(convert(ba,nodal_flag_dir[d]), dmap, 1, ang);
            touched[d].define(convert(ba,nodal_flag_dir[d]), dmap, 1, ang);
            umacM[d].define(convert(ba,nodal_flag_dir[d]), dmap, 1, 1);
            umac [d].setVal(0.);
            umacM[d].setVal(0.);
        }

        pres.define(ba,dmap,1,1);
        pres.setVal(0.);

        bc = ba;
        bp = ba;

        // particle grid_refine: <1 = refine, >1 = coarsen.
        // assume only powers of 2 for now
        if (particle_grid_refine < 1) {
            int sizeRatio = (int)(1.0/particle_grid_refine);
            bc.refine(sizeRatio);
        }
        else {
            int sizeRatio = (int)(particle_grid_refine);
            bc.coarsen(sizeRatio);
        }

        if (es_grid_refine < 1) {
            int sizeRatio = (int)(1.0/es_grid_refine);
            bp.refine(sizeRatio);
        }
        else {
            int sizeRatio = (int)(es_grid_refine);
            bp.coarsen(sizeRatio);
        }

        // Variables (C++ index)
        // ( 0) Members
        // ( 1) Density
        // ( 2) velx
        // ( 3) vely
        // ( 4) velz
        // ( 5) Temperature
        // ( 6) jx
        // ( 7) jy
        // ( 8) jz
        // ( 9) energyDensity
        // (10) pressure
        // (11) Cx
        // (12) Cy
        // (13) Cz
        particleMeans.define(bc, dmap, 8+nspecies, 0);
        particleMeans.setVal(0.);

        // Variables (C++ index)
        // ( 0) Members
        // ( 1) Density
        // ( 2) velx
        // ( 3) vely
        // ( 4) velz
        // ( 5) Temperature
        // ( 6) jx
        // ( 7) jy
        // ( 8) jz
        // ( 9) energyDensity
        // (10) pressure
        // (11) GVar
        // (12) KGCross
        // (13) KRhoCross
        // (14) RhoGCross
        // (15) Cx
        // (16) Cy
        // (17) Cz

        //Cell centred es potential
        potential.define(bp, dmap, 1, ngp);
        potentialM.define(bp, dmap, 1, 1);
        potential.setVal(0.);
        potentialM.setVal(0.);

        chargeM.define(bp, dmap, 1, 1);  // mean
        chargeM.setVal(0);

        struct_cc_numdens0_real.define(bc, dmap, 1, 0);
        struct_cc_numdens0_imag.define(bc, dmap, 1, 0);
    }
    else {

        // restart from checkpoint
        ReadCheckPoint(step,time,statsCount,umac,umacM,pres,
                       particleMeans,particleVars,chargeM,
                       potential,potentialM,
                       struct_cc_numdens0_real,struct_cc_numdens0_imag);

        // grab DistributionMap from umac
        dmap = umac[0].DistributionMap();

        // grab fluid BoxArray from umac and convert to cell-centered
        ba = umac[0].boxArray();
        ba.enclosedCells();

        // grab particle BoxArray from particleMeans
        bc = particleMeans.boxArray();

        // grab electrostatic potential BoxArray from potential
        bp = potential.boxArray();
    }

    // Domain boxes for particle and electrostatic grids
    Box domainC = domain;
    Box domainP = domain;

    // particle grid and es grid_refine: <1 = refine, >1 = coarsen.
    // assume only powers of 2 for now
    // note particle grid BoxArray was handled above
    if (particle_grid_refine < 1) {
        int sizeRatio = (int)(1.0/particle_grid_refine);
        domainC.refine(sizeRatio);
    }
    else {
        int sizeRatio = (int)(particle_grid_refine);
        domainC.coarsen(sizeRatio);
    }
    if (es_grid_refine < 1) {
        int sizeRatio = (int)(1.0/es_grid_refine);
        domainP.refine(sizeRatio);
    }
    else {
        int sizeRatio = (int)(es_grid_refine);
        domainP.coarsen(sizeRatio);
    }

    // is the problem periodic?
    Vector<int> is_periodic  (AMREX_SPACEDIM,0);  // set to 0 (not periodic) by default
    Vector<int> is_periodic_c(AMREX_SPACEDIM,0);  // set to 0 (not periodic) by default
    Vector<int> is_periodic_p(AMREX_SPACEDIM,0);  // set to 0 (not periodic) by default
    for (int i=0; i<AMREX_SPACEDIM; ++i) {
        if (bc_vel_lo[i] == -1 && bc_vel_hi[i] == -1) {
            is_periodic  [i] = 1;
            is_periodic_c[i] = 1;
        }
        if (bc_es_lo[i] == -1 && bc_es_hi[i] == -1) {
            is_periodic_p[i] = 1;
        }
    }

    // This defines a Geometry object
    RealBox realDomain({AMREX_D_DECL(prob_lo[0],prob_lo[1],prob_lo[2])},
                       {AMREX_D_DECL(prob_hi[0],prob_hi[1],prob_hi[2])});

    Geometry geom (domain ,&realDomain,CoordSys::cartesian,is_periodic.  data());
    Geometry geomC(domainC,&realDomain,CoordSys::cartesian,is_periodic_c.data());
    Geometry geomP(domainP,&realDomain,CoordSys::cartesian,is_periodic_p.data());

    const Real* dx = geom.CellSize();
    const Real* dxc = geomC.CellSize();
    const Real* dxp = geomP.CellSize();


    Real dt = fixed_dt;
    Real dtinv = 1.0/dt;

    // AJN - get rid of collision stuff?
#if (AMREX_SPACEDIM == 2)
    int totalCollisionCells = (domainC.hiVect()[0]+1)*(domainC.hiVect()[1]+1);
#elif (AMREX_SPACEDIM == 3)
    int totalCollisionCells = (domainC.hiVect()[0]+1)*(domainC.hiVect()[1]+1)*(domainC.hiVect()[2]+1);
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
    double wetRad[nspecies];
    double dxAv = (dx[0] + dx[1] + dx[2])/3.0; //This is probably the wrong way to do this.
    std::string line;

    for(int j=0;j<nspecies;j++) {
       if (pkernel_fluid[j] == 3) {
           wetRad[j] = 0.91*dxAv;
       }
       else if (pkernel_fluid[j] == 4) {
           wetRad[j] = 1.255*dxAv;
       }
       else if (pkernel_fluid[j] == 6) {
           wetRad[j] = 1.481*dxAv;
       }
       else if (pkernel_fluid[j] == 1) {
           wetRad[j] = 1.255*dxAv;
       }
       else if (eskernel_fluid[j] == 4) {
           if (eskernel_beta[j] < 4 || eskernel_beta[j] > 12) {
              Abort("Please provide eskernel_beta within the range [1,3]*eskernel_fluid.");
           }

           int targetLine = ((int)(10*eskernel_beta[j])-10*eskernel_fluid[j])+1;
           std::ifstream wetRad_w4("wetRad_w4.dat");
           for (int lineCount=0; lineCount < targetLine-1; lineCount++) {
               wetRad_w4.ignore(100000, '\n');
           }
           wetRad_w4 >> wetRad[j];
           Print() << "wetRad read from file is " << wetRad[j] << std::endl;
           wetRad[j] *= dxAv;
           wetRad_w4.close();
           //wetRad[j] = 1.300*dxAv; // With beta = 5.22
       }
       else if (eskernel_fluid[j] == 5) {
           if (eskernel_beta[j] < 5 || eskernel_beta[j] > 15) {
              Abort("Please provide eskernel_beta within the range [1,3]*eskernel_fluid.");
           }

           int targetLine = ((int)(10*eskernel_beta[j])-10*eskernel_fluid[j])+1;
           //Print() << targetLine << std::endl;
           std::ifstream wetRad_w5("wetRad_w5.dat");
           for (int lineCount=0; lineCount < targetLine-1; lineCount++) {
               wetRad_w5.ignore(100000, '\n');
           }
           wetRad_w5 >> wetRad[j];
           Print() << "wetRad read from file is " << wetRad[j] << std::endl;
           wetRad[j] *= dxAv;
           wetRad_w5.close();
       }
       else if (eskernel_fluid[j] == 6) {
           if (eskernel_beta[j] < 6 || eskernel_beta[j] > 18) {
              Abort("Please provide eskernel_beta within the range [1,3]*eskernel_fluid.");
           }

           int targetLine = ((int)(10*eskernel_beta[j])-10*eskernel_fluid[j])+1;
           //Print() << targetLine << std::endl;
           std::ifstream wetRad_w6("wetRad_w6.dat");
           for (int lineCount=0; lineCount < targetLine-1; lineCount++) {
               wetRad_w6.ignore(100000, '\n');
           }
           wetRad_w6 >> wetRad[j];
           Print() << "wetRad read from file is " << wetRad[j] << std::endl;
           wetRad[j] *= dxAv;
           wetRad_w6.close();
           //wetRad[j] = 1.478*dxAv; // With beta = 8.64
       }
       else if (eskernel_fluid[j] == 3) {
           if (eskernel_beta[j] < 3 || eskernel_beta[j] > 9) {
              Abort("Please provide eskernel_beta within the range [1,3]*eskernel_fluid.");
           }

           int targetLine = ((int)(10*eskernel_beta[j])-10*eskernel_fluid[j])+1;
           //Print() << targetLine << std::endl;
           std::ifstream wetRad_w3("wetRad_w3.dat");
           for (int lineCount=0; lineCount < targetLine-1; lineCount++) {
               wetRad_w3.ignore(100000, '\n');
           }
           wetRad_w3 >> wetRad[j];
           Print() << "wetRad read from file is " << wetRad[j] << std::endl;
           wetRad[j] *= dxAv;
           wetRad_w3.close();
       }
       else if (eskernel_fluid[j] == 7) {
           if (eskernel_beta[j] < 7 || eskernel_beta[j] > 21) {
              Abort("Please provide eskernel_beta within the range [1,3]*eskernel_fluid.");
           }

           int targetLine = ((int)(10*eskernel_beta[j])-10*eskernel_fluid[j])+1;
           Print() << "TARGET: " << targetLine << std::endl;
           std::ifstream wetRad_w7("wetRad_w7.dat");
           for (int lineCount=0; lineCount < targetLine-1; lineCount++) {
               wetRad_w7.ignore(100000, '\n');
           }
           wetRad_w7 >> wetRad[j];
           Print() << "wetRad read from file is " << wetRad[j] << std::endl;
           wetRad[j] *= dxAv;
           wetRad_w7.close();
       }
       else {
           Abort("Currently the code only supports pkernel_fluid = 1,3,4,6 or eskernel_fluid = 3,4,5,6,7.");
           //wetRad[j] = 1.255*dxAv;
       }
    }

    for(int i=0;i<nspecies;i++) {

        ionParticle[i].m = mass[i];
        ionParticle[i].q = qval[i];

        if (diameter[i] > 0) { // positive diameter

            // set diameter from inputs
            ionParticle[i].d = diameter[i];

            // compute total diffusion from input diameter
            ionParticle[i].totalDiff = (k_B*T_init[0])/(6*M_PI*(diameter[i]/2.0)*visc_coef);

            // compute wet diffusion from wetRad
            ionParticle[i].wetDiff = (k_B*T_init[0])/(6*M_PI*wetRad[i]*visc_coef);

            if (all_dry == 1) {
                ionParticle[i].dryDiff = ionParticle[i].totalDiff;
                ionParticle[i].wetDiff = 0;
            }
            else {
                // dry = total - wet
                ionParticle[i].dryDiff = ionParticle[i].totalDiff - ionParticle[i].wetDiff;
            }

        }
        else { // zero or negative diameter

            // set total diffusion from inputs
            ionParticle[i].totalDiff = diff[i];

            // set diameter from total diffusion (Stokes Einsten)
            ionParticle[i].d = 2.0*(k_B*T_init[0])/(6*M_PI*(ionParticle[i].totalDiff)*visc_coef);

               // std::cout << "Species " << i << " radius: " << ionParticle[i].d << std::endl;

            // compute wet diffusion from wetRad
            ionParticle[i].wetDiff = (k_B*T_init[0])/(6*M_PI*wetRad[i]*visc_coef);

            if (all_dry == 1) {
                ionParticle[i].dryDiff = ionParticle[i].totalDiff;
                ionParticle[i].wetDiff = 0;
            }
            else {
                // dry = total - wet
                ionParticle[i].dryDiff = ionParticle[i].totalDiff - ionParticle[i].wetDiff;
            }
        }

        if (all_dry == 1 && fluid_tog != 0) {
            Abort("Abort: Don't use all_dry=1 and fluid_tog!=0 together");
        }

        Print() << "Species " << i << "\n"
                << " total diffusion: " << ionParticle[i].totalDiff << "\n"
                << " wet diffusion: " << ionParticle[i].wetDiff
                << " percent wet: " << 100.*ionParticle[i].wetDiff/ionParticle[i].totalDiff << "\n"
                << " dry diffusion: " << ionParticle[i].dryDiff
                << " percent dry: " << 100.*ionParticle[i].dryDiff/ionParticle[i].totalDiff << "\n"
                << " total radius: " << ionParticle[i].d/2.0 << "\n"
                << " wet radius: " << wetRad[i] << "\n"
                << " dry radius: " << (k_B*T_init[0])/(6*3.14159265359*(ionParticle[i].dryDiff)*visc_coef) << "\n";

        //if (ionParticle[i].dryDiff < 0) {
        //    Print() << "Negative dry diffusion in species " << i << "\n";
        //    Abort();
        //}

        ionParticle[i].Neff = particle_neff; // From DSMC, this will be set to 1 for electolyte calcs
        ionParticle[i].R = k_B/ionParticle[i].m; //used a lot in kinetic stats cals, bu not otherwise necessary for electrolytes

        ionParticle[i].sigma = sigma[i];
        ionParticle[i].eepsilon = eepsilon[i];

        // round up particles so there are the same number in each box;
        // we have to divide them into whole numbers of particles somehow.
        if (particle_count[i] >= 0) {
            ionParticle[i].ppb = (double)particle_count[i]/(double)ba.size();
            ionParticle[i].total = particle_count[i];
            ionParticle[i].n0 = ionParticle[i].total/domainVol;

            Print() << "Species " << i << " count adjusted to " << ionParticle[i].total << "\n";
        }
        else {
            // if particle count is negative, we instead compute the number of particles based on particle density and particle_neff
            ionParticle[i].total = (int)amrex::Math::ceil(particle_n0[i]*domainVol/particle_neff);
            // adjust number of particles up so there is the same number per box
            ionParticle[i].ppb = (double)ionParticle[i].total/(double)ba.size();
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

    // see the variable list used above above for particleMeans
    MultiFab particleInstant(bc, dmap, 8+nspecies, 0);
    // also initialize a MultiFab to store variables at t0, which is used in dynamic structure factor
    MultiFab particleInstant0(bc, dmap, 8+nspecies, 0);

    //-----------------------------
    //  Hydro setup

    ///////////////////////////////////////////
    // rho, alpha, beta, gamma:
    ///////////////////////////////////////////

    // this only needs 1 ghost cell
    MultiFab rho(ba, dmap, 1, 1);
    rho.setVal(1.);

    // the number of ghost cells needed for the GMRES solve for alpha, beta, etc., is a fixed number
    // and not dependent on the Peskin kernels.
    // alpha_fc -> 1 ghost cell
    // beta -> 1
    // beta_ed -> 1
    // gamma -> 1

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

    // the number of ghost cells needed here is
    // eta_cc -> 1
    // temp_cc -> 1
    // eta_ed -> 0
    // temp_ed -> 0

    // eta and temperature; cell-centered
    MultiFab  eta_cc(ba, dmap, 1, 1);
    MultiFab temp_cc(ba, dmap, 1, 1);

    // eta and temperature; nodal
    std::array< MultiFab, NUM_EDGE >  eta_ed;
    std::array< MultiFab, NUM_EDGE > temp_ed;
#if (AMREX_SPACEDIM == 2)
    eta_ed [0].define(convert(ba,nodal_flag),    dmap, 1, 0);
    temp_ed[0].define(convert(ba,nodal_flag),    dmap, 1, 0);
#elif (AMREX_SPACEDIM == 3)
    eta_ed [0].define(convert(ba,nodal_flag_xy), dmap, 1, 0);
    eta_ed [1].define(convert(ba,nodal_flag_xz), dmap, 1, 0);
    eta_ed [2].define(convert(ba,nodal_flag_yz), dmap, 1, 0);
    temp_ed[0].define(convert(ba,nodal_flag_xy), dmap, 1, 0);
    temp_ed[1].define(convert(ba,nodal_flag_xz), dmap, 1, 0);
    temp_ed[2].define(convert(ba,nodal_flag_yz), dmap, 1, 0);
#endif

    // Initalize eta & temperature multifabs
    eta_cc .setVal(eta_const);
    temp_cc.setVal(temp_const);
    for (int d=0; d<NUM_EDGE; ++d) {
        eta_ed[d] .setVal(eta_const);
        temp_ed[d].setVal(temp_const);
    }
    ///////////////////////////////////////////

    ///////////////////////////////////////////
    // random fluxes:
    ///////////////////////////////////////////

    // mflux divergence, staggered in x,y,z
    // Define mfluxdiv predictor/corrector multifabs
    std::array< MultiFab, AMREX_SPACEDIM >  stochMfluxdiv;
    std::array< MultiFab, AMREX_SPACEDIM >  stochMfluxdivC;
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        stochMfluxdiv [d].define(convert(ba,nodal_flag_dir[d]), dmap, 1, 0);
        stochMfluxdivC[d].define(convert(ba,nodal_flag_dir[d]), dmap, 1, 0);
        stochMfluxdiv [d].setVal(0.0);
        stochMfluxdivC[d].setVal(0.0);
    }

    // Declare object of StochMomFlux class
    int n_rngs = 1; // we only need 1 stage of random numbers
    StochMomFlux sMflux (ba,dmap,geom,n_rngs);

    // weights for random number stages
    Vector< amrex::Real> weights;
    weights = {1.0};

    ///////////////////////////////////////////

    // additional staggered velocity MultiFabs
    std::array< MultiFab, AMREX_SPACEDIM > umacNew;
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        umacNew[d].define(convert(ba,nodal_flag_dir[d]), dmap, 1, ang);
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

    std::array< MultiFab, AMREX_SPACEDIM > sourceRFD;
    AMREX_D_TERM(sourceRFD[0].define(convert(ba,nodal_flag_x), dmap, 1, ang);,
                 sourceRFD[1].define(convert(ba,nodal_flag_y), dmap, 1, ang);,
                 sourceRFD[2].define(convert(ba,nodal_flag_z), dmap, 1, ang););

    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        source    [d].setVal(0.0);
        sourceTemp[d].setVal(0.0);
        sourceRFD[d].setVal(0.0);
        touched[d].setVal(0.0);
    }

    // simple shear
    std::array< MultiFab, AMREX_SPACEDIM > externalV;
    AMREX_D_TERM(externalV[0].define(convert(ba,nodal_flag_x), dmap, 1, ang);,
                 externalV[1].define(convert(ba,nodal_flag_y), dmap, 1, ang);,
                 externalV[2].define(convert(ba,nodal_flag_z), dmap, 1, ang););
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        externalV[d].setVal(0.0);  // simple shear
    }

    // uniform velocity due to shear
    std::array< MultiFab, AMREX_SPACEDIM > externalVU;
    AMREX_D_TERM(externalVU[0].define(convert(ba,nodal_flag_x), dmap, 1, ang);,
                 externalVU[1].define(convert(ba,nodal_flag_y), dmap, 1, ang);,
                 externalVU[2].define(convert(ba,nodal_flag_z), dmap, 1, ang););
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        externalVU[d].setVal(0.0);  // simple shear
    }
    //Define parametric paramplanes for particle interaction - declare array for paramplanes and then define properties in BuildParamplanes


    int paramPlaneCount = 6;
    // Make paramPlaneList GPU compatible
    std::ifstream planeFile("paramplanes.dat");
    int fileCount = 0;
    if(planeFile.good())
    {
        planeFile >> fileCount;
    }
    paramPlaneCount = paramPlaneCount + fileCount;

    Gpu::ManagedVector<paramPlane> paramPlaneList(paramPlaneCount);
    // Set up a pointer to access data of paramPlaneList, which is used as a normal vector everywhere
    paramPlane* pparamPlaneList = paramPlaneList.data();
    BuildParamplanes(pparamPlaneList,paramPlaneCount,realDomain.lo(),realDomain.hi());

    // IBMarkerContainerBase default behaviour is to do tiling. Turn off here:

    //----------------------
    // Particle tile size
    //----------------------
    Vector<int> ts(BL_SPACEDIM);

    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        if (max_particle_tile_size[d] > 0) {
            ts[d] = max_particle_tile_size[d];
        }
        else {
            ts[d] = max_grid_size[d];
        }
    }

    ParmParse pp ("particles");
    pp.addarr("tile_size", ts);

    //int num_neighbor_cells = 4; replaced by input var
    //Particles! Build on geom & box array for collision cells/ poisson grid?
    double relRefine = particle_grid_refine/es_grid_refine;
    Real max_es_range = 0;
    Real max_sr_range = 0;
    Real max_range = 0;

    for(int i=0;i<nspecies;i++) {
        Real range = (pkernel_es[i] + 0.5)*dxp[0];
        if(range > max_es_range)
        {
           max_es_range = range ;
        }
    }

    for(int i=0;i<(nspecies*nspecies);i++) {
        Real range = sigma[i]*rmax[i];

        if(range > max_sr_range)
        {
           max_sr_range = range ;
        }
    }
//    for(int i=0;i<nspecies;i++) {
//        Real range = sigma_wall[i]*rmax_wall[i];

//        if(range > max_sr_range)
//        {
//           max_sr_range = range ;
//        }
//    }

    if(max_sr_range > max_es_range)
    {
        max_range = max_sr_range;
    }else
    {
        max_range = max_es_range;
    }

    int cRange = (int)ceil(max_range/dxc[0]);

    FhdParticleContainer particles(geomC, geom, dmap, bc, ba, cRange, ang);

    if (restart < 0 && particle_restart < 0) {
        // create particles
        if (false) {
            //particles.InitParticlesFromFile(ionParticle, dxp);
        }
        else {
            particles.InitParticles(ionParticle, dxp);
        }
    }
    else {
        ReadCheckPointParticles(particles, ionParticle, dxp);
    }

    particles.UpdatePIDMap();

    for (int i=0; i<simParticles; i++) {
        Print() << "id_global and id_pulldown " << i << " " << particles.id_global_map[i] << std::endl;
    }

    //Find coordinates of cell faces (fluid grid). May be used for interpolating fields to particle locations
    FindFaceCoords(RealFaceCoords, geom); //May not be necessary to pass Geometry?

    //Calculate external simple shear flow, used to update umac due to shear
    //AllPrint() << "pin_y " << particles.getPinnedY() << std::endl;
    SimpleShear(externalV, geom, RealFaceCoords, particles.getPinnedY());
    UniformVel(externalVU, geom);
    //for (MFIter mfi(externalV[0]); mfi.isValid(); ++mfi)
    //{
    //    Array4<Real> const& xvel = (externalV[0]).array(mfi);
    //    AllPrint() << "xvel below pinned particles from MF is " << xvel(4,16,4) << std::endl;
    //    AllPrint() << "xvel above pinned particles from MF is " << xvel(4,48,4) << std::endl;
    //}
    Print() << "Pinned particles y location is " << particles.getPinnedY() << std::endl;
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        // record actual shear speed, used to update particle velocity due to shear
        particles.wallspeed_x_lo_real[d] = wallspeed_x_lo[d];
        particles.wallspeed_y_lo_real[d] = wallspeed_y_lo[d];
        particles.wallspeed_z_lo_real[d] = wallspeed_z_lo[d];
        particles.wallspeed_x_hi_real[d] = wallspeed_x_hi[d];
        particles.wallspeed_y_hi_real[d] = wallspeed_y_hi[d];
        particles.wallspeed_z_hi_real[d] = wallspeed_z_hi[d];
        // reset shear speed to 0 before GMRES solver
        wallspeed_x_lo[d] = 0.;
        wallspeed_y_lo[d] = 0.;
        wallspeed_z_lo[d] = 0.;
        wallspeed_x_hi[d] = 0.;
        wallspeed_y_hi[d] = 0.;
        wallspeed_z_hi[d] = 0.;
    }

    //----------------------
    // Electrostatic setup
    //----------------------

    // cell centered real coordinates - es grid
    MultiFab RealCenteredCoords;
    RealCenteredCoords.define(bp, dmap, AMREX_SPACEDIM, ngp);

    FindCenterCoords(RealCenteredCoords, geomP);

    //charage density for RHS of Poisson Eq.
    MultiFab charge(bp, dmap, 1, ngp);
    charge.setVal(0);

    // temporary used to help collect charge
    MultiFab chargeTemp(bp, dmap, 1, ngp);
    chargeTemp.setVal(0);

    //mass density on ES grid - not necessary
    MultiFab massFrac(bp, dmap, 1, 1);
    MultiFab massFracTemp(bp, dmap, 1, 1);
    massFrac.setVal(0);
    massFracTemp.setVal(0);

    //Staggered electric fields
    std::array< MultiFab, AMREX_SPACEDIM > efield;
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        efield[d].define(convert(bp,nodal_flag_dir[d]), dmap, 1, ngp);
    }

    //Centred electric fields (using an array instead of a multi-component MF)
    std::array< MultiFab, AMREX_SPACEDIM > efieldCC;
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        efieldCC[d].define(bp, dmap, 1, ngp);
        efieldCC[d].setVal(0.);
    }

    // external field
    std::array< MultiFab, AMREX_SPACEDIM > external;
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        external[d].define(bp, dmap, 1, ngp);
    }

    ///////////////////////////////////////////
    // structure factor for numdens-numdens
    ///////////////////////////////////////////

    // names of variables in struct_cc_dens
    Vector< std::string > var_names_dens(1);
    var_names_dens[0] = "dens";

    // multifab that stores real t0 data on which FFT will be performed
    MultiFab struct_cc_numdens0;
    struct_cc_numdens0.define(bc, dmap, 1, 0);

    // variables that store FFT at current time t
    MultiFab struct_cc_numdens;
    struct_cc_numdens.define(bc, dmap, 1, 0);

    // these are the number of pairs we want to report
    int nvar_sf_numdens = 1;
    amrex::Vector< int > s_pairA_numdens(nvar_sf_numdens);
    amrex::Vector< int > s_pairB_numdens(nvar_sf_numdens);

    // Select which variable pairs to include in structure factor:
    s_pairA_numdens[0] = 0; // numdens-numdens
    s_pairB_numdens[0] = 0;

    Vector<Real> scaling_numdens(nvar_sf_numdens);
    for (int i=0; i<nvar_sf_numdens; ++i) {
        scaling_numdens[i] = 1.;
    }

    StructFact structFact_numdens(bc,dmap,var_names_dens,scaling_numdens,
                                 s_pairA_numdens,s_pairB_numdens);

    //// number of pairs for numdens-numdens, used for copying particle stats multifab
    //// declaration of StructFact is inside FhdParticleContainer::DynStructFact
    //int nvar_sf_numdens = 1;


    ///////////////////////////////////////////
    // structure factor for charge-charge
    ///////////////////////////////////////////

    // names of variables in struct_cc_charge
    Vector< std::string > var_names_charge(1);
    var_names_charge[0] = "charge";

    // variables we want to FFT
    MultiFab struct_cc_charge;
    struct_cc_charge.define(bp, dmap, 1, 0);

    // these are the number of pairs we want to report
    int nvar_sf_charge = 1;
    amrex::Vector< int > s_pairA_charge(nvar_sf_charge);
    amrex::Vector< int > s_pairB_charge(nvar_sf_charge);

    // Select which variable pairs to include in structure factor:
    s_pairA_charge[0] = 0; // charge-charge
    s_pairB_charge[0] = 0;

    Vector<Real> scaling_charge(nvar_sf_charge);
    for (int i=0; i<nvar_sf_charge; ++i) {
        scaling_charge[i] = 1.;
    }

    StructFact structFact_charge(bp,dmap,var_names_charge,scaling_charge,
                                 s_pairA_charge,s_pairB_charge);

    ///////////////////////////////////////////
    // structure factor for vel-vel
    ///////////////////////////////////////////

    // names of variables in struct_cc_vel
    Vector< std::string > var_names_vel(AMREX_SPACEDIM);
    var_names_vel[0] = "xvel";
    var_names_vel[1] = "yvel";
    var_names_vel[2] = "zvel";

    // variables we want to FFT
    MultiFab struct_cc_vel;
    struct_cc_vel.define(bp, dmap, AMREX_SPACEDIM, 0);

    // these are the number of pairs we want to report
    int nvar_sf_vel = 6;
    amrex::Vector< int > s_pairA_vel(nvar_sf_vel);
    amrex::Vector< int > s_pairB_vel(nvar_sf_vel);

    // Select which variable pairs to include in structure factor:
    s_pairA_vel[0] = 0; // xvel-xvel
    s_pairB_vel[0] = 0;
    s_pairA_vel[1] = 0; // xvel-yvel
    s_pairB_vel[1] = 1;
    s_pairA_vel[2] = 0; // xvel-zvel
    s_pairB_vel[2] = 2;
    s_pairA_vel[3] = 1; // yvel-yvel
    s_pairB_vel[3] = 1;
    s_pairA_vel[4] = 1; // yvel-zvel
    s_pairB_vel[4] = 2;
    s_pairA_vel[5] = 2; // zvel-zvel
    s_pairB_vel[5] = 2;

    Vector<Real> scaling_vel(nvar_sf_vel);
    for (int i=0; i<nvar_sf_vel; ++i) {
        scaling_vel[i] = 1.;
    }

    StructFact structFact_vel(ba,dmap,var_names_vel,scaling_vel,
                              s_pairA_vel,s_pairB_vel);

    // vectors used for dsf
    Gpu::ManagedVector<Real> posx;
    posx.resize(simParticles);
    Real * posxPtr = posx.dataPtr();
    Vector<Real> posxVec(simParticles);
    Vector<Real> axVec(simParticles);

    Gpu::ManagedVector<Real> posy;
    posy.resize(simParticles);
    Real * posyPtr = posy.dataPtr();
    Vector<Real> posyVec(simParticles);
    Vector<Real> ayVec(simParticles);

    Gpu::ManagedVector<Real> posz;
    posz.resize(simParticles);
    Real * poszPtr = posz.dataPtr();
    Vector<Real> poszVec(simParticles);
    Vector<Real> azVec(simParticles);


//    // Writes instantaneous flow field and some other stuff? Check with Guy.
    //WritePlotFileHydro(0, time, geom, umac, pres, umacM);
    WritePlotFile(0, time, geom, geomC, geomP,
                  particleInstant, particleMeans, particles,
                  charge, chargeM, potential, potentialM, efieldCC);

    remove("bulkFlowEst");
    remove("partPos");
    //Time stepping loop


    int init_step = step;
    if(ramp_step==2){
        dt = dt*1e-7;
        if(step < 100 && step > 50) dt = dt*2;
        if(step < 150 && step > 100) dt = dt*4;
        if(step < 200 && step > 150) dt = dt*8;
        if(step < 250 && step > 200) dt = dt*16;
        if(step < 300 && step > 250) dt = dt*32;
        if(step < 350 && step > 300) dt = dt*64;
        if(step < 400 && step > 350) dt = dt*128;
        if(step < 500 && step > 400) dt = dt*128*sqrt(5);
        if(step < 600 && step > 500) dt = dt*640;
        if(step < 700 && step > 600) dt = dt*640*sqrt(5);
        if(step < 800 && step > 700) dt = dt*3200;
        if(step < 900 && step > 800) dt = dt*3200*sqrt(5);
        if(step < 1000 && step > 900) dt = dt*16000;
        if(step < 1100 && step > 1000) dt = dt*16000*sqrt(5);
        if(step < 1200 && step > 1100) dt = dt*80000;
        if(step < 1300 && step > 1200) dt = dt*80000*sqrt(5);
        if(step < 1400 && step > 1300) dt = dt*400000;
        if(step < 1500 && step > 1400) dt = dt*400000*sqrt(5);
        if(step < 1600 && step > 1500) dt = dt*2000000;
        if(step < 1700 && step > 1600) dt = dt*2000000*sqrt(5);
        if(step > 1700) dt = dt*10000000;

    }else if(ramp_step==1){
        dt = dt*1e-6;
        if(step < 40 && step > 20) dt = dt*2;
        if(step < 80 && step > 40) dt = dt*10;
        if(step < 160 && step > 80) dt = dt*50;
        if(step < 240 && step > 160) dt = dt*100;
        if(step < 320 && step > 240) dt = dt*500;
        if(step < 400 && step > 320) dt = dt*1000;
        if(step < 480 && step > 400) dt = dt*5000;
        if(step < 560 && step > 480) dt = dt*10000;
        if(step < 640 && step > 560) dt = dt*50000;
        if(step < 720 && step > 640) dt = dt*100000;
        if(step < 800 && step > 720) dt = dt*500000;
        if(step > 800) dt = dt*1000000;

    }else{
        dt = dt*1e-5;
        if(step < 40 && step > 20) dt = dt*10;
        if(step < 60 && step > 40) dt = dt*100;
        if(step < 80 && step > 60) dt = dt*1000;
        if(step < 100 && step > 80) dt = dt*10000;
        if(step > 100) dt = dt*100000;
    }



    particles.initRankLists(simParticles);

    Real init_time = ParallelDescriptor::second() - strt_time;
    ParallelDescriptor::ReduceRealMax(init_time);
    amrex::Print() << "Initialization time = " << init_time << " seconds " << std::endl;

    for (int istep=step; istep<=max_step; ++istep) {

        // timer for time step
        Real time1 = ParallelDescriptor::second();

        if(ramp_step==2)
        {
            if(istep == 50)
            {
                dt = dt*2;

                Print() << "\n\nNew dt: " << dt << std::endl<< std::endl<< std::endl;
            }

            if(istep == 100)
            {
                dt = dt*2;

                Print() << "\n\nNew dt: " << dt << std::endl<< std::endl<< std::endl;
            }

            if(istep == 150)
            {
                dt = dt*2;

                Print() << "\n\nNew dt: " << dt << std::endl<< std::endl<< std::endl;
            }


            if(istep == 200)
            {
                dt = dt*2;

                Print() << "\n\nNew dt: " << dt << std::endl<< std::endl<< std::endl;
            }

            if(istep == 250)
            {
                dt = dt*2;

                Print() << "\n\nNew dt: " << dt << std::endl<< std::endl<< std::endl;
            }

            if(istep == 300)
            {
                dt = dt*2;

                Print() << "\n\nNew dt: " << dt << std::endl<< std::endl<< std::endl;
            }
            if(istep == 350)
            {
                dt = dt*2;

                Print() << "\n\nNew dt: " << dt << std::endl<< std::endl<< std::endl;
            }

            if(istep == 400)
            {
                dt = dt*sqrt(5);

                Print() << "\n\nNew dt: " << dt << std::endl<< std::endl<< std::endl;

            }

            if(istep == 500)
            {
                dt = dt*sqrt(5);

                Print() << "\n\nNew dt: " << dt << std::endl<< std::endl<< std::endl;

            }

            if(istep == 600)
            {
                dt = dt*sqrt(5);

                Print() << "\n\nNew dt: " << dt << std::endl<< std::endl<< std::endl;

            }

            if(istep == 700)
            {
                dt = dt*sqrt(5);

                Print() << "\n\nNew dt: " << dt << std::endl<< std::endl<< std::endl;

            }

            if(istep == 800)
            {
                dt = dt*sqrt(5);

                Print() << "\n\nNew dt: " << dt << std::endl<< std::endl<< std::endl;
            }

            if(istep == 900)
            {
                dt = dt*sqrt(5);

                Print() << "\n\nNew dt: " << dt << std::endl<< std::endl<< std::endl;
            }
            if(istep == 1000)
            {
                dt = dt*sqrt(5);

                Print() << "\n\nNew dt: " << dt << std::endl<< std::endl<< std::endl;
            }
            if(istep == 1100)
            {
                dt = dt*sqrt(5);

                Print() << "\n\nNew dt: " << dt << std::endl<< std::endl<< std::endl;
            }

            if(istep == 1200)
            {
                dt = dt*sqrt(5);

                Print() << "\n\nNew dt: " << dt << std::endl<< std::endl<< std::endl;
            }if(istep == 1300)
            {
                dt = dt*sqrt(5);

                Print() << "\n\nNew dt: " << dt << std::endl<< std::endl<< std::endl;
            }if(istep == 1400)
            {
                dt = dt*sqrt(5);
                Print() << "\n\nNew dt: " << dt << std::endl<< std::endl<< std::endl;
            }if(istep == 1500)
            {
                dt = dt*sqrt(5);
                Print() << "\n\nNew dt: " << dt << std::endl<< std::endl<< std::endl;
            }if(istep == 1600)
            {
                dt = dt*sqrt(5);
                Print() << "\n\nNew dt: " << dt << std::endl<< std::endl<< std::endl;
            }if(istep == 1700)
            {
                dt = dt*sqrt(5);
                Print() << "\n\nNew dt: " << dt << std::endl<< std::endl<< std::endl;
            }
        }else if(ramp_step==1)
        {
            if(istep == 20)
            {
                dt = dt*2;
                Print() << "\n\nNew dt: " << dt << std::endl<< std::endl<< std::endl;
            }

            if(istep == 40)
            {
                dt = dt*5;
                Print() << "\n\nNew dt: " << dt << std::endl<< std::endl<< std::endl;
            }

            if(istep == 80)
            {
                dt = dt*5;
                Print() << "\n\nNew dt: " << dt << std::endl<< std::endl<< std::endl;
            }


            if(istep == 160)
            {
                dt = dt*2;
                Print() << "\n\nNew dt: " << dt << std::endl<< std::endl<< std::endl;
            }

            if(istep == 240)
            {
                dt = dt*5;
                Print() << "\n\nNew dt: " << dt << std::endl<< std::endl<< std::endl;
            }

            if(istep == 320)
            {
                dt = dt*2;
                Print() << "\n\nNew dt: " << dt << std::endl<< std::endl<< std::endl;
            }

            if(istep == 400)
            {
                dt = dt*5;
                Print() << "\n\nNew dt: " << dt << std::endl<< std::endl<< std::endl;

            }

            if(istep == 480)
            {
                dt = dt*2;
                Print() << "\n\nNew dt: " << dt << std::endl<< std::endl<< std::endl;

            }


            if(istep == 560)
            {
                dt = dt*5;
                Print() << "\n\nNew dt: " << dt << std::endl<< std::endl<< std::endl;

            }

            if(istep == 640)
            {
                dt = dt*2;
                Print() << "\n\nNew dt: " << dt << std::endl<< std::endl<< std::endl;

            }

            if(istep == 720)
            {
                dt = dt*5;
                Print() << "\n\nNew dt: " << dt << std::endl<< std::endl<< std::endl;
            }

            if(istep == 800)
            {
                dt = dt*2;
                Print() << "\n\nNew dt: " << dt << std::endl<< std::endl<< std::endl;
            }
          }else
        {

            if(istep == 20)
            {
                dt = dt*10;
                Print() << "\n\nNew dt: " << dt << std::endl<< std::endl<< std::endl;
            }
            if(istep == 40)
            {
                dt = dt*10;
                Print() << "\n\nNew dt: " << dt << std::endl<< std::endl<< std::endl;
            }

            if(istep == 60)
            {
                dt = dt*10;
                Print() << "\n\nNew dt: " << dt << std::endl<< std::endl<< std::endl;
            }

            if(istep == 80)
            {
                dt = dt*10;
                Print() << "\n\nNew dt: " << dt << std::endl<< std::endl<< std::endl;
            }

            if(istep == 100)
            {
                dt = dt*10;
                Print() << "\n\nNew dt: " << dt << std::endl<< std::endl<< std::endl;
            }
        }


//            particles.SetPosition(1, prob_hi[0]*0.501, prob_hi[1]*0.501, prob_hi[2]*0.501);
//              particles.SetForce(1,1,0,0);
//          int kk =1;
//          //Print() << "Moving " << ionParticle[0].total << " particles.\n";
//          while(kk<ionParticle[0].total)
//          {
////            Real x1 = prob_lo[0] + 0.25*(prob_hi[0]-prob_lo[0]) + amrex::Random()*(prob_hi[0]-prob_lo[0])/2.0;
//            Real x1 = prob_lo[0] + amrex::Random()*(prob_hi[0]-prob_lo[0]);
//            Real y1 = prob_lo[1] + amrex::Random()*(prob_hi[1]-prob_lo[1]);
//            Real z1 = prob_lo[2] + amrex::Random()*(prob_hi[2]-prob_lo[2]);
////
////            Real seprad = ionParticle[0].d/2.0;
////
////            Real x2 = x1 + 0.0*ionParticle[0].d/2.0;
////            Real y2 = y1;
////            Real z2 = z1;
//////
//            particles.SetPosition(kk,x1 ,y1, z1);
////            particles.SetPosition(kk+1,x2 ,y2, z2);
//
//           // Print() << "Seperating particles " << kk << " and " << kk+1 << " by " << (x1 -x2)/seprad << " radii.\n";
//
//            kk = kk+2;
//

//          }

//          kk =1;
//          while(kk<ionParticle[0].total)
//          {
//            Real xPos[ionParticle[0].total];
//            Real yPos[ionParticle[0].total];
//            Real zPos[ionParticle[0].total];
//
//            particles.PullDown(0, xPos, -1, ionParticle[0].total);
//            particles.PullDown(0, yPos, -2, ionParticle[0].total);
//            particles.PullDown(0, zPos, -3, ionParticle[0].total);
//
//            Real x1 = xPos[kk-1];
//            Real y1 = yPos[kk-1];
//            Real z1 = zPos[kk-1];
////
//            Real seprad = ionParticle[0].d/2.0;
//
//            Real x2 = x1 + 6.0*ionParticle[0].d/2.0;
//            Real y2 = y1;
//            Real z2 = z1;
//////
////            particles.SetPosition(kk,x1 ,y1, z1);
//            particles.SetPosition(kk+1,x2 ,y2, z2);
//
//           // Print() << "Seperating particles " << kk << " and " << kk+1 << " by " << (x1 -x2)/seprad << " radii.\n";
//
//            kk = kk+2;
//

//          }
//
//            x1 = 0.5*prob_hi[0] + (amrex::Random()-0.5)*(prob_hi[0]-prob_lo[0])*0.25;
//            y1 = 0.1875*dxp[0];
//            z1 = 0.5*prob_hi[2] + (amrex::Random()-0.5)*(prob_hi[2]-prob_lo[2])*0.25;


//            Real costheta = 2.*amrex::Random() - 1.;
//            Real sintheta = sqrt(1. - costheta*costheta);

//            Real phi = amrex::Random() * 2. * 3.14159265358979;
//            Real cosphi = std::cos(phi);
//            Real sinphi = std::sin(phi);

//            Real dr = 2*dxp[0];

//            Real x2 = x1 + dr*cosphi;
//            Real y2 = y1;
//            Real z2 = z1 + dr*sinphi;;
//
//            particles.SetPosition(1,x1 ,y1, z1);
//            particles.SetPosition(2,x2 ,y2, z2);


        //Most of these functions are sensitive to the order of execution. We can fix this, but for now leave them in this order.

        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            external[d].setVal(eamp[d]*cos(efreq[d]*time + ephase[d]));  // external field
            source    [d].setVal(body_force_density[d]);      // reset source terms
            sourceTemp[d].setVal(0.0);      // reset source terms
            sourceRFD[d].setVal(0.0);      // reset source terms
            particles.ResetMarkers(0);
            umac[d].setVal(0.0);
        }

        //// subtract shear velocity to previous umac
        //if (istep > 1) {
        //   MultiFab::Subtract(umac[0],externalV[0],0,0,externalV[0].nComp(),externalV[0].nGrow());
        //   MultiFab::Subtract(umac[1],externalV[1],0,0,externalV[1].nComp(),externalV[1].nGrow());
        //   MultiFab::Subtract(umac[2],externalV[2],0,0,externalV[2].nComp(),externalV[2].nGrow());
        //}

        //particles.BuildCorrectionTable(dxp,0);

        if (rfd_tog==1) {
            // Apply RFD force to fluid
            particles.RFD(0, dx, sourceRFD, RealFaceCoords);
            particles.ResetMarkers(0);
        //            particles.DoRFD(dt, dx, dxp, geom, umac, efieldCC, RealFaceCoords, RealCenteredCoords,
        //                            source, sourceTemp, paramPlaneList, paramPlaneCount, 3 /*this number currently does nothing, but we will use it later*/);
        }
        else {
            // set velx/y/z and forcex/y/z for each particle to zero
            particles.ResetMarkers(0);
        }
//        particles.SetForce(1,0.00001,0,0);
//        Real origin[3];
//        origin[0] = prob_hi[0]/2.0;
//        origin[1] = prob_hi[1]/2.0;
//        origin[2] = prob_hi[2]/2.0;

        //particles.forceFunction(dt);

        // sr_tog is short range forces
        // es_tog is electrostatic solve (0=off, 1=Poisson, 2=Pairwise, 3=P3M)


        if (sr_tog != 0 || es_tog==3) {
            // compute short range forces (if sr_tog=1)
            // compute P3M short range correction (if es_tog=3)
            particles.computeForcesNLGPU(charge, RealCenteredCoords, dxp);
        }

        if (es_tog==1 || es_tog==3) {
            // spreads charge density from ions onto multifab 'charge'.
            particles.collectFieldsGPU(dt, dxp, RealCenteredCoords, geomP, charge, chargeTemp, massFrac, massFracTemp);
        }

        // do Poisson solve using 'charge' for RHS, and put potential in 'potential'.
        // Then calculate gradient and put in 'efieldCC', then add 'external'.
        esSolve(potential, charge, efieldCC, external, geomP);

        if (es_tog==2) {
            // compute pairwise Coulomb force (currently hard-coded to work with y-wall).
            particles.computeForcesCoulombGPU(simParticles);
         }

        //particles.computeForcesSpringGPU(simParticles);
        //particles.computeForcesFENEGPU(simParticles);
        if (bond_tog != 0) {
            particles.computeForcesBondGPU(simParticles);
        }

        //particles.SetForce(1,0,0,-1, particles.id_global_map);
        //particles.SetForce(0,0,0,1, particles.id_global_map);

        // compute other forces and spread to grid
        particles.SpreadIonsGPU(dx, dxp, geom, umac, RealFaceCoords, efieldCC, source, sourceTemp);

        //particles.BuildCorrectionTable(dxp,1);

        if ((variance_coef_mom != 0.0) && fluid_tog != 0) {
            // compute the random numbers needed for the stochastic momentum forcing
            sMflux.fillMomStochastic();

            // compute stochastic momentum force
            //sMflux.StochMomFluxDivWideSplit(stochMfluxdiv,0,eta_cc,eta_ed,temp_cc,temp_ed,weights,dt);
            //sMflux.StochMomFluxDivOrder3(stochMfluxdiv,0,eta_cc,eta_ed,temp_cc,temp_ed,weights,dt);
            sMflux.StochMomFluxDiv(stochMfluxdiv,0,eta_cc,eta_ed,temp_cc,temp_ed,weights,dt);

            // integrator containing inertial terms and predictor/corrector requires 2 RNG stages
            if (fluid_tog ==2) {
                sMflux.StochMomFluxDiv(stochMfluxdivC,0,eta_cc,eta_ed,temp_cc,temp_ed,weights,dt);
            }
        }


        MultiFab::Add(source[0],sourceRFD[0],0,0,sourceRFD[0].nComp(),sourceRFD[0].nGrow());
        MultiFab::Add(source[1],sourceRFD[1],0,0,sourceRFD[1].nComp(),sourceRFD[1].nGrow());
        MultiFab::Add(source[2],sourceRFD[2],0,0,sourceRFD[2].nComp(),sourceRFD[2].nGrow());


        if (fluid_tog == 1) {

            if(particles.getTotalPinnedMarkers() != 0)
            {

                Real check;

                /* */
                // Uncomment this section to calculate mobility matrix for pinned particles; this should only run for 1 step
                //   if discos-particle wall is regular, we can just calculate mobility matrix on one particle and shift it around.
                if (pinMatrix_tog == 1)
                {
                    particles.clearMobilityMatrix();
                    for(int ii=0;ii<particles.getTotalPinnedMarkers();ii++)
                    {
                        int id_global_pinned = particles.pinnedParticlesIDGlobal[ii];

                        particles.SetForce(id_global_pinned,1,0,0,particles.id_global_map);
                        for (int d=0; d<AMREX_SPACEDIM; ++d) {
                            source    [d].setVal(0.0);      // reset source terms
                            sourceTemp[d].setVal(0.0);      // reset source terms
                        }
                        particles.SpreadIonsGPU(dx, dxp, geom, umac, RealFaceCoords, efieldCC, source, sourceTemp);
                        advanceStokes(umac,pres,stochMfluxdiv,source,alpha_fc,beta,gamma,beta_ed,geom,dt);
                        particles.InterpolateMarkersGpu(0, dx, umac, RealFaceCoords, check);
                        particles.fillMobilityMatrix(id_global_pinned,0,particles.id_global_map);

                        //particles.PrintParticles();

                        particles.SetForce(id_global_pinned,0,1,0,particles.id_global_map);
                        for (int d=0; d<AMREX_SPACEDIM; ++d) {
                            source    [d].setVal(0.0);      // reset source terms
                            sourceTemp[d].setVal(0.0);      // reset source terms
                        }
                        particles.SpreadIonsGPU(dx, dxp, geom, umac, RealFaceCoords, efieldCC, source, sourceTemp);
                        advanceStokes(umac,pres,stochMfluxdiv,source,alpha_fc,beta,gamma,beta_ed,geom,dt);
                        particles.InterpolateMarkersGpu(0, dx, umac, RealFaceCoords, check);
                        particles.fillMobilityMatrix(id_global_pinned,1,particles.id_global_map);

                        particles.SetForce(id_global_pinned,0,0,1,particles.id_global_map);
                        for (int d=0; d<AMREX_SPACEDIM; ++d) {
                            source    [d].setVal(0.0);      // reset source terms
                            sourceTemp[d].setVal(0.0);      // reset source terms
                        }
                        particles.SpreadIonsGPU(dx, dxp, geom, umac, RealFaceCoords, efieldCC, source, sourceTemp);
                        advanceStokes(umac,pres,stochMfluxdiv,source,alpha_fc,beta,gamma,beta_ed,geom,dt);
                        particles.InterpolateMarkersGpu(0, dx, umac, RealFaceCoords, check);
                        particles.fillMobilityMatrix(id_global_pinned,2,particles.id_global_map);


                    }
                    particles.writeMat();

                    particles.invertMatrix();

                    if(ParallelDescriptor::MyProc() == 0) {
                        Abort("Finish calculating pinned mobility matrix, thus program aborts. To use pinned mobility matrix, set pinMatrix_tog=0 and rerun");
                    }
                }
                /* */

                advanceStokes(umac,pres,stochMfluxdiv,source,alpha_fc,beta,gamma,beta_ed,geom,dt);
                //MultiFab::Add(umac[0],externalVU[0],0,0,externalVU[0].nComp(),externalVU[0].nGrow());
                //MultiFab::Add(umac[1],externalVU[1],0,0,externalVU[1].nComp(),externalVU[1].nGrow());
                //MultiFab::Add(umac[2],externalVU[2],0,0,externalVU[2].nComp(),externalVU[2].nGrow());
                particles.InterpolateMarkersGpu(0, dx, umac, RealFaceCoords, check);
                particles.velNorm();

                //particles.PrintParticles();

                particles.pinnedParticleInversion(particles.id_global_map);

                //particles.PrintParticles();

                for (int d=0; d<AMREX_SPACEDIM; ++d) {
                        source    [d].setVal(0.0);      // reset source terms
                        sourceTemp[d].setVal(0.0);      // reset source terms
                    }

                particles.SpreadIonsGPU(dx, geom, umac, RealFaceCoords, source, sourceTemp);

                MultiFab::Add(source[0],sourceRFD[0],0,0,sourceRFD[0].nComp(),sourceRFD[0].nGrow());
                MultiFab::Add(source[1],sourceRFD[1],0,0,sourceRFD[1].nComp(),sourceRFD[1].nGrow());
                MultiFab::Add(source[2],sourceRFD[2],0,0,sourceRFD[2].nComp(),sourceRFD[2].nGrow());

                advanceStokes(umac,pres,stochMfluxdiv,source,alpha_fc,beta,gamma,beta_ed,geom,dt);
                //MultiFab::Add(umac[0],externalV[0],0,0,externalV[0].nComp(),externalV[0].nGrow());
                //MultiFab::Add(umac[1],externalV[1],0,0,externalV[1].nComp(),externalV[1].nGrow());
                //MultiFab::Add(umac[2],externalV[2],0,0,externalV[2].nComp(),externalV[2].nGrow());
                particles.InterpolateMarkersGpu(0, dx, umac, RealFaceCoords, check);
                particles.velNorm();

                //particles.PrintParticles();

            }else
            {
                // TODO: still missing RFD here?
                advanceStokes(umac,pres,stochMfluxdiv,source,alpha_fc,beta,gamma,beta_ed,geom,dt);
                //MultiFab::Add(umac[0],externalV[0],0,0,externalV[0].nComp(),externalV[0].nGrow());
                //MultiFab::Add(umac[1],externalV[1],0,0,externalV[0].nComp(),externalV[1].nGrow());
                //MultiFab::Add(umac[2],externalV[2],0,0,externalV[0].nComp(),externalV[2].nGrow());

            }

        }
        else if (fluid_tog == 2) {
            Abort("Don't use fluid_tog=2 (inertial Low Mach solver)");
        }

        //Real check;
        //particles.InterpolateMarkersGpu(0, dx, umac, RealFaceCoords,check);
        //particles.TwoParticleCorrelation();

        // total particle move (1=single step, 2=midpoint)
        if (move_tog != 0)
        {
            //Calls wet ion interpolation and movement.

            particles.MoveIonsCPP(dt, dx, dxp, geom, umac, efield, RealFaceCoords, source, sourceTemp, pparamPlaneList,
                               paramPlaneCount, 3 /*this number currently does nothing, but we will use it later*/);


            //particles.PrintParticles();

            // reset statistics after step n_steps_skip
            // if n_steps_skip is negative, we use it as an interval
            if ((n_steps_skip > 0 && istep == n_steps_skip) ||
                (n_steps_skip < 0 && istep%n_steps_skip == 0) ) {

                particles.MeanSqrCalcCM(0, 1);
                particles.stressP(0, 1);
            }
            else {
                particles.MeanSqrCalcCM(0, 0);
                particles.stressP(0, 0);
            }

            Print() << "Finish move.\n";
        }

        /* uncomment this section if using Delong et al method: add shear velocity to umac
         *   also need to add shear velocity in MoveIonsCPP (1-step or 2-step)
        MultiFab::Add(umac[0],externalV[0],0,0,externalV[0].nComp(),externalV[0].nGrow());
        MultiFab::Add(umac[1],externalV[1],0,0,externalV[1].nComp(),externalV[1].nGrow());
        MultiFab::Add(umac[2],externalV[2],0,0,externalV[2].nComp(),externalV[2].nGrow());
        */

        // collect particle positions onto one processor
        particles.GetAllParticlePositions(posxVec,posyVec,poszVec,axVec,ayVec,azVec);
        if (dsf_flag == 1)
        {
           if (ParallelDescriptor::MyProc()==0){
              std::string filename = Concatenate("partPos_restart", init_step, 9);
              std::ofstream ofs2(filename, std::ofstream::app);
              //ofstream ofs( filename, ios::binary );

              for (int i=0;i<simParticles;i++)
              {
                  ofs2 << istep << " " << i << " " << posxVec[i] << " " << posyVec[i] << " " << poszVec[i] << " " << axVec[i] << " " << ayVec[i] << " " << azVec[i] << std::endl;
              }

              ofs2.close();
              //ofs.close();
              Print() << "Finished writing particle positions.\n";
           }
        }

        /*
        // FIXME - AJN
        // NOTE: this stats resetting should eventually be moved to *after* the statsCount++ line below
        // this way, e.g., plot 10 will contain the average of steps 1-10
        // instead of the instantaneous value at step 10
        // however this has the same effect on currentEst so the diagnostics
        // in immersedIons/postprocessing will have to be updated to account for this
        */
        // reset statistics after step n_steps_skip
        // if n_steps_skip is negative, we use it as an interval
        if ((n_steps_skip > 0 && istep == n_steps_skip) ||
            (n_steps_skip < 0 && istep%n_steps_skip == 0) ) {

            particleMeans.setVal(0.0);

            for (int d=0; d<AMREX_SPACEDIM; ++d) {
                umacM[d].setVal(0.);
            }

            Print() << "Resetting stat collection.\n";

            statsCount = 1;
        }

        // g(r)
        if(radialdist_int>0 && istep%radialdist_int == 0) {

            // timer
            Real time_PC1 = ParallelDescriptor::second();

            // compute g(r)
            particles.RadialDistribution(simParticles, istep, ionParticle);
            //particles.potentialDistribution(simParticles, istep, ionParticle);

            // timer
            Real time_PC2 = ParallelDescriptor::second() - time_PC1;
            ParallelDescriptor::ReduceRealMax(time_PC2);
            amrex::Print() << "Time spend computing radial distribution = " << time_PC2 << " seconds" << std::endl;
        }

        // g(x), g(y), g(z)
        if(cartdist_int>0 && istep%cartdist_int == 0) {

            // timer
            Real time_PC1 = ParallelDescriptor::second();

            // compute g(x), g(y), g(z)
            particles.CartesianDistribution(simParticles, istep, ionParticle);

            // timer
            Real time_PC2 = ParallelDescriptor::second() - time_PC1;
            ParallelDescriptor::ReduceRealMax(time_PC2);
            amrex::Print() << "Time spend computing Cartesian distribution = " << time_PC2 << " seconds" << std::endl;
        }

        // compute particle fields, means, anv variances
        // also write out time-averaged current to currentEst
        particles.EvaluateStats(particleInstant, particleMeans, ionParticle[0], dt,statsCount);

        if (dsf_fft) {
           // save t0 stats
           if ((n_steps_skip > 0 && istep == n_steps_skip) ||
               (n_steps_skip < 0 && istep%n_steps_skip == 0) ||
               istep == 1) {

               Print() << "resetting dsf stats at " << istep << " step.\n";
               MultiFab::Copy(struct_cc_numdens0, particleInstant, 0, 0, nvar_sf_numdens, 0);
               structFact_numdens.ComputeFFT(struct_cc_numdens0, struct_cc_numdens0_real, struct_cc_numdens0_imag,geomC);

           }

           // compute dynamic structure factor, using particleInstant and particleInstant0
           MultiFab::Copy(struct_cc_numdens, particleInstant, 0, 0, nvar_sf_numdens, 0);
           structFact_numdens.DynStructureDens(struct_cc_numdens,
                                    struct_cc_numdens0_real, struct_cc_numdens0_imag,
                                geomC, ktarg);

           Print() << "Finish dynamic structure factor.\n";
        }


        // compute the mean and variance of umac
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            ComputeBasicStats(umac[d], umacM[d], 0, 0, statsCount);
        }
        ComputeBasicStats(potential, potentialM, 0, 0, statsCount);
        ComputeBasicStats(charge   , chargeM   , 0, 0, statsCount);

        //Don't forget to add a remove(filename) so it doesn't append to old data
        OutputVolumeMean(umac[0], 0, domainVol, "bulkFlowEst", geom);

        statsCount++;

        //_______________________________________________________________________
        // Update structure factor
        if (struct_fact_int > 0 &&
            istep > amrex::Math::abs(n_steps_skip) &&
            (istep-amrex::Math::abs(n_steps_skip)-1)%struct_fact_int == 0) {

            // charge
            MultiFab::Copy(struct_cc_charge, charge, 0, 0, nvar_sf_charge, 0);
            structFact_charge.FortStructure(struct_cc_charge);

            // velocity
            for (int d=0; d<AMREX_SPACEDIM; ++d) {
                ShiftFaceToCC(umac[d],0,struct_cc_vel,d,1);
            }
            structFact_vel.FortStructure(struct_cc_vel);

            // plot structure factor on plot_int
            if (istep%plot_int == 0) {
                structFact_charge.WritePlotFile(istep,time,"plt_SF_charge");
                structFact_vel   .WritePlotFile(istep,time,"plt_SF_vel");
            }
        }

        // FIXME - AJN: at the moment we are writing out plotfile plot_int-1 also
        // because the time-averaging for the fields resets at n_steps_skip
        // see the FIXME - AJN note above
        bool writePlt = false;
        if (plot_int > 0) {
            if (n_steps_skip >= 0) { // for positive n_steps_skip, write out at plot_int
                writePlt = (istep%plot_int == 0);
            }
            else if (n_steps_skip < 0) { // for negative n_steps_skip, write out at plot_int-1
                writePlt = ((istep+1)%plot_int == 0);
            }
        }
        if (writePlt) {
            // This write particle data and associated fields and electrostatic fields
            WritePlotFile(istep, time, geom, geomC, geomP,
                          particleInstant, particleMeans, particles,
                          charge, chargeM, potential, potentialM, efieldCC);

            // Writes instantaneous flow field and some other stuff? Check with Guy.
            WritePlotFileHydro(istep, time, geom, umac, pres, umacM);
        }

        if (chk_int > 0 && istep%chk_int == 0) {
            WriteCheckPoint(istep, time, statsCount, umac, umacM, pres,
                            particles, particleMeans, particleVars, chargeM,
                            potential, potentialM,
                            struct_cc_numdens0_real, struct_cc_numdens0_imag);
        }

//        particles.PrintParticles();

        // timer for time step
        Real time2 = ParallelDescriptor::second() - time1;
        ParallelDescriptor::ReduceRealMax(time2);
        amrex::Print() << "Advanced step " << istep << " in " << time2 << " seconds\n";

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
    ///////////////////////////////////////////
        //test change
    // timer for total simulation time
    Real stop_time = ParallelDescriptor::second() - strt_time;
    ParallelDescriptor::ReduceRealMax(stop_time);
    amrex::Print() << "Run time = " << stop_time << " seconds" << std::endl;

}
