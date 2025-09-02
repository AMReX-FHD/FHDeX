#include "LocalFunctions.H"
#include "species.H"
#include "paramPlane.H"
#include "StructFact.H"
#include "particle_functions.H"
#include "Checkpoint.H"
#include "chrono"
#include "iostream"
#include "fstream"
#include "DsmcParticleContainer.H"
#include "paramplane_functions_K.H"
#include <AMReX_MultiFab.H>
#include <AMReX_PlotFileUtil.H>

using namespace std::chrono;
using namespace std;

void main_driver(const char* argv)
{
    // timer for total simulation time
    Real strt_time = ParallelDescriptor::second();
    std::string inputs_file = argv;

    InitializeCommonNamespace();

    BoxArray ba;
    IntVect dom_lo(AMREX_D_DECL(           0,            0,            0));
    IntVect dom_hi(AMREX_D_DECL(n_cells[0]-1, n_cells[1]-1, n_cells[2]-1));
    Box domain(dom_lo, dom_hi);
    DistributionMapping dmap;

    int step = 1;
    Real dt = fixed_dt;
    Real time = 0.;
    int statsCount = 1;

    MultiFab cuInst, cuMeans, cuVars;
    MultiFab primInst, primMeans, primVars;
    MultiFab coVars;

    // For long-range temperature-related correlations
    MultiFab cvlMeans, cvlInst, QMeans;

    int ncross = 38+nspecies*nspecies + nspecies;
    MultiFab spatialCross1D;

    int iprim = 10;
    int nprim = (nspecies+1)*iprim;


    if (seed > 0)
    {
        InitRandom(seed+ParallelDescriptor::MyProc(),
                   ParallelDescriptor::NProcs(),
                   seed+ParallelDescriptor::MyProc());
    }
    else if (seed == 0)
    {
        auto now = time_point_cast<nanoseconds>(system_clock::now());
        int randSeed = now.time_since_epoch().count();
        ParallelDescriptor::Bcast(&randSeed,1,ParallelDescriptor::IOProcessorNumber());
        InitRandom(randSeed+ParallelDescriptor::MyProc(),
                   ParallelDescriptor::NProcs(),
                   randSeed+ParallelDescriptor::MyProc());
    }
    else
    {
        Abort("Must supply non-negative seed");
    }


    if (restart < 0)
    {
        if (seed > 0)
        {
            InitRandom(seed+ParallelDescriptor::MyProc(),
                       ParallelDescriptor::NProcs(),
                       seed+ParallelDescriptor::MyProc());
        }
        else if (seed == 0)
        {
            auto now = time_point_cast<nanoseconds>(system_clock::now());
            int randSeed = now.time_since_epoch().count();
            // broadcast the same root seed to all processors
            ParallelDescriptor::Bcast(&randSeed,1,ParallelDescriptor::IOProcessorNumber());
            InitRandom(randSeed+ParallelDescriptor::MyProc(),
                       ParallelDescriptor::NProcs(),
                       randSeed+ParallelDescriptor::MyProc());
        }
        else
        {
            Abort("Must supply non-negative seed");
        }

        ba.define(domain);
        ba.maxSize(IntVect(max_grid_size));
        dmap.define(ba);
        //////////////////////////////////////
        // Conserved/Primitive Var Setup
        //////////////////////////////////////
        /*
            Conserved Vars:
            0  - rho = (1/V) += m
            1  - Jx  = (1/V) += mu
            2  - Jy  = (1/V) += mv
            3  - Jz  = (1/V) += mw
            4  - K   = (1/V) += m|v|^2
            ... (repeat for each species)
        */

        int ncon  = (nspecies+1)*5;
        cuInst.define(ba, dmap, ncon, 0); cuInst.setVal(0.);
        cuMeans.define(ba, dmap, ncon, 0); cuMeans.setVal(0.);
        cuVars.define(ba,dmap, ncon, 0); cuVars.setVal(0.);

        /*
           Primitive Vars:
            0 - n   (n_ns)
            1 - Yk  (Y_ns)
            2 - u   (u_ns)
            3 - v   (v_ns)
            4 - w   (w_ns)
            5 - G   (G_ns) = dot(u_mean,dJ)
            6 - T   (T_ns)
            7 - P   (P_ns)
            8 - E   (E_ns)
            ... (repeat for each species)
        */


        primInst.define(ba, dmap, nprim, 0); primInst.setVal(0.);
        primMeans.define(ba, dmap, nprim, 0); primMeans.setVal(0.);
        primVars.define(ba, dmap, ncon+nprim, 0); primVars.setVal(0.);

        // Covariances
        /*
            // Conserved
            0  - drho.dJx
            1  - drho.dJy
            2  - drho.dJz
            3  - drho.dK
            4  - dJx.dJy
            5  - dJx.dJz
            6  - dJx.dK
            7  - dJy.dJz
            8  - dJy.dK
            9  - dJz.dk

            // Energy
            10 - drho.dG
            11 - dJx.dG
            12 - dJy.dG
            13 - dJz.dG
            14 - dK.dG

            // Hydro
            15 - drho.du
            16 - drho.dv
            17 - drho.dw
            18 - du.dv
            19 - du.dw
            20 - dv.dw
            21 - drho.dT
            22 - du.dT
            23 - dv.dT
            24 - dw.dT
        */

        int ncovar = 25;
        coVars.define(ba, dmap, ncovar, 0); coVars.setVal(0.);


        spatialCross1D.define(ba,dmap,ncross,0); spatialCross1D.setVal(0.);

    }
    else
    {
        ReadCheckPoint(step, time, dt, statsCount,
            cuInst, cuMeans, cuVars,
            primInst, primMeans, primVars,
            coVars, spatialCross1D, ncross);
        dmap = cuInst.DistributionMap();
        ba = cuInst.boxArray();

    }

    // Specific Heat
    int ncvl = nspecies+1;
    cvlInst.define(ba, dmap, ncvl, 0);  cvlInst.setVal(0.);
    cvlMeans.define(ba, dmap, ncvl, 0); cvlMeans.setVal(0.);
    QMeans.define(ba, dmap, ncvl, 0); QMeans.setVal(0.);

    Vector<int> is_periodic (AMREX_SPACEDIM,0);
    for (int i=0; i<AMREX_SPACEDIM; ++i)
    {
        if (bc_vel_lo[i] == -1 && bc_vel_hi[i] == -1)
        {
            is_periodic [i] = -1;
        }
    }

    // This defines a Geometry object
    RealBox realDomain({AMREX_D_DECL(prob_lo[0],prob_lo[1],prob_lo[2])},
        {AMREX_D_DECL(prob_hi[0],prob_hi[1],prob_hi[2])});

    Geometry geom (domain ,&realDomain,CoordSys::cartesian,is_periodic.data());

    std::ifstream planeFile("paramplanes.dat");
    int fileCount = 0;
    if(planeFile.good())
    {
        planeFile >> fileCount;
    }

    int paramPlaneCount = 6 + fileCount;
    paramPlane paramPlaneList[paramPlaneCount];
    BuildParamplanes(paramPlaneList,paramPlaneCount,realDomain.lo(),realDomain.hi());

    // Particle tile size
    Vector<int> ts(BL_SPACEDIM);

    for (int d=0; d<AMREX_SPACEDIM; ++d)
    {
        if (max_particle_tile_size[d] > 0)
        {
            ts[d] = max_particle_tile_size[d];
        }
        else
        {
            ts[d] = max_grid_size[d];
        }
    }

    ParmParse pp ("particles");
    pp.addarr("tile_size", ts);

    int cRange = 0;
    FhdParticleContainer particles(geom, dmap, ba, cRange);

    // Output all primitives for structure factor
    int nvarstruct = 6+nspecies*2;
    const Real* dx = geom.CellSize();
    int nstruct = std::ceil((double)nvarstruct*(nvarstruct+1)/2);
    // scale SF results by inverse cell volume
    Vector<Real> var_scaling(nstruct);
    for (int d=0; d<var_scaling.size(); ++d)
    {
        var_scaling[d] = 1./(dx[0]*dx[1]*dx[2]);
    }

    // Standard 3D structure factors

    MultiFab structFactPrimMF;
    int structVarsPrim = (nspecies+1)*4;

    Vector< std::string > primNames;
    primNames.resize(structVarsPrim);
    std::string x;

    int cnt=0;

    primNames[cnt++] = "rhoInstant";
    primNames[cnt++] = "uInstant";
    primNames[cnt++] = "vInstant";
    primNames[cnt++] = "cInstant";

    for(int ispec=0;ispec<nspecies;ispec++) {

        primNames[cnt++] = amrex::Concatenate("rhoInstant_",ispec,2);
        primNames[cnt++] = amrex::Concatenate("uInstant_",ispec,2);
        primNames[cnt++] = amrex::Concatenate("vInstant_",ispec,2);
        primNames[cnt++] = amrex::Concatenate("cInstant_",ispec,2);
    }

    Vector<Real> var_scaling_prim;
    var_scaling_prim.resize(structVarsPrim*(structVarsPrim+1)/2);
    for (int d=0; d<var_scaling_prim.size(); ++d) {
        var_scaling_prim[d] = 1./(dx[0]*dx[1]*dx[2]);
    }


    structFactPrimMF.define(ba, dmap, structVarsPrim, 0);
    StructFact structFactPrim(ba,dmap,primNames,var_scaling_prim);



    // Collision Cell Vars
    particles.mfselect.define(ba, dmap, nspecies*nspecies, 0);
    particles.mfselect.setVal(0.);
    particles.mfphi.define(ba, dmap, nspecies, 0);
    particles.mfphi.setVal(0.);
    particles.mfvrmax.define(ba, dmap, nspecies*nspecies, 0);
    particles.mfvrmax.setVal(0.);

    Print() << "init particles\n";
    particles.InitParticles(dt);
    Print() << "init cells\n";
    particles.InitCollisionCells();

    Real init_time = ParallelDescriptor::second() - strt_time;
    ParallelDescriptor::ReduceRealMax(init_time);
    amrex::Print() << "Initialization time = " << init_time << " seconds " << std::endl;

    max_step += step;
    n_steps_skip += step;
    Real tbegin, tend;


    particles.zeroCells();
    zeroMassFlux(paramPlaneList, paramPlaneCount);

    //Initial condition
    spatialCross1D.setVal(0.);
    cuMeans.setVal(0.);
    primMeans.setVal(0.);
    cuVars.setVal(0.);
    primVars.setVal(0.);
    coVars.setVal(0.);

    //particles.SortParticlesDB();
    particles.EvaluateStats(cuInst,cuMeans,cuVars,primInst,primMeans,primVars,
        cvlInst,cvlMeans,QMeans,coVars,spatialCross1D,statsCount++,time);
    int stepTemp = 0;

    writePlotFile(cuInst,cuMeans,cuVars,primInst,primMeans,primVars,
            coVars,spatialCross1D,particles,geom,time,ncross,stepTemp);



    for (int istep=step; istep<=max_step; ++istep)
    {
        tbegin = ParallelDescriptor::second();

        particles.CalcSelections(dt);
        particles.CollideParticles(dt);

//      if(istep%2!=0)
        if(false)
        {
            particles.EvaluateStats(cuInst,cuMeans,cuVars,primInst,primMeans,primVars,
                    cvlInst,cvlMeans,QMeans,coVars,spatialCross1D,statsCount++,time);

            if (istep > n_steps_skip && struct_fact_int > 0 && (istep-n_steps_skip)%struct_fact_int == 0) {


                for(int l=0;l<=nspecies;l++)
                {
                    MultiFab::Copy(structFactPrimMF, primInst, l*iprim+1, l*4+0, 1, 0);
                    MultiFab::Copy(structFactPrimMF, primInst, l*iprim+2, l*4+1, 1, 0);
                    MultiFab::Copy(structFactPrimMF, primInst, l*iprim+3, l*4+2, 1, 0);
                    MultiFab::Copy(structFactPrimMF, primInst, l*iprim+9, l*4+3, 1, 0);

                    //Print() << l*iprim+1 << " -> " << l*4+0 << endl;
                }

                //PrintMF(structFactPrimMF,0,-1);
                //PrintMF(primInst,1,1);

                //structFactPrim.FortStructure(structFactPrimMF);

            }
        }

        //particles.externalForce(dt);
        particles.Source(dt, paramPlaneList, paramPlaneCount, cuInst);
        particles.MoveParticlesCPP(dt, paramPlaneList, paramPlaneCount);
        //particles.updateTimeStep(geom,dt);
                //reduceMassFlux(paramPlaneList, paramPlaneCount);

        if(true)
        //if(istep%2==0)
        {
            particles.EvaluateStats(cuInst,cuMeans,cuVars,primInst,primMeans,primVars,
                    cvlInst,cvlMeans,QMeans,coVars,spatialCross1D,statsCount++,time);

            if (istep > n_steps_skip && struct_fact_int > 0 && (istep-n_steps_skip)%struct_fact_int == 0) {


                for(int l=0;l<=nspecies;l++)
                {
                    MultiFab::Copy(structFactPrimMF, primInst, l*iprim+1, l*4+0, 1, 0);
                    MultiFab::Copy(structFactPrimMF, primInst, l*iprim+2, l*4+1, 1, 0);
                    MultiFab::Copy(structFactPrimMF, primInst, l*iprim+3, l*4+2, 1, 0);
                    MultiFab::Copy(structFactPrimMF, primInst, l*iprim+9, l*4+3, 1, 0);

                    //Print() << l*iprim+1 << " -> " << l*4+0 << endl;
                }

                //PrintMF(structFactPrimMF,0,-1);
                //PrintMF(primInst,1,1);

                //structFactPrim.FortStructure(structFactPrimMF);

            }
        }
        //////////////////////////////////////
        // PlotFile
        //////////////////////////////////////

        bool writePlt = false;
        if (plot_int > 0 && istep>0 && istep>=n_steps_skip)
        {
            if (n_steps_skip >= 0) // for positive n_steps_skip, write out at plot_int
            {
                writePlt = (istep%plot_int == 0);
            }
            else if (n_steps_skip < 0) // for negative n_steps_skip, write out at plot_int-1
            {
                writePlt = ((istep+1)%plot_int == 0);
            }
        }

        if ((plot_int > 0 && istep%plot_int==0))
        {
            writePlotFile(cuInst,cuMeans,cuVars,primInst,primMeans,primVars,
                coVars,spatialCross1D,particles,geom,time,ncross,istep);

            structFactPrim.WritePlotFile(istep,time,"plt_SF_prim");
        }

        if ((n_steps_skip > 0 && istep == n_steps_skip) || (n_steps_skip < 0 && istep%n_steps_skip == 0) ) {
            //reset stats
            statsCount = 1;
            spatialCross1D.setVal(0.);
            cuMeans.setVal(0.);
            primMeans.setVal(0.);
            cuVars.setVal(0.);
            primVars.setVal(0.);
            coVars.setVal(0.);

        }

        if (chk_int > 0 && istep%chk_int == 0 && istep > step)
        {
            WriteCheckPoint(istep, time, dt, statsCount,
                cuInst, cuMeans, cuVars, primInst, primMeans, primVars, coVars,
                particles, spatialCross1D, ncross);
        }
        tend = ParallelDescriptor::second() - tbegin;
        ParallelDescriptor::ReduceRealMax(tend);
        if(istep%100==0)
        {
            amrex::Print() << "Advanced step " << istep << " of " << max_step << " in " << tend << " seconds. " << particles.simParticles << " particles.\n";
        }

        time += dt;
    }

    Real stop_time = ParallelDescriptor::second() - strt_time;
    ParallelDescriptor::ReduceRealMax(stop_time);
    amrex::Print() << "Run time = " << stop_time << " seconds" << std::endl;
}
