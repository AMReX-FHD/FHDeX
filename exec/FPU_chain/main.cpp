#include <AMReX.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParmParse.H>

#include "FPU.H"

// for clock-based random seed
#include "chrono"
#include <array>
#include <vector>
using namespace std::chrono;

using namespace amrex;

int main (int argc, char* argv[])
{
Initialize(argc,argv);
{

    if (AMREX_SPACEDIM != 2) {
        Abort("Must build with DIM=2");
    }

    // **********************************
    // DECLARE SIMULATION PARAMETERS
    // **********************************

    int n_particles = 8000;
    int n_steps = 1000000;

    int n_ensembles = 10;
    int max_ensembles_per_rank = 1000;  // for parallelization purposes

    Real a_coef = 1.0;
    Real b_coef = 0.0;
    Real c_coef = 1.0;
    Real rest   = 1.0;

    Real beta = 1.0;
    Real pressure = 1.0;

    Real r_eq = 1.;
    Real p_eq = 0.;
    Real e_eq = 1.;

    Real R_00 = 1.;
    Real R_01 = 0.;
    Real R_02 = 0.;
    Real R_10 = 0.;
    Real R_11 = 1.;
    Real R_12 = 0.;
    Real R_20 = 0.;
    Real R_21 = 0.;
    Real R_22 = 1.;

    int seed = 1;

    Real dt = 1.e-3;

    int diag_int = 1000;
    int plot_int = 10000;

    int hf_nevery  = 10;
    int hf_nrepeat = 100;
    int hf_nfreq   = 1000;

    // ***************************************************************
    // READ INPUT PARAMETER VALUES FROM INPUT FILE AND/OR COMMAND LINE
    // ***************************************************************
    {
        // ParmParse is way of reading inputs from the inputs file
        // pp.get means we require the inputs file to have it
        // pp.query means we optionally need the inputs file to have it - but we must supply a default here
        ParmParse pp;

        pp.query("n_particles",n_particles);
        pp.query("n_steps",n_steps);

        pp.query("n_ensembles",n_ensembles);
        pp.query("max_ensembles_per_rank",max_ensembles_per_rank);

        pp.query("a_coef",a_coef);
        pp.query("b_coef",b_coef);
        pp.query("c_coef",c_coef);
        pp.query("rest",rest);

        pp.query("beta",beta);
        pp.query("pressure",pressure);

        pp.query("r_eq",r_eq);
        pp.query("p_eq",p_eq);
        pp.query("e_eq",e_eq);

        pp.query("R_00",R_00);
        pp.query("R_01",R_01);
        pp.query("R_02",R_02);
        pp.query("R_10",R_10);
        pp.query("R_11",R_11);
        pp.query("R_12",R_12);
        pp.query("R_20",R_20);
        pp.query("R_21",R_21);
        pp.query("R_22",R_22);

        pp.query("seed",seed);

        pp.query("dt",dt);

        pp.query("diag_int",diag_int);
        pp.query("plot_int",plot_int);

        pp.query("hf_nevery",hf_nevery);
        pp.query("hf_nrepeat",hf_nrepeat);
        pp.query("hf_nfreq",hf_nfreq);
    }

    if (hf_nevery <= 0 || hf_nrepeat <= 0 || hf_nfreq <= 0) {
        Abort("hf_nevery, hf_nrepeat, and hf_nfreq must be positive");
    }

    if (hf_nfreq % hf_nevery != 0) {
        Abort("hf_nfreq must be an integer multiple of hf_nevery");
    }

    const long long nsamples_first_output = static_cast<long long>(hf_nfreq) / hf_nevery + 1;
    if (hf_nrepeat > nsamples_first_output) {
        Abort("hf_nrepeat is too large to be fully populated by the first output time");
    }

    int the_seed;
    if (seed > 0) {
        // fixed seed
        the_seed = seed;
    } else if (seed == 0) {
        // initializes the seed for C++ random number calls based on the clock
        // broadcast the same seed to all processors
        auto now = time_point_cast<nanoseconds>(system_clock::now());
        the_seed = now.time_since_epoch().count();
        ParallelDescriptor::Bcast(&the_seed,1,ParallelDescriptor::IOProcessorNumber());
    } else {
        Abort("Must supply non-negative seed");
    }
    InitRandom(the_seed+ParallelDescriptor::MyProc(),
               ParallelDescriptor::NProcs(),
               seed+ParallelDescriptor::MyProc());

    // define lower and upper indices
    IntVect dom_lo(0,0);
    IntVect dom_hi(n_particles-1, n_ensembles-1);

    // Make a single box that is the entire domain
    Box domain(dom_lo, dom_hi);

    // Initialize the boxarray "ba" from the single box "domain"
    BoxArray ba(domain);

    IntVect max_grid_size(1024000,max_ensembles_per_rank);

    // Break up boxarray "ba" into chunks no larger than "max_grid_size" along a direction
    ba.maxSize(max_grid_size);

    // How Boxes are distributed among MPI processes
    DistributionMapping dm(ba);

    // physical box size in this FPU_chain implementation is not relevant
    RealBox real_box({ 0., 0.,}, { 1., 1.,});

    // periodic in x only
    Array<int,AMREX_SPACEDIM> is_periodic{1,0};

    // This defines a Geometry object
    Geometry geom(domain, real_box, CoordSys::cartesian, is_periodic);

    // time = starting time in the simulation
    Real time = 0.0;

    // we only need a ghost cell for the x-direction
    IntVect ng_vect(1,0);

    // components are r, p, and e
    MultiFab state(ba,dm,3,ng_vect);

    // heat flux: 2 components: 1. convective; 2. conductive; 3. total
    MultiFab heat_flux(ba,dm,3,ng_vect);

    // storage for phi = R*state
    MultiFab phi (ba,dm,3,0);
    MultiFab phi0(ba,dm,3,0);

    // diagonal elements of C_alphaalpha
    MultiFab C_alphaalpha(ba,dm,3,0);

    Gpu::HostVector<Real> C_alphaalpha_00(n_particles);
    Gpu::HostVector<Real> C_alphaalpha_11(n_particles);
    Gpu::HostVector<Real> C_alphaalpha_22(n_particles);

    // Running average of the pointwise heat-flux autocorrelation, 
    // sampled every hf_nevery steps and averaged over the domain.
    // hf_nevery steps and averaged in the LAMMPS "ave running" sense.
    Vector<std::unique_ptr<MultiFab>> heat_flux_history(hf_nrepeat);
    for (int i = 0; i < hf_nrepeat; ++i) {
        heat_flux_history[i] = std::make_unique<MultiFab>(ba,dm,1,0);
        heat_flux_history[i]->setVal(0.0);
    }
    MultiFab heat_flux_corr(ba,dm,1,0);
    std::vector<Real> hf_corr_sum(hf_nrepeat, 0.0_rt);
    std::vector<Real> hf_corr_running(hf_nrepeat, 0.0_rt);
    std::vector<long long> hf_corr_count(hf_nrepeat, 0);

    int hf_lastindex = -1;
    int hf_history_size = 0;

    const Real inv_npts = 1.0_rt / static_cast<Real>(n_particles*n_ensembles);

    auto sample_heat_flux_corr = [&](int step)
    {
        if (step % hf_nevery != 0) {
            return;
        }

        hf_lastindex++;
        if (hf_lastindex == hf_nrepeat) {
            hf_lastindex = 0;
        }

        MultiFab::Copy(*heat_flux_history[hf_lastindex], heat_flux, 2, 0, 1, 0);

        if (hf_history_size < hf_nrepeat) {
            ++hf_history_size;
        }

        int history_index = hf_lastindex;
        for (int lag = 0; lag < hf_history_size; ++lag) {
            const MultiFab& heat_flux_old = *heat_flux_history[history_index];

            for (MFIter mfi(heat_flux); mfi.isValid(); ++mfi) {
                const Box& bx = mfi.validbox();

                const Array4<Real const>& heat_flux_now = heat_flux.const_array(mfi);
                const Array4<Real const>& heat_flux_prev = heat_flux_old.const_array(mfi);
                const Array4<Real>& heat_flux_corr_fab = heat_flux_corr.array(mfi);

                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    heat_flux_corr_fab(i,j,k) = heat_flux_now(i,j,k,2) * heat_flux_prev(i,j,k,0);
                });
            }

            hf_corr_sum[lag] += heat_flux_corr.sum(0)*inv_npts;
            ++hf_corr_count[lag];

            --history_index;
            if (history_index < 0) {
                history_index = hf_nrepeat - 1;
            }
        }
    };

    // ******************************
    // SAMPLE TO OBTAIN INITIAL STATE
    // ******************************
    init(state, beta, pressure, a_coef, b_coef, c_coef, 0., 10000, 1.e-3, n_particles, n_ensembles, geom);
    compute_energy(state,a_coef,b_coef,c_coef);
    compute_heat_flux(heat_flux,state,a_coef,b_coef,c_coef,rest,geom);
    sample_heat_flux_corr(0);

    Copy(phi0,state,0,0,3,0);
    ComputePhiFromState(phi0,r_eq,p_eq,e_eq,R_00,R_01,R_02,R_10,R_11,R_12,R_20,R_21,R_22);

/*
    // save g_alpha(0,0)
    g_alpha_zero.ParallelCopy(state, 0, 0, 3);
*/

    // write out diagnostics (means)
    if (diag_int > 0) {
        compute_means(state,n_particles,n_ensembles,0);
    }

    // write the initial state (r,p,e) to a plotfile
    if (plot_int > 0) {
        MultiFab plotdata(ba,dm,6,0);
        MultiFab::Copy(plotdata,state,0,0,3,0);
        MultiFab::Copy(plotdata,heat_flux,0,3,3,0);
        const std::string& pltfile = amrex::Concatenate("plt",0,7);
        amrex::Print() << "Writing plotfile " << pltfile << std::endl;
        WriteSingleLevelPlotfile(pltfile, plotdata, {"r","p","e","heat_flux_conv","heat_flux_cond","heat_flux_total"}, geom, time, 0);
    }

    for (int step=1; step<=n_steps; ++step) {

        Real step_strt_time = ParallelDescriptor::second();

        time += dt;

        // ****************
        // INTEGRATE A STEP
        // ****************
        FPU_RK4(state,a_coef,b_coef,c_coef,dt,n_particles,n_ensembles,geom);
        compute_energy(state,a_coef,b_coef,c_coef);
        compute_heat_flux(heat_flux,state,a_coef,b_coef,c_coef,rest,geom);
        sample_heat_flux_corr(step);

        if (step % hf_nfreq == 0) {
            for (int lag = 0; lag < hf_nrepeat; ++lag) {
                if (hf_corr_count[lag] > 0) {
                    hf_corr_running[lag] =
                        hf_corr_sum[lag] / static_cast<Real>(hf_corr_count[lag]);
                } else {
                    hf_corr_running[lag] = 0.0_rt;
                }
            }

            const std::string hf_corr_file = amrex::Concatenate("heat_flux_corr",step,7);
            amrex::Print() << "Writing heat-flux correlation file "
                           << hf_corr_file << std::endl;

            if (ParallelDescriptor::IOProcessor()) {
                std::ofstream hfout(hf_corr_file, std::ios::out);
                for (int lag = 0; lag < hf_nrepeat; ++lag) {
                    hfout << lag << " "
                          << lag*hf_nevery << " "
                          << hf_corr_count[lag] << " "
                          << hf_corr_running[lag] << "\n";
                }
            }
        }

        Real step_stop_time = ParallelDescriptor::second() - step_strt_time;
        ParallelDescriptor::ReduceRealMax(step_stop_time);

        amrex::Print() << "Completed step " << step << " in " << step_stop_time << " seconds " << std::endl;

        if (plot_int > 0 && step%plot_int == 0) {

            Real plot_strt_time = ParallelDescriptor::second();

            // write the current state (r,p,e,heat_flux_conv,heat_flux_cond) to a plotfile
            MultiFab plotdata(ba,dm,6,0);
            MultiFab::Copy(plotdata,state,0,0,3,0);
            MultiFab::Copy(plotdata,heat_flux,0,3,3,0);
            const std::string& pltfile = amrex::Concatenate("plt",step,7);
            amrex::Print() << "Writing plotfile " << pltfile << std::endl;
            WriteSingleLevelPlotfile(pltfile, plotdata, {"r","p","e","heat_flux_conv","heat_flux_cond","heat_flux_total"}, geom, time, step);

            Copy(phi,state,0,0,3,0);
            ComputePhiFromState(phi,r_eq,p_eq,e_eq,R_00,R_01,R_02,R_10,R_11,R_12,R_20,R_21,R_22);

            ComputeCalphaalpha(C_alphaalpha,phi,phi0);
            C_alphaalpha_00 = sumToLine(C_alphaalpha, 0, 1, domain, 0);
            C_alphaalpha_11 = sumToLine(C_alphaalpha, 1, 1, domain, 0);
            C_alphaalpha_22 = sumToLine(C_alphaalpha, 2, 1, domain, 0);

            const std::string C_alphaalphafile = amrex::Concatenate("C_alphaalpha",step,7);
            amrex::Print() << "Writing C_alphaalphafile " << C_alphaalphafile << std::endl;
            std::ofstream C_alphaalphaout;
            if (ParallelDescriptor::IOProcessor()) {
                C_alphaalphaout.open(C_alphaalphafile, std::ios::out);
                for (int i=0; i<n_particles; ++i) {

                    C_alphaalphaout << " C_alphaalpha_00/11/22 = " << i << " "
                                    << C_alphaalpha_00[(i+n_particles/2)%n_particles] << " "
                                    << C_alphaalpha_11[(i+n_particles/2)%n_particles] << " "
                                    << C_alphaalpha_22[(i+n_particles/2)%n_particles] << "\n";
                }
            }

            Real plot_stop_time = ParallelDescriptor::second() - plot_strt_time;
            ParallelDescriptor::ReduceRealMax(plot_stop_time);

            amrex::Print() << "Wrote out plot data for step " << step << " in " << plot_stop_time << " seconds\n";
        }

        // write out diagnostics (means)
        if (diag_int > 0 && step%diag_int == 0) {

            Real diag_strt_time = ParallelDescriptor::second();

            compute_means(state,n_particles,n_ensembles,step);

            Real diag_stop_time = ParallelDescriptor::second() - diag_strt_time;
            ParallelDescriptor::ReduceRealMax(diag_stop_time);

            amrex::Print() << "Wrote out diag data for step " << step << " in " << diag_stop_time << " seconds\n";
        }

    }



}
Finalize();
return 0;
}
