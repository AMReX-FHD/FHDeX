#include "FPU.H"

// number of cells in each spatial direction
int FPU::n_cell_x;
int FPU::n_ensembles;
int FPU::max_ensembles_per_rank; // for parallelization purposes

// total steps in simulation and time step
int FPU::nsteps;
amrex::Real FPU::dt;

// enable fluctuations
int FPU::enable_fluctuations;

// how many steps to skip before defining t=0
int FPU::n_steps_skip;

// how often to write a plotfile
int FPU::plot_int;

// how often to write out correlation diagnostics
int FPU::diag_int;

// random number seed (positive integer=fixed seed; 0=clock-based seed)
int FPU::seed;

// size of each finite volume cell - all 3 must be defined regardless of dimensionality
amrex::Real FPU::cell_dx;
amrex::Real FPU::cell_dy;
amrex::Real FPU::cell_dz;

AMREX_GPU_MANAGED amrex::Real FPU::r0;
AMREX_GPU_MANAGED amrex::Real FPU::p0;
AMREX_GPU_MANAGED amrex::Real FPU::e0;

AMREX_GPU_MANAGED amrex::Real FPU::A_00;
AMREX_GPU_MANAGED amrex::Real FPU::A_01;
AMREX_GPU_MANAGED amrex::Real FPU::A_02;
AMREX_GPU_MANAGED amrex::Real FPU::A_10;
AMREX_GPU_MANAGED amrex::Real FPU::A_11;
AMREX_GPU_MANAGED amrex::Real FPU::A_12;
AMREX_GPU_MANAGED amrex::Real FPU::A_20;
AMREX_GPU_MANAGED amrex::Real FPU::A_21;
AMREX_GPU_MANAGED amrex::Real FPU::A_22;

AMREX_GPU_MANAGED amrex::Real FPU::D_00;
AMREX_GPU_MANAGED amrex::Real FPU::D_01;
AMREX_GPU_MANAGED amrex::Real FPU::D_02;
AMREX_GPU_MANAGED amrex::Real FPU::D_10;
AMREX_GPU_MANAGED amrex::Real FPU::D_11;
AMREX_GPU_MANAGED amrex::Real FPU::D_12;
AMREX_GPU_MANAGED amrex::Real FPU::D_20;
AMREX_GPU_MANAGED amrex::Real FPU::D_21;
AMREX_GPU_MANAGED amrex::Real FPU::D_22;

AMREX_GPU_MANAGED amrex::Real FPU::B_00;
AMREX_GPU_MANAGED amrex::Real FPU::B_01;
AMREX_GPU_MANAGED amrex::Real FPU::B_02;
AMREX_GPU_MANAGED amrex::Real FPU::B_10;
AMREX_GPU_MANAGED amrex::Real FPU::B_11;
AMREX_GPU_MANAGED amrex::Real FPU::B_12;
AMREX_GPU_MANAGED amrex::Real FPU::B_20;
AMREX_GPU_MANAGED amrex::Real FPU::B_21;
AMREX_GPU_MANAGED amrex::Real FPU::B_22;

AMREX_GPU_MANAGED amrex::Real FPU::R_00;
AMREX_GPU_MANAGED amrex::Real FPU::R_01;
AMREX_GPU_MANAGED amrex::Real FPU::R_02;
AMREX_GPU_MANAGED amrex::Real FPU::R_10;
AMREX_GPU_MANAGED amrex::Real FPU::R_11;
AMREX_GPU_MANAGED amrex::Real FPU::R_12;
AMREX_GPU_MANAGED amrex::Real FPU::R_20;
AMREX_GPU_MANAGED amrex::Real FPU::R_21;
AMREX_GPU_MANAGED amrex::Real FPU::R_22;

// dynamic structure factor control
int FPU::DFS_n_steps_skip;       // how many steps before we start recording (t0)
int FPU::DFS_stats_int;          // how often to take a snapshop
int FPU::DFS_num_snapshots;      // how many total snapshots
int FPU::DFS_num_padding;        // how much zero padding in time

void ComputePhiFromState(MultiFab& phi) {

    // convert state to state-state_eq
    phi.plus(-r0,0,1,0);
    phi.plus(-p0,1,1,0);
    phi.plus(-e0,2,1,0);

    // phi = R(state-state_eq)
    for (MFIter mfi(phi); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.tilebox();

        const Array4<Real>& phi_fab = phi.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real R0 = R_00*phi_fab(i,j,k,0) + R_01*phi_fab(i,j,k,1) + R_02*phi_fab(i,j,k,2);
            Real R1 = R_10*phi_fab(i,j,k,0) + R_11*phi_fab(i,j,k,1) + R_12*phi_fab(i,j,k,2);
            Real R2 = R_20*phi_fab(i,j,k,0) + R_21*phi_fab(i,j,k,1) + R_22*phi_fab(i,j,k,2);

            phi_fab(i,j,k,0) = R0;
            phi_fab(i,j,k,1) = R1;
            phi_fab(i,j,k,2) = R2;

        });
    }

}

void ComputeCalphaalpha(MultiFab& C_alphaalpha,
                        const MultiFab& phi,
                        const MultiFab& phi0) {

    for (MFIter mfi(C_alphaalpha); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.tilebox();

        const Array4<      Real>& C_alphaalpha_fab = C_alphaalpha.array(mfi);
        const Array4<const Real>& phi_fab          = phi.array(mfi);
        const Array4<const Real>& phi0_fab         = phi0.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            C_alphaalpha_fab(i,j,k,0) = phi_fab(i,j,k,0) * phi0_fab(0,j,k,0);
            C_alphaalpha_fab(i,j,k,1) = phi_fab(i,j,k,1) * phi0_fab(0,j,k,1);
            C_alphaalpha_fab(i,j,k,2) = phi_fab(i,j,k,2) * phi0_fab(0,j,k,2);
        });
    }

}


