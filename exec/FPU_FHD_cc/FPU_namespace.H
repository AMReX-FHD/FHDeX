namespace FPU {

    //////////////////////////////////////////////////////
    // See FPU.H for a description of these parameters
    //////////////////////////////////////////////////////

    extern int n_cell_x;
    extern int n_ensembles;
    extern int max_ensembles_per_rank;

    extern int nsteps;
    extern AMREX_GPU_MANAGED amrex::Real dt;

    extern AMREX_GPU_MANAGED int enable_fluctuations;

    extern int n_steps_skip;

    extern int plot_int;

    extern int diag_int;

    extern int seed;

    extern amrex::Real cell_dx;
    extern amrex::Real cell_dy;
    extern amrex::Real cell_dz;

    extern AMREX_GPU_MANAGED amrex::Real r0;
    extern AMREX_GPU_MANAGED amrex::Real p0;
    extern AMREX_GPU_MANAGED amrex::Real e0;

    extern AMREX_GPU_MANAGED amrex::Real A_00;
    extern AMREX_GPU_MANAGED amrex::Real A_01;
    extern AMREX_GPU_MANAGED amrex::Real A_02;
    extern AMREX_GPU_MANAGED amrex::Real A_10;
    extern AMREX_GPU_MANAGED amrex::Real A_11;
    extern AMREX_GPU_MANAGED amrex::Real A_12;
    extern AMREX_GPU_MANAGED amrex::Real A_20;
    extern AMREX_GPU_MANAGED amrex::Real A_21;
    extern AMREX_GPU_MANAGED amrex::Real A_22;

    extern AMREX_GPU_MANAGED amrex::Real D_00;
    extern AMREX_GPU_MANAGED amrex::Real D_01;
    extern AMREX_GPU_MANAGED amrex::Real D_02;
    extern AMREX_GPU_MANAGED amrex::Real D_10;
    extern AMREX_GPU_MANAGED amrex::Real D_11;
    extern AMREX_GPU_MANAGED amrex::Real D_12;
    extern AMREX_GPU_MANAGED amrex::Real D_20;
    extern AMREX_GPU_MANAGED amrex::Real D_21;
    extern AMREX_GPU_MANAGED amrex::Real D_22;

    extern AMREX_GPU_MANAGED amrex::Real B_00;
    extern AMREX_GPU_MANAGED amrex::Real B_01;
    extern AMREX_GPU_MANAGED amrex::Real B_02;
    extern AMREX_GPU_MANAGED amrex::Real B_10;
    extern AMREX_GPU_MANAGED amrex::Real B_11;
    extern AMREX_GPU_MANAGED amrex::Real B_12;
    extern AMREX_GPU_MANAGED amrex::Real B_20;
    extern AMREX_GPU_MANAGED amrex::Real B_21;
    extern AMREX_GPU_MANAGED amrex::Real B_22;

    extern AMREX_GPU_MANAGED amrex::Real R_00;
    extern AMREX_GPU_MANAGED amrex::Real R_01;
    extern AMREX_GPU_MANAGED amrex::Real R_02;
    extern AMREX_GPU_MANAGED amrex::Real R_10;
    extern AMREX_GPU_MANAGED amrex::Real R_11;
    extern AMREX_GPU_MANAGED amrex::Real R_12;
    extern AMREX_GPU_MANAGED amrex::Real R_20;
    extern AMREX_GPU_MANAGED amrex::Real R_21;
    extern AMREX_GPU_MANAGED amrex::Real R_22;

    extern int DFS_n_steps_skip;
    extern int DFS_stats_int;
    extern int DFS_num_snapshots;
    extern int DFS_num_padding;
}
