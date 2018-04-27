namespace common {

    // Problem specification
    extern amrex::Vector<amrex::Real> prob_lo;       // physical lo coordinate
    extern amrex::Vector<amrex::Real> prob_hi;       // physical hi coordinate
    extern amrex::Vector<int>         n_cells;       // number of cells in domain
    extern amrex::Vector<int>         max_grid_size; // max number of cells in a box

    // Time-step control
    extern amrex::Real                fixed_dt;
    extern amrex::Real                cfl;

    // Controls for number of steps between actions
    extern int                        max_step;
    extern int                        plot_int;
    extern std::string                plot_base_name;
    extern int                        chk_int;
    extern std::string                chk_base_name;
    extern int                        prob_type;
    extern int                        restart;
    extern int                        print_int;
    extern int                        project_eos_int;

    // Physical parameters
    extern amrex::Vector<amrex::Real> grav;
    extern int                        nspecies;
    extern amrex::Vector<amrex::Real> molmass;
    extern amrex::Vector<amrex::Real> rhobar;
    extern amrex::Real                rho0;

    // stochastic forcing amplitudes (1 for physical values, 0 to run them off)
    extern amrex::Real                variance_coef_mom;
    extern amrex::Real                variance_coef_mass;
    extern amrex::Real                k_B;
    extern amrex::Real                Runiv;

    // Algorithm control / selection
    extern int                        algorithm_type;
    extern int                        barodiffusion_type;
    extern bool                       use_bl_rng;

    // random number seed (for HydroGrid RNGs)
    // 0        = unpredictable seed based on clock
    // positive = fixed seed
    extern int                        seed;

    // Random number seeds for each physical process for use_bl_rng=T
    // for positive value, the value is assigned as seed value
    // for 0, a positive value is randomly chosen
    // if -1 (only for restart), RNGs status is restored from checkpoint data
    extern int                        seed_momentum;
    extern int                        seed_diffusion;
    extern int                        seed_reaction;
    extern int                        seed_init_mass;
    extern int                        seed_init_momentum;

    // Viscous friction L phi operator
    // if abs(visc_type) = 1, L = div beta grad
    // if abs(visc_type) = 2, L = div [ beta (grad + grad^T) ]
    // if abs(visc_type) = 3, L = div [ beta (grad + grad^T) + I (gamma - (2/3)*beta) div ]
    // positive = assume constant coefficients
    // negative = assume spatially-varying coefficients
    extern amrex::Real                visc_coef;
    extern int                        visc_type;

    extern int                        advection_type;

    // Stochastic momentum flux controls:
    extern int                        filtering_width;
    extern int                        stoch_stress_form;

    // Initial conditions
    extern amrex::Vector<amrex::Real> u_init;
    extern amrex::Real                perturb_width;
    extern amrex::Real                smoothing_width;
    extern amrex::Real                initial_variance_mom;
    extern amrex::Real                initial_variance_mass;

    // Boundary conditions
    // ----------------------
    // BC specifications:
    // -1 = periodic
    // 100 = no-slip wall      (Dir condition for normal vel; Dir velocity condition for trans vel)
    // 101 = no-slip reservoir (Dir condition for normal vel; Dir velocity condition for trans vel)
    // 200 = slip wall         (Dir condition for normal vel; Dir traction condition for trans vel)
    // 201 = slip reservoir    (Dir condition for normal vel; Dir traction condition for trans vel)
    // For a complete list see bc.f90
    extern amrex::Vector<int>         bc_lo;
    extern amrex::Vector<int>         bc_hi;

    // Each no-slip wall may be moving with a specified tangential 
    // velocity along the tangential directions
    // In 2D:
    // wallspeed_lo/hi(1,1) - yvel on x-face
    // wallspeed_lo/hi(1,2) - xvel on y-face
    // In 3D:
    // wallspeed_lo/hi(1,1) - yvel on x-face
    // wallspeed_lo/hi(2,1) - zvel on x-face
    // wallspeed_lo/hi(1,2) - xvel on y-face
    // wallspeed_lo/hi(2,2) - zvel on y-face
    // wallspeed_lo/hi(1,3) - xvel on z-face
    // wallspeed_lo/hi(2,3) - yvel on z-face
    extern amrex::Vector<amrex::Real> wallspeed_lo_x;
    extern amrex::Vector<amrex::Real> wallspeed_hi_x;
    extern amrex::Vector<amrex::Real> wallspeed_lo_y;
    extern amrex::Vector<amrex::Real> wallspeed_hi_y;
    extern amrex::Vector<amrex::Real> wallspeed_lo_z;
    extern amrex::Vector<amrex::Real> wallspeed_hi_z;


    // Hydrogrid analysis (not available yet)
/*
  extern int                        hydro_grid_int;
  extern int                        project_dir;
  extern amrex::Vector<int>         max_grid_projection;
  extern int                        stats_int;
  extern int                        n_steps_save_stats;
  extern int                        n_steps_skip;
  extern bool                       analyze_conserved;
  extern bool                       center_snapshots;
*/

    // These are mostly used for reaction-diffusion: 
    extern int                        histogram_unit;
    extern amrex::Vector<amrex::Real> density_weights;

    extern amrex::Vector<int>         shift_cc_to_boundary_lo;
    extern amrex::Vector<int>         shift_cc_to_boundary_hi;
}

