namespace common {

//////////////
#define ADJ 0.99999
#define ADJALT 0.00002
#define WRITE_BUFFER 1000
#define MAX_BONDS 12
//////////////

    // misc global variables
    extern amrex::IntVect nodal_flag;
    extern amrex::Vector<amrex::IntVect> nodal_flag_dir;
    extern amrex::Vector<amrex::IntVect> nodal_flag_edge;
    extern amrex::IntVect nodal_flag_x;
    extern amrex::IntVect nodal_flag_y;
    extern amrex::IntVect nodal_flag_z;
    extern amrex::IntVect nodal_flag_xy;
    extern amrex::IntVect nodal_flag_xz;
    extern amrex::IntVect nodal_flag_yz;

    // Below are items that are read in from the common namelist

    // Problem specification
    extern AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, 3> prob_lo; // physical lo coordinate
    extern AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, 3> prob_hi; // physical hi coordinate
    extern AMREX_GPU_MANAGED amrex::GpuArray<int, AMREX_SPACEDIM> n_cells;              // number of cells in domain
    extern amrex::Vector<int>         max_grid_size;            // max number of cells in a box
    extern amrex::Vector<int>         max_particle_tile_size;        // max number of cells in a box
    extern AMREX_GPU_MANAGED amrex::Real cell_depth;

    extern amrex::IntVect             ngc;           //number of chost cells
    extern AMREX_GPU_MANAGED int      nvars;         //number of conserved variables
    extern AMREX_GPU_MANAGED int      nprimvars;     //number of primitive variables

    extern AMREX_GPU_MANAGED int      cross_cell;    //cell to compute spatial correlation
    extern AMREX_GPU_MANAGED int      do_slab_sf;    //whether to compute SF for two slabs separated by cross_cell

    extern amrex::Vector<amrex::Real> qval;
    extern AMREX_GPU_MANAGED amrex::GpuArray<int, MAX_SPECIES> pkernel_fluid; // GALEN - FLUID KERNEL
    extern AMREX_GPU_MANAGED amrex::GpuArray<int, MAX_SPECIES> pkernel_es;
    extern AMREX_GPU_MANAGED amrex::GpuArray<int, MAX_SPECIES> eskernel_fluid;
    extern AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> eskernel_beta;

    // Time-step control
    extern amrex::Real                fixed_dt;
    extern amrex::Real                cfl;

    // Random finite difference size, fraction of cell
    extern amrex::Real                rfd_delta;

    // Controls for number of steps between actions
    extern int                        max_step;
    extern int                        ramp_step;
    extern int                        plot_int;
    extern int                        plot_stag;
    extern std::string                plot_base_name;
    extern int                        chk_int;
    extern std::string                chk_base_name;
    extern std::string                plot_init_file;
    extern AMREX_GPU_MANAGED int      prob_type;
    extern int                        restart;
    extern int                        stats_int;
    extern int                        reset_stats;
    extern int                        particle_restart;
    extern int                        print_int;
    extern int                        project_eos_int;

    // Physical parameters
    extern AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> grav;
    extern AMREX_GPU_MANAGED int nspecies;
    extern AMREX_GPU_MANAGED int nbonds;
    extern AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> molmass;
    extern AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> rhobar;
    extern AMREX_GPU_MANAGED amrex::Real rho0;
    extern AMREX_GPU_MANAGED amrex::Real mach0; // initial mach number
    extern AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> diameter;
    extern AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> dof;
    extern AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> e0;
    extern AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> hcv;
    extern AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> hcp;

    // stochastic forcing amplitudes (1 for physical values, 0 to run them off)
    extern AMREX_GPU_MANAGED amrex::Real variance_coef_mom;
    extern AMREX_GPU_MANAGED amrex::Real variance_coef_mass;
    extern AMREX_GPU_MANAGED amrex::Real variance_coef_ener;
    extern AMREX_GPU_MANAGED amrex::Real k_B;
    extern AMREX_GPU_MANAGED amrex::Real h_bar;
    extern AMREX_GPU_MANAGED amrex::Real Runiv;
    extern AMREX_GPU_MANAGED amrex::Real avogadro;
    extern AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> T_init;

    // Algorithm control / selection
    extern AMREX_GPU_MANAGED int      algorithm_type;
    // 0 = centered
    // 1 = unlimited bilinear bds
    // 2 = limited bilinear bds
    extern AMREX_GPU_MANAGED int      advection_type;
    extern int                        barodiffusion_type;

    // random number seed
    // 0        = unpredictable seed based on clock
    // positive = fixed seed
    extern int                        seed;

    // Viscous friction L phi operator
    // if abs(visc_type) = 1, L = div beta grad
    // if abs(visc_type) = 2, L = div [ beta (grad + grad^T) ]
    // if abs(visc_type) = 3, L = div [ beta (grad + grad^T) + I (gamma - (2/3)*beta) div ]
    // positive = assume constant coefficients
    // negative = assume spatially-varying coefficients
    extern AMREX_GPU_MANAGED amrex::Real visc_coef;
    extern AMREX_GPU_MANAGED int      visc_type;

    // Stochastic momentum flux controls:
    extern int                        filtering_width;
    extern int                        stoch_stress_form;

    // Initial conditions
    extern amrex::Vector<amrex::Real> u_init;
    extern amrex::Real                perturb_width;
    extern AMREX_GPU_MANAGED amrex::Real smoothing_width;
    extern AMREX_GPU_MANAGED amrex::Real radius_cyl;
    extern AMREX_GPU_MANAGED amrex::Real radius_outer;
    extern AMREX_GPU_MANAGED amrex::Real film_thickness;
    extern amrex::Real                initial_variance_mom;
    extern amrex::Real                initial_variance_mass;
    extern amrex::Real                domega;

    // Boundary conditions
    // ----------------------
    // BC specifications:
    // -1 = periodic
    //  1 = slip
    //  2 = no-slip
    extern AMREX_GPU_MANAGED amrex::GpuArray<int, AMREX_SPACEDIM>          bc_vel_lo;
    extern AMREX_GPU_MANAGED amrex::GpuArray<int, AMREX_SPACEDIM>          bc_vel_hi;

    // BC specifications:
    // -1 = periodic
    //  1 = Dirichlet
    //  2 = Neumann
    extern AMREX_GPU_MANAGED amrex::GpuArray<int, AMREX_SPACEDIM>          bc_es_lo;
    extern AMREX_GPU_MANAGED amrex::GpuArray<int, AMREX_SPACEDIM>          bc_es_hi;

    // BC specifications:
    // -1 = periodic
    //  1 = wall
    //  2 = concentration (for compressible codes), reservoir (for low Mach multispecies codes)
    //  3 = reservoir (for compressible codes)
    extern AMREX_GPU_MANAGED amrex::GpuArray<int, AMREX_SPACEDIM>          bc_mass_lo;
    extern AMREX_GPU_MANAGED amrex::GpuArray<int, AMREX_SPACEDIM>          bc_mass_hi;

    // BC specifications:
    // -1 = periodic
    //  1 = adiabatic
    //  2 = isothermal
    extern AMREX_GPU_MANAGED amrex::GpuArray<int, AMREX_SPACEDIM>          bc_therm_lo;
    extern AMREX_GPU_MANAGED amrex::GpuArray<int, AMREX_SPACEDIM>          bc_therm_hi;

    // Pressure drop are periodic inflow/outflow walls (bc_[hi,lo]=-2).
    extern AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> p_lo;
    extern AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> p_hi;

    extern AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> t_lo;
    extern AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> t_hi;

    extern AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> rho_lo;
    extern AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> rho_hi;

    extern AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> bc_Yk_x_lo;
    extern AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> bc_Yk_x_hi;
    extern AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> bc_Yk_y_lo;
    extern AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> bc_Yk_y_hi;
    extern AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> bc_Yk_z_lo;
    extern AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> bc_Yk_z_hi;

    extern AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> n_lo;
    extern AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> n_hi;

    extern AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> bc_Xk_x_lo;
    extern AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> bc_Xk_x_hi;
    extern AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> bc_Xk_y_lo;
    extern AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> bc_Xk_y_hi;
    extern AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> bc_Xk_z_lo;
    extern AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> bc_Xk_z_hi;

    extern AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> wallspeed_x_lo;
    extern AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> wallspeed_y_lo;
    extern AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> wallspeed_z_lo;
    extern AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> wallspeed_x_hi;
    extern AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> wallspeed_y_hi;
    extern AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> wallspeed_z_hi;

    extern AMREX_GPU_MANAGED amrex::Real bc_rhotot_x_lo;
    extern AMREX_GPU_MANAGED amrex::Real bc_rhotot_x_hi;
    extern AMREX_GPU_MANAGED amrex::Real bc_rhotot_y_lo;
    extern AMREX_GPU_MANAGED amrex::Real bc_rhotot_y_hi;
    extern AMREX_GPU_MANAGED amrex::Real bc_rhotot_z_lo;
    extern AMREX_GPU_MANAGED amrex::Real bc_rhotot_z_hi;

    // Each no-slip wall may be moving with a specified tangential
    // velocity along the tangential directions
    // In 2D:
    // wallspeed_lo/hi_x - yvel
    // wallspeed_lo/hi_y - xvel
    // In 3D:
    // wallspeed_lo/hi_x - yvel,zvel
    // wallspeed_lo/hi_y - xvel,zvel
    // wallspeed_lo/hi_z - xvel,yvel
    extern amrex::Vector<amrex::Real> wallspeed_lo;
    extern amrex::Vector<amrex::Real> wallspeed_hi;

    // specify the Dirichlet or Neumann BC value
    extern AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> potential_lo;
    extern AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> potential_hi;

    extern int                        n_burn;
    extern int                        dsmc_boundaries;
    extern amrex::Real                phonon_sound_speed;
    extern amrex::Real                tau_i;
    extern amrex::Real                tau_la;
    extern amrex::Real                tau_ta;
    extern int                        toggleTimeFrac;

    // structure factor and radial/cartesian pair correlation function analysis
    extern int                        struct_fact_int;
    extern int                        radialdist_int;
    extern int                        cartdist_int;
    extern int                        n_steps_skip;
    extern AMREX_GPU_MANAGED amrex::Real binSize;
    extern AMREX_GPU_MANAGED amrex::Real searchDist;

    // projection
    extern int                        project_dir;
    extern int                        slicepoint;
    extern amrex::Vector<int>         max_grid_projection;

    // These are mostly used for reaction-diffusion:
    extern int                        histogram_unit;
    extern amrex::Vector<amrex::Real> density_weights;
    extern amrex::Vector<int>         shift_cc_to_boundary;

    extern int                        particle_placement;
    extern int                        particle_input;
    extern amrex::Vector<int>         particle_count;
    extern amrex::Vector<int>         p_move_tog;
    extern amrex::Vector<int>         p_force_tog;
    extern amrex::Real                particle_neff;
    extern amrex::Vector<amrex::Real> particle_n0;
    extern amrex::Vector<amrex::Real> mass;
    extern amrex::Vector<amrex::Real> nfrac;

    extern AMREX_GPU_MANAGED amrex::Real permittivity;
    extern AMREX_GPU_MANAGED int      wall_mob;

    extern AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES*MAX_SPECIES> rmin;
    extern AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES*MAX_SPECIES> rmax;
    extern AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES*MAX_SPECIES> eepsilon;
    extern AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES*MAX_SPECIES> sigma;
    extern AMREX_GPU_MANAGED amrex::GpuArray<int, MAX_SPECIES*MAX_SPECIES>         p_int_tog;

    extern AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES*MAX_SPECIES> alpha_pp;
    extern AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES*MAX_SPECIES> alpha_pw;
    extern AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES*MAX_SPECIES> friction_pp;
    extern AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES*MAX_SPECIES> friction_pw;

    extern AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> phi_domain;

    extern AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> Yk0;

    extern AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> offset_wall;
    extern AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> rmin_wall;
    extern AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> rmax_wall;
    extern AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> eepsilon_wall;
    extern AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> sigma_wall;
    extern AMREX_GPU_MANAGED amrex::GpuArray<int, MAX_SPECIES> p_int_tog_wall;

    extern AMREX_GPU_MANAGED amrex::GpuArray<int, MAX_SPECIES> msd_int;
    extern AMREX_GPU_MANAGED amrex::GpuArray<int, MAX_SPECIES> msd_len;

    extern int                        poisson_verbose;
    extern int                        poisson_bottom_verbose;
    extern int                        poisson_max_iter;
    extern amrex::Real                poisson_rel_tol;

    extern amrex::Real                particle_grid_refine;
    extern amrex::Real                es_grid_refine;

    // total diffusion coefficient for immersed particle
    extern amrex::Vector<amrex::Real> diff;
    extern int                        all_dry;

    extern int                        fluid_tog;
    extern int                        es_tog;
    extern int                        drag_tog;
    extern int                        move_tog;
    extern int                        rfd_tog;
    extern AMREX_GPU_MANAGED int      dry_move_tog;
    extern AMREX_GPU_MANAGED int      sr_tog;
    extern int                        graphene_tog;
    extern int                        crange;
    extern int                        thermostat_tog;
    extern int                        zero_net_force;

    extern AMREX_GPU_MANAGED int      images;
    extern amrex::Vector<amrex::Real> eamp;
    extern amrex::Vector<amrex::Real> efreq;
    extern amrex::Vector<amrex::Real> ephase;

    extern amrex::Vector<amrex::Real> body_force_density;

    extern int                        plot_ascii;
    extern int                        plot_means;
    extern int                        plot_vars;
    extern int                        plot_mom3;
    extern int                        plot_covars;
    extern int                        plot_cross;
    extern int                        plot_deltaY_dir;
    extern int                        particle_motion;

    // parameters for turbulent forcing example
    extern AMREX_GPU_MANAGED amrex::Real turb_a;
    extern AMREX_GPU_MANAGED amrex::Real turb_b;
    extern AMREX_GPU_MANAGED amrex::Real turb_c;
    extern AMREX_GPU_MANAGED amrex::Real turb_d;
    extern AMREX_GPU_MANAGED amrex::Real turb_alpha;
    extern AMREX_GPU_MANAGED int         turbForcing;

}