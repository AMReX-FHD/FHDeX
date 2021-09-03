#include "common_functions.H"

#include <AMReX_ParmParse.H>

amrex::IntVect                common::nodal_flag;
amrex::Vector<amrex::IntVect> common::nodal_flag_dir;
amrex::Vector<amrex::IntVect> common::nodal_flag_edge;
amrex::IntVect                common::nodal_flag_x;
amrex::IntVect                common::nodal_flag_y;
amrex::IntVect                common::nodal_flag_z;
amrex::IntVect                common::nodal_flag_xy;
amrex::IntVect                common::nodal_flag_xz;
amrex::IntVect                common::nodal_flag_yz;

AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, 3> common::prob_lo;
AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, 3> common::prob_hi;
AMREX_GPU_MANAGED amrex::GpuArray<int, AMREX_SPACEDIM> common::n_cells;
amrex::Vector<int>         common::max_grid_size;
amrex::Vector<int>         common::max_particle_tile_size;
AMREX_GPU_MANAGED amrex::Real common::cell_depth;
amrex::IntVect             common::ngc;
AMREX_GPU_MANAGED int      common::nvars;
AMREX_GPU_MANAGED int      common::nprimvars;

AMREX_GPU_MANAGED int      common::membrane_cell;
AMREX_GPU_MANAGED int      common::cross_cell;
AMREX_GPU_MANAGED int      common::do_slab_sf;
AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> common::transmission;

AMREX_GPU_MANAGED amrex::GpuArray<int, MAX_SPECIES> common::pkernel_fluid; // GALEN - FLUID KERNEL
AMREX_GPU_MANAGED amrex::GpuArray<int, MAX_SPECIES> common::pkernel_es;
AMREX_GPU_MANAGED amrex::GpuArray<int, MAX_SPECIES> common::eskernel_fluid; // EXPONENTIAL OF SEMICIRCLE KERNEL; use this as w
AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> common::eskernel_beta; // EXPONENTIAL OF SEMICIRCLE KERNEL
amrex::Vector<amrex::Real> common::qval;

amrex::Real                common::fixed_dt;
amrex::Real                common::cfl;
amrex::Real                common::rfd_delta;
int                        common::max_step;
int                        common::plot_int;
int                        common::plot_stag;
std::string                common::plot_base_name;
int                        common::chk_int;
std::string                common::chk_base_name;
AMREX_GPU_MANAGED int      common::prob_type;
int                        common::restart;
int                        common::reset_stats;
int                        common::particle_restart;
int                        common::print_int;
int                        common::project_eos_int;
AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> common::grav;
AMREX_GPU_MANAGED int      common::nspecies;
AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> common::molmass;
AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> common::rhobar;
AMREX_GPU_MANAGED amrex::Real common::rho0;
AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> common::diameter;
AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> common::dof;
AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> common::e0;
AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> common::hcv;
AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> common::hcp;

AMREX_GPU_MANAGED amrex::Real common::variance_coef_mom;
AMREX_GPU_MANAGED amrex::Real common::variance_coef_mass;
AMREX_GPU_MANAGED amrex::Real common::k_B;
AMREX_GPU_MANAGED amrex::Real common::Runiv;
AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> common::T_init;
AMREX_GPU_MANAGED int      common::algorithm_type;
int                        common::barodiffusion_type;
int                        common::use_bl_rng;
int                        common::seed;
int                        common::seed_momentum;
int                        common::seed_diffusion;
int                        common::seed_reaction;
int                        common::seed_init_mass;
int                        common::seed_init_momentum;
AMREX_GPU_MANAGED amrex::Real common::visc_coef;
AMREX_GPU_MANAGED int      common::visc_type;
int                        common::advection_type;
int                        common::filtering_width;
int                        common::stoch_stress_form;
amrex::Vector<amrex::Real> common::u_init;
amrex::Real                common::perturb_width;
AMREX_GPU_MANAGED amrex::Real common::smoothing_width;
amrex::Real                common::initial_variance_mom;
amrex::Real                common::initial_variance_mass;
amrex::Real                common::domega;


AMREX_GPU_MANAGED amrex::GpuArray<int, AMREX_SPACEDIM>         common::bc_vel_lo;
AMREX_GPU_MANAGED amrex::GpuArray<int, AMREX_SPACEDIM>         common::bc_vel_hi;
AMREX_GPU_MANAGED amrex::GpuArray<int, AMREX_SPACEDIM>         common::bc_es_lo;
AMREX_GPU_MANAGED amrex::GpuArray<int, AMREX_SPACEDIM>         common::bc_es_hi;
AMREX_GPU_MANAGED amrex::GpuArray<int, AMREX_SPACEDIM>         common::bc_mass_lo;
AMREX_GPU_MANAGED amrex::GpuArray<int, AMREX_SPACEDIM>         common::bc_mass_hi;
AMREX_GPU_MANAGED amrex::GpuArray<int, AMREX_SPACEDIM>         common::bc_therm_lo;
AMREX_GPU_MANAGED amrex::GpuArray<int, AMREX_SPACEDIM>         common::bc_therm_hi;

AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> common::p_lo;
AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> common::p_hi;

AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> common::t_lo;
AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> common::t_hi;

AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> common::rho_lo;
AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> common::rho_hi;

AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> common::Yk0;

AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> common::bc_Yk_x_lo;
AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> common::bc_Yk_x_hi;
AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> common::bc_Yk_y_lo;
AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> common::bc_Yk_y_hi;
AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> common::bc_Yk_z_lo;
AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> common::bc_Yk_z_hi;

AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> common::n_lo;
AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> common::n_hi;

AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> common::bc_Xk_x_lo;
AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> common::bc_Xk_x_hi;
AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> common::bc_Xk_y_lo;
AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> common::bc_Xk_y_hi;
AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> common::bc_Xk_z_lo;
AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> common::bc_Xk_z_hi;

AMREX_GPU_MANAGED amrex::Real common::bc_rhotot_x_lo;
AMREX_GPU_MANAGED amrex::Real common::bc_rhotot_x_hi;
AMREX_GPU_MANAGED amrex::Real common::bc_rhotot_y_lo;
AMREX_GPU_MANAGED amrex::Real common::bc_rhotot_y_hi;
AMREX_GPU_MANAGED amrex::Real common::bc_rhotot_z_lo;
AMREX_GPU_MANAGED amrex::Real common::bc_rhotot_z_hi;

amrex::Vector<amrex::Real> common::wallspeed_lo;
amrex::Vector<amrex::Real> common::wallspeed_hi;

AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> common::potential_lo;
AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> common::potential_hi;

int                        common::dsmc_boundaries;

int                           common::struct_fact_int;
int                           common::radialdist_int;
int                           common::cartdist_int;
int                           common::n_steps_skip;
AMREX_GPU_MANAGED amrex::Real common::binSize;
AMREX_GPU_MANAGED amrex::Real common::searchDist;
int                           common::project_dir;
int                           common::slicepoint;
amrex::Vector<int>            common::max_grid_projection;
int                           common::histogram_unit;
amrex::Vector<amrex::Real>    common::density_weights;
amrex::Vector<int>            common::shift_cc_to_boundary;

int                           common::particle_placement;
int			      common::particle_input;
amrex::Vector<int>            common::particle_count;
amrex::Vector<int>            common::p_move_tog;
amrex::Vector<int>            common::p_force_tog;

amrex::Real                   common::particle_neff;
amrex::Vector<amrex::Real>    common::particle_n0;
amrex::Vector<amrex::Real>    common::mass;
amrex::Vector<amrex::Real>    common::nfrac;

AMREX_GPU_MANAGED amrex::GpuArray<int, MAX_SPECIES*MAX_SPECIES>         common::p_int_tog;
AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES*MAX_SPECIES> common::eepsilon;
AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES*MAX_SPECIES> common::sigma;
AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES*MAX_SPECIES> common::rmin;
AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES*MAX_SPECIES> common::rmax;

AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES*MAX_SPECIES> common::alpha_pp;
AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES*MAX_SPECIES> common::alpha_pw;
AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES*MAX_SPECIES> common::friction_pp;
AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES*MAX_SPECIES> common::friction_pw;

AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> common::phi_domain;

AMREX_GPU_MANAGED amrex::GpuArray<int, MAX_SPECIES>         common::p_int_tog_wall;
AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> common::eepsilon_wall;
AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> common::sigma_wall;
AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> common::rmin_wall;
AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> common::rmax_wall;

int                        common::poisson_verbose;
int                        common::poisson_bottom_verbose;
int                        common::poisson_max_iter;

amrex::Real                common::poisson_rel_tol;
AMREX_GPU_MANAGED amrex::Real common::permittivity;
AMREX_GPU_MANAGED int      common::wall_mob;


amrex::Real                common::particle_grid_refine;
amrex::Real                common::es_grid_refine;
amrex::Vector<amrex::Real> common::diff;
int                        common::all_dry;

int                        common::fluid_tog;
int                        common::es_tog;
int                        common::drag_tog;
int                        common::move_tog;
int                        common::rfd_tog;
AMREX_GPU_MANAGED int      common::dry_move_tog;
AMREX_GPU_MANAGED int      common::sr_tog;
int                        common::graphene_tog;
int	                   common::thermostat_tog;
int	                   common::zero_net_force;

int                        common::crange;

AMREX_GPU_MANAGED int      common::images;
amrex::Vector<amrex::Real> common::eamp;
amrex::Vector<amrex::Real> common::efreq;
amrex::Vector<amrex::Real> common::ephase;

int                        common::plot_ascii;
int                        common::plot_means;
int                        common::plot_vars;
int                        common::plot_covars;
int                        common::plot_cross;
int                        common::particle_motion;
int                        common::solve_chem;
amrex::Real                common::diffcoeff;
amrex::Real                common::scaling_factor;
amrex::Real                common::source_strength;
int                        common::regrid_int;
int                        common::do_reflux;

amrex::Real                common::turb_a;
amrex::Real                common::turb_b;
int                        common::turbForcing;

AMREX_GPU_MANAGED int      common::do_1D;



void InitializeCommonNamespace() {

    BL_PROFILE_VAR("InitializeCommonNamespace()",InitializeCommonNameSpace);
    
    nodal_flag_dir.resize(AMREX_SPACEDIM);
    nodal_flag_edge.resize(AMREX_SPACEDIM);

    for (int i=0; i<AMREX_SPACEDIM; ++i) {

        //_______________________________________________________________________
        // Designates data on nodes
        nodal_flag[i] = 1;

        //_______________________________________________________________________
        // Designates data on faces
        nodal_flag_x[i] = int(i==0);
        nodal_flag_y[i] = int(i==1);
        nodal_flag_z[i] = int(i==2);

        // Enable indexing flags above in loops
        AMREX_D_TERM(nodal_flag_dir[0][i] = nodal_flag_x[i];,
                     nodal_flag_dir[1][i] = nodal_flag_y[i];,
                     nodal_flag_dir[2][i] = nodal_flag_z[i];);

        //_______________________________________________________________________
        // Designates data on edges
        nodal_flag_xy[i] = int(i==0 || i==1);
        nodal_flag_xz[i] = int(i==0 || i==2);
        nodal_flag_yz[i] = int(i==1 || i==2);

        // Enable indexing flags above in loops
        AMREX_D_TERM(nodal_flag_edge[0][i] = nodal_flag_xy[i];,
                     nodal_flag_edge[1][i] = nodal_flag_xz[i];,
                     nodal_flag_edge[2][i] = nodal_flag_yz[i];);
    }

    max_grid_size.resize(AMREX_SPACEDIM);
    max_particle_tile_size.resize(AMREX_SPACEDIM);
    u_init.resize(2);

    wallspeed_lo.resize((AMREX_SPACEDIM-1)*AMREX_SPACEDIM);
    wallspeed_hi.resize((AMREX_SPACEDIM-1)*AMREX_SPACEDIM);

    max_grid_projection.resize(AMREX_SPACEDIM-1);

    density_weights.resize(MAX_SPECIES);
    shift_cc_to_boundary.resize(AMREX_SPACEDIM*2);

    mass.resize(MAX_SPECIES);
    nfrac.resize(MAX_SPECIES);
    particle_count.resize(MAX_SPECIES);
    p_move_tog.resize(MAX_SPECIES);
    p_force_tog.resize(MAX_SPECIES);
    particle_n0.resize(MAX_SPECIES);

    qval.resize(MAX_SPECIES);

    diff.resize(MAX_SPECIES);

    eamp.resize(3);    
    efreq.resize(3);
    ephase.resize(3);

    // specify default values first, then read in values from inputs file

    // Problem specification
    for (int i=0; i<3; ++i) {
        prob_lo[i] = 0.; // physical lo coordinate
        prob_hi[i] = 0.; // physical hi coordinate
    }

    for (int i=0; i<AMREX_SPACEDIM; ++i) {
        n_cells[1] = 1;                // number of cells in domain
        max_grid_size[i] = 1;          // max number of cells in a box
        max_particle_tile_size[i] = 0;
    }
    
    cell_depth = 1.;

    for (int i=0; i<AMREX_SPACEDIM; ++i) {
        ngc[i] = 1;           // number of ghost cells
    }
    
    // nvars - number of conserved variables (no default)
    // primvars - number of primative variables (no default)

    membrane_cell = -1; // location of membrane
    cross_cell = 0;     // cell to compute spatial correlation
    do_slab_sf = 0;     // whether to compute SF in two slabs separated by cross_cell
    // transmission - no default

    for (int i=0; i<MAX_SPECIES; ++i) {
        qval[i] = 0.;                // charge on an ion
        pkernel_fluid[i] = 4;        // peskin kernel for fluid
        pkernel_es[i] = 4;           // peskin kernel for es
        eskernel_fluid[i] = -1;      // ES kernel for fluid
        eskernel_beta[i] = -1;       // ES kernel for fluid: beta
    }

    // mass (no default)
    // nfrac (no default)

    // particle_placement (no default)
    particle_input = -1;
    for (int i=0; i<MAX_SPECIES; ++i) {
        particle_count[i] = -1.;
        p_move_tog[i] = 1.;
        p_force_tog[i] = 1.;
        p_int_tog[i] = 1.;
        particle_n0[i] = -1.;
    }
    
    // p_int_tog_wall (no default)
    particle_neff = 1;

    // Time-step control
    fixed_dt = 1.;
    cfl = 0.5;

    //random finite difference size, fraction of cell size
    rfd_delta = 1.e-5;

    // Controls for number of steps between actions
    max_step = 1;
    plot_int = 0;
    plot_stag = 0;
    plot_base_name = "plt";
    chk_int = 0;
    chk_base_name = "chk";
    prob_type = 1;
    restart = -1;
    reset_stats = 0;
    particle_restart = -1;
    print_int = 0;
    project_eos_int = -1;

    // Physical parameters
    for (int i=0; i<AMREX_SPACEDIM; ++i) {
        grav[i] = 0.;
    }
    for (int i=0; i<MAX_SPECIES; ++i) {
        rhobar[i] = 1.;
    }
    rho0 = 1.;

    // Kinetic parameters
    nspecies = 2;
    for (int i=0; i<MAX_SPECIES; ++i) {
        molmass[i] = 1.;
        diameter[i] = 1.;
    }

    // dof (no default)
    for (int i=0; i<MAX_SPECIES; ++i) {
        e0[i] = 0.;
    }
    // hcv (no default)
    // hcp (no default)

    // stochastic forcing amplitudes (1 for physical values, 0 to run them off)
    variance_coef_mom = 1.;
    variance_coef_mass = 1.;
    k_B = 1.38064852e-16;
    Runiv = 8.314462175e7;
    for (int i=0; i<MAX_SPECIES; ++i) {
        T_init[i] = 1.;
    }

    // Algorithm control / selection
    algorithm_type = 0;
    advection_type = 0;
    barodiffusion_type = 0;
    use_bl_rng = 0;

    // random number seed
    // 0        = unpredictable seed based on clock
    // positive = fixed seed
    seed = 0;

    // as assortment of other seeds in case one needs different engines
    // implementation is problem-dependent
    // 0        = unpredictable seed based on clock
    // positive = fixed seed
    seed_momentum = 1;
    seed_diffusion = 1;
    seed_reaction = 1;
    seed_init_mass = 1;
    seed_init_momentum = 1;

    // Viscous friction L phi operator
    // if abs(visc_type) = 1, L = div beta grad
    // if abs(visc_type) = 2, L = div [ beta (grad + grad^T) ]
    // if abs(visc_type) = 3, L = div [ beta (grad + grad^T) + I (gamma - (2/3)*beta) div ]
    // positive = assume constant coefficients
    // negative = assume spatially-varying coefficients
    visc_coef = 1.;
    visc_type = 1;

    // Stochastic momentum flux controls:
    filtering_width = 0;
    stoch_stress_form = 1;

    // Initial conditions
    u_init[0] = 0.;
    u_init[1] = 0.;
    perturb_width = 0.;
    smoothing_width = 1.;
    initial_variance_mom = 0.;
    initial_variance_mass = 0.;
    domega = 0.;

    // Boundary conditions
    for (int i=0; i<AMREX_SPACEDIM; ++i) {
        bc_vel_lo[i] = 0;
        bc_vel_hi[i] = 0;
        bc_es_lo[i] = 0;
        bc_es_hi[i] = 0;
        bc_mass_lo[i] = 0;
        bc_mass_hi[i] = 0;
        bc_therm_lo[i] = 0;
        bc_therm_hi[i] = 0;

        // Pressure drop are periodic inflow/outflow walls (bc_[hi,lo]=-2).
        p_lo[i] = 0.;
        p_hi[i] = 0.;

        t_lo[i] = 0.;
        t_hi[i] = 0.;
  
        rho_lo[i] = 0.;
        rho_hi[i] = 0.;
    } 

    // c_i boundary conditions
    for (int i=0; i<MAX_SPECIES; ++i) {
        bc_Yk_x_lo[i] = -1.;
        bc_Yk_x_hi[i] = -1.;
        bc_Yk_y_lo[i] = -1.;
        bc_Yk_y_hi[i] = -1.;
        bc_Yk_z_lo[i] = -1.;
        bc_Yk_z_hi[i] = -1.;

        bc_Xk_x_lo[i] = -1.;
        bc_Xk_x_hi[i] = -1.;
        bc_Xk_y_lo[i] = -1.;
        bc_Xk_y_hi[i] = -1.;
        bc_Xk_z_lo[i] = -1.;
        bc_Xk_z_hi[i] = -1.;
    }

    for (int i=0; i<AMREX_SPACEDIM; ++i) {
        n_lo[i] = -1.;
        n_hi[i] = -1.;
    }

    // Each no-slip wall may be moving with a specified tangential
    for (int i=0; i<(AMREX_SPACEDIM-1)*AMREX_SPACEDIM; ++i) {
        wallspeed_lo[i] = 0.;
        wallspeed_hi[i] = 0.;
    }

    for (int i=0; i<AMREX_SPACEDIM; ++i) {
        potential_lo[i] = 0.;
        potential_hi[i] = 0.;
    }

    dsmc_boundaries = 0;

    // structure factor and radial/cartesian pair correlation function analysis
    struct_fact_int = 0;
    radialdist_int = 0;
    cartdist_int = 0;
    n_steps_skip = 0;
    binSize = 0.;
    searchDist = 0.;

    // projection
    project_dir = -1;
    slicepoint = -1;
    for (int i=0; i<AMREX_SPACEDIM-1; ++i) {
        max_grid_projection[i] = 4096;
    }

    // These are mostly used for reaction-diffusion:
    histogram_unit = 0;
    for (int i=0; i<MAX_SPECIES; ++i) {
        density_weights[i] = 0.;
    }
    for (int i=0; i<2*AMREX_SPACEDIM; ++i) {
        shift_cc_to_boundary[i] = 0;
    }

    // permittivity (no default)
    wall_mob = 1;
    // rmin (no default)
    // rmax (no default)
    // eepsilon (no default)
    // sigma (no default)

    // rmin_wall (no default)
    // rmax_wall (no default)
    // eepsilon_wall (no default)
    // sigma_wall (no default)
  
    poisson_verbose = 1;
    poisson_bottom_verbose = 0;
    poisson_max_iter = 100;
    poisson_rel_tol = 1.e-10;

    particle_grid_refine = 1;
    es_grid_refine = 1;
    // diff (no default)
    all_dry = 0;

    // fluid_tog (no default)
    // es_tog (no default)
    drag_tog = 0;
    // move_tog (no default)
    // rfd_tog (no default)
    // dry_move_tog (no default)
    // sr_tog (no default)
    graphene_tog = 0;
    crange = 5;
    thermostat_tog = 0;
    zero_net_force = 0;

    // images (no default)
    for (int i=0; i<3; ++i) {
        eamp[i] = 0.;
        efreq[i] = 0.;
        ephase[i] = 0.;
    }

    // plot_ascii (no default)
    plot_means = 0;
    plot_vars = 0;
    plot_covars = 0;
    plot_cross = 0;
    particle_motion = 0;

    // chemistry
    solve_chem = 0.;
    diffcoeff = 1.e-3;
    source_strength = 0.1;
    scaling_factor = 0.1;
    regrid_int = 25;
    do_reflux = 0;

    // turblent forcing parameters
    turb_a = 1.;
    turb_b = 1.;
    turbForcing = 0;

  // DSMC Granular
    for (int i=0; i<MAX_SPECIES*MAX_SPECIES; ++i) {
        alpha_pp[i] = 1.;
        alpha_pw[i] = 1.;
        friction_pp[i] = 0.;
        friction_pw[i] = 0.;
    }
    for (int i=0; i<MAX_SPECIES; ++i) {
        phi_domain[i] = -1.;
        Yk0[i] = 0.;
    }
  
    do_1D = 0;
    
    char temp_plot_base_name[128];
    char temp_chk_base_name[128];

    initialize_common_namespace(prob_lo.begin(), prob_hi.begin(), n_cells.data(),
                                max_grid_size.dataPtr(),
                                max_particle_tile_size.dataPtr(), &cell_depth, ngc.getVect(),
                                &nvars, &nprimvars,
                                &membrane_cell, &cross_cell, &do_slab_sf, transmission.data(),
                                qval.dataPtr(), pkernel_fluid.begin(), pkernel_es.begin(),
                                eskernel_fluid.begin(), eskernel_beta.begin(),
                                &fixed_dt, &cfl, &rfd_delta, &max_step,
                                &plot_int, &plot_stag, temp_plot_base_name, 128,
                                &chk_int, temp_chk_base_name, 128,
                                &prob_type, &restart, &reset_stats, &particle_restart, &print_int, &project_eos_int,
                                grav.data(), &nspecies, molmass.data(), diameter.data(),
                                dof.data(), e0.data(), hcv.data(), hcp.data(),
                                rhobar.data(),
                                &rho0, &variance_coef_mom, &variance_coef_mass, &k_B, &Runiv,
                                T_init.begin(),
                                &algorithm_type,  &advection_type,
                                &barodiffusion_type, &use_bl_rng, &seed,
                                &seed_momentum, &seed_diffusion, &seed_reaction,
                                &seed_init_mass,
                                &seed_init_momentum, &visc_coef, &visc_type,
                                &filtering_width, &stoch_stress_form, u_init.dataPtr(),
                                &perturb_width, &smoothing_width, &initial_variance_mom,
                                &initial_variance_mass, &domega,
                                bc_vel_lo.data(), bc_vel_hi.data(),
                                bc_es_lo.data(), bc_es_hi.data(),
                                bc_mass_lo.data(), bc_mass_hi.data(),
                                bc_therm_lo.data(), bc_therm_hi.data(),
                                p_lo.data(), p_hi.data(),
                                t_lo.data(), t_hi.data(),
                                rho_lo.data(), rho_hi.data(),
                                bc_Yk_x_lo.data(), bc_Yk_x_hi.data(),
                                bc_Yk_y_lo.data(), bc_Yk_y_hi.data(),
                                bc_Yk_z_lo.data(), bc_Yk_z_hi.data(),
                                n_lo.data(), n_hi.data(),
                                bc_Xk_x_lo.data(), bc_Xk_x_hi.data(),
                                bc_Xk_y_lo.data(), bc_Xk_y_hi.data(),
                                bc_Xk_z_lo.data(), bc_Xk_z_hi.data(),
                                wallspeed_lo.dataPtr(), wallspeed_hi.dataPtr(),
                                potential_lo.data(), potential_hi.data(),
                                &struct_fact_int, &dsmc_boundaries, &radialdist_int, &cartdist_int,
                                &n_steps_skip, &binSize, &searchDist,
                                &project_dir, &slicepoint, max_grid_projection.dataPtr(),
                                &histogram_unit,
                                density_weights.dataPtr(), shift_cc_to_boundary.dataPtr(),
                                &particle_placement, &particle_input, particle_count.dataPtr(),
                                p_move_tog.dataPtr(), p_force_tog.dataPtr(),
                                p_int_tog.begin(),p_int_tog_wall.begin(), &particle_neff,
                                particle_n0.dataPtr(), mass.dataPtr(), nfrac.dataPtr(),
                                &permittivity,
                                &wall_mob,rmin.begin(),rmax.begin(), eepsilon.begin(), 
                                sigma.begin(),rmin_wall.begin(),rmax_wall.begin(),
                                eepsilon_wall.begin(), sigma_wall.begin(),
                                &poisson_verbose, &poisson_bottom_verbose, &poisson_max_iter,
                                &poisson_rel_tol, &particle_grid_refine, &es_grid_refine,
                                diff.dataPtr(), &all_dry, &fluid_tog, &es_tog, &drag_tog, &move_tog, &rfd_tog,
                                &dry_move_tog, &sr_tog, &graphene_tog, &crange, &thermostat_tog, &zero_net_force,
                                &images, eamp.dataPtr(), efreq.dataPtr(), ephase.dataPtr(),
                                &plot_ascii, &plot_means, &plot_vars, &plot_covars, &plot_cross,
                                &solve_chem, &diffcoeff, &scaling_factor,
                                &source_strength, &regrid_int, &do_reflux, &particle_motion,
                                &turb_a, &turb_b, &turbForcing,
                                alpha_pp.begin(), alpha_pw.begin(),
                                friction_pp.begin(), friction_pw.begin(),
                                phi_domain.begin(), Yk0.begin(),
                                &do_1D);

    plot_base_name = temp_plot_base_name;
    chk_base_name = temp_chk_base_name;

    if (nspecies > MAX_SPECIES) {
        Abort("InitializeCommonNamespace: nspecies > MAX_SPECIES");
    }
    
    ParmParse pp;

    // read in from command line
    pp.query("max_step",max_step);

    // copy value into fortran namelist
    set_max_step(&max_step);

    // read in from command line
    pp.query("domega",domega);

    // copy value into fortran namelist
    set_domega(&domega);
    
}
