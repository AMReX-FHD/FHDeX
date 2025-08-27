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

AMREX_GPU_MANAGED int      common::cross_cell;
AMREX_GPU_MANAGED int      common::do_slab_sf;

AMREX_GPU_MANAGED amrex::GpuArray<int, MAX_SPECIES> common::pkernel_fluid; // GALEN - FLUID KERNEL
AMREX_GPU_MANAGED amrex::GpuArray<int, MAX_SPECIES> common::pkernel_es;
AMREX_GPU_MANAGED amrex::GpuArray<int, MAX_SPECIES> common::eskernel_fluid; // EXPONENTIAL OF SEMICIRCLE KERNEL; use this as w
AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> common::eskernel_beta; // EXPONENTIAL OF SEMICIRCLE KERNEL
amrex::Vector<amrex::Real> common::qval;

amrex::Real                common::fixed_dt;
amrex::Real                common::cfl;
amrex::Real                common::rfd_delta;
int                        common::max_step;
int                        common::ramp_step;
int                        common::plot_int;
int                        common::plot_stag;
std::string                common::plot_base_name;
int                        common::chk_int;
std::string                common::chk_base_name;
std::string                common::plot_init_file;
AMREX_GPU_MANAGED int      common::prob_type;
int                        common::restart;
int                        common::stats_int;
int                        common::reset_stats;
int                        common::particle_restart;
int                        common::print_int;
int                        common::project_eos_int;
AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> common::grav;
AMREX_GPU_MANAGED int      common::nspecies;
AMREX_GPU_MANAGED int      common::nbonds;
AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> common::molmass;
AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> common::rhobar;
AMREX_GPU_MANAGED amrex::Real common::rho0;
AMREX_GPU_MANAGED amrex::Real common::mach0;
AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> common::diameter;
AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> common::dof;
AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> common::e0;
AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> common::hcv;
AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> common::hcp;

AMREX_GPU_MANAGED amrex::Real common::variance_coef_mom;
AMREX_GPU_MANAGED amrex::Real common::variance_coef_mass;
AMREX_GPU_MANAGED amrex::Real common::variance_coef_ener;
AMREX_GPU_MANAGED amrex::Real common::k_B;
AMREX_GPU_MANAGED amrex::Real common::h_bar;
AMREX_GPU_MANAGED amrex::Real common::Runiv;
AMREX_GPU_MANAGED amrex::Real common::avogadro;
AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> common::T_init;
AMREX_GPU_MANAGED int      common::algorithm_type;
int                        common::barodiffusion_type;
int                        common::seed;
AMREX_GPU_MANAGED amrex::Real common::visc_coef;
AMREX_GPU_MANAGED int      common::visc_type;
AMREX_GPU_MANAGED int      common::advection_type;
int                        common::filtering_width;
int                        common::stoch_stress_form;
amrex::Vector<amrex::Real> common::u_init;
amrex::Real                common::perturb_width;
AMREX_GPU_MANAGED amrex::Real common::smoothing_width;
AMREX_GPU_MANAGED amrex::Real common::radius_cyl;
AMREX_GPU_MANAGED amrex::Real common::radius_outer;
AMREX_GPU_MANAGED amrex::Real common::film_thickness;
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

AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> common::wallspeed_x_lo;
AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> common::wallspeed_y_lo;
AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> common::wallspeed_z_lo;
AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> common::wallspeed_x_hi;
AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> common::wallspeed_y_hi;
AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> common::wallspeed_z_hi;

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

int                           common::n_burn;

int                           common::dsmc_boundaries;
amrex::Real                   common::phonon_sound_speed;
amrex::Real                   common::tau_ta;
amrex::Real                   common::tau_la;
amrex::Real                   common::tau_i;
int                           common::toggleTimeFrac;

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
int                           common::particle_input;
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
AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> common::offset_wall;
AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> common::rmax_wall;

AMREX_GPU_MANAGED amrex::GpuArray<int, MAX_SPECIES>         common::msd_int;
AMREX_GPU_MANAGED amrex::GpuArray<int, MAX_SPECIES>         common::msd_len;

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
int                        common::thermostat_tog;
int                        common::zero_net_force;

int                        common::crange;

AMREX_GPU_MANAGED int      common::images;
amrex::Vector<amrex::Real> common::eamp;
amrex::Vector<amrex::Real> common::efreq;
amrex::Vector<amrex::Real> common::ephase;
amrex::Vector<amrex::Real> common::body_force_density;

int                        common::plot_ascii;
int                        common::plot_means;
int                        common::plot_vars;
int                        common::plot_mom3;
int                        common::plot_covars;
int                        common::plot_cross;
int                        common::plot_deltaY_dir;
int                        common::particle_motion;

AMREX_GPU_MANAGED amrex::Real common::turb_a;
AMREX_GPU_MANAGED amrex::Real common::turb_b;
AMREX_GPU_MANAGED amrex::Real common::turb_c;
AMREX_GPU_MANAGED amrex::Real common::turb_d;
AMREX_GPU_MANAGED amrex::Real common::turb_alpha;
AMREX_GPU_MANAGED int         common::turbForcing;


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
    body_force_density.resize(3);


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

    cross_cell = 0;     // cell to compute spatial correlation
    do_slab_sf = 0;     // whether to compute SF in two slabs separated by membrane_cell

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


    for (int i=0; i<MAX_SPECIES; ++i) {
        msd_int[i] = 0;
        msd_len[i] = 0;
    }

    // Time-step control
    fixed_dt = 1.;
    cfl = 0.5;

    //random finite difference size, fraction of cell size
    rfd_delta = 1.e-5;

    // Controls for number of steps between actions
    max_step = 1;
    ramp_step = 0;
    plot_int = 0;
    plot_stag = 0;
    plot_base_name = "plt";
    chk_int = 0;
    chk_base_name = "chk";
    plot_init_file = "";
    prob_type = 1;
    restart = -1;
    stats_int = 1;
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
    mach0 = -1.0;

    // Kinetic parameters
    nspecies = 2;
    for (int i=0; i<MAX_SPECIES; ++i) {
        molmass[i] = 1.;
        diameter[i] = 1.;
    }
    nbonds = 0;

    // dof (no default)
    for (int i=0; i<MAX_SPECIES; ++i) {
        e0[i] = 0.;
    }
    // hcv (no default)
    // hcp (no default)

    // stochastic forcing amplitudes (1 for physical values, 0 to run them off)
    variance_coef_mom = 1.;
    variance_coef_mass = 1.;
    variance_coef_ener = 1.;
    k_B = 1.38064852e-16;
    h_bar = 1.0546e-27;
    Runiv = 8.314462175e7;
    for (int i=0; i<MAX_SPECIES; ++i) {
        T_init[i] = 1.;
    }

    // Algorithm control / selection
    algorithm_type = 0;
    // 0 = centered
    // 1 = unlimited bilinear bds
    // 2 = limited bilinear bds
    advection_type = 0;
    barodiffusion_type = 0;

    // random number seed
    // 0        = unpredictable seed based on clock
    // positive = fixed seed
    seed = 0;

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
    radius_cyl = 0.;
    radius_outer = 0.;
    film_thickness = 0.;
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

        rho_lo[i] = -1.;
        rho_hi[i] = -1.;
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
        wallspeed_x_lo[i] = 0.;
        wallspeed_y_lo[i] = 0.;
        wallspeed_z_lo[i] = 0.;
        wallspeed_x_hi[i] = 0.;
        wallspeed_y_hi[i] = 0.;
        wallspeed_z_hi[i] = 0.;
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
    n_burn = 1000;
    phonon_sound_speed = 600000.0;
    tau_i = 2.95e45;
    tau_ta = 9.3e13;
    tau_la = 2.0e24;
    toggleTimeFrac = 1;

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
        body_force_density[i] = 0.;

    }

    // plot_ascii (no default)
    plot_means = 0;
    plot_vars = 0;
    plot_mom3 = 0;
    plot_covars = 0;
    plot_cross = 0;
    plot_deltaY_dir = -1;
    particle_motion = 0;

    // turblent forcing parameters
    turb_a = 1.;
    turb_b = 1.;
    turb_c = 1.;
    turb_d = 1.;
    turb_alpha = 1.;
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

    ParmParse pp;

    int temp_max = std::max(3,MAX_SPECIES*MAX_SPECIES);

    amrex::Vector<amrex::Real> temp    (temp_max,0.);
    amrex::Vector<int>         temp_int(temp_max,0 );

    // pp.query searches for optional parameters
    // pp.get aborts if the parameter is not found
    // pp.getarr and queryarr("string",inputs,start_indx,count); can be used for arrays

    pp.query("nspecies",nspecies);
    if (nspecies > MAX_SPECIES) {
        Abort("nspecies > MAX_SPECIES; recompile with a new MAX_SPEC in the GNUmakefile");
    }
    pp.query("nbonds",nbonds);

    if (pp.queryarr("prob_lo",temp)) {
        for (int i=0; i<3; ++i) {
            prob_lo[i] = temp[i];
        }
    }
    if (pp.queryarr("prob_hi",temp)) {
        for (int i=0; i<3; ++i) {
            prob_hi[i] = temp[i];
        }
    }
    if (pp.queryarr("n_cells",temp_int,0,AMREX_SPACEDIM)) {
        for (int i=0; i<AMREX_SPACEDIM; ++i) {
            n_cells[i] = temp_int[i];
        }
    }
    pp.queryarr("max_grid_size",max_grid_size,0,AMREX_SPACEDIM);
    pp.queryarr("max_particle_tile_size",max_particle_tile_size,0,AMREX_SPACEDIM);
    pp.query("cell_depth",cell_depth);
    if (pp.queryarr("ngc",temp_int,0,AMREX_SPACEDIM)) {
        for (int i=0; i<AMREX_SPACEDIM; ++i) {
            ngc[i] = temp_int[i];
        }
    }
    pp.query("nvars",nvars);
    pp.query("nprimvars",nprimvars);
    pp.query("cross_cell",cross_cell);
    pp.query("do_slab_sf",do_slab_sf);

    if (pp.queryarr("pkernel_fluid",temp_int,0,nspecies)) {
        for (int i=0; i<nspecies; ++i) {
            pkernel_fluid[i] = temp_int[i];
        }
    }
    if (pp.queryarr("pkernel_es",temp_int,0,nspecies)) {
        for (int i=0; i<nspecies; ++i) {
            pkernel_es[i] = temp_int[i];
        }
    }
    if (pp.queryarr("eskernel_fluid",temp_int,0,nspecies)) {
        for (int i=0; i<nspecies; ++i) {
            eskernel_fluid[i] = temp_int[i];
        }
    }
    if (pp.queryarr("eskernel_beta",temp,0,nspecies)) {
        for (int i=0; i<nspecies; ++i) {
            eskernel_beta[i] = temp[i];
        }
    }
    pp.queryarr("qval",qval,0,nspecies);
    pp.query("fixed_dt",fixed_dt);
    pp.query("cfl",cfl);
    pp.query("rfd_delta",rfd_delta);
    pp.query("max_step",max_step);
    pp.query("ramp_step",ramp_step);
    pp.query("plot_int",plot_int);
    pp.query("plot_stag",plot_stag);
    pp.query("plot_base_name",plot_base_name);
    pp.query("chk_int",chk_int);
    pp.query("chk_base_name",chk_base_name);
    pp.query("plot_init_file",plot_init_file);
    pp.query("prob_type",prob_type);
    pp.query("restart",restart);
    pp.query("stats_int",stats_int);
    pp.query("reset_stats",reset_stats);
    pp.query("particle_restart",particle_restart);
    pp.query("print_int",print_int);
    pp.query("project_eos_int",project_eos_int);
    if (pp.queryarr("grav",temp,0,AMREX_SPACEDIM)) {
        for (int i=0; i<AMREX_SPACEDIM; ++i) {
            grav[i] = temp[i];
        }
    }
    if (pp.queryarr("molmass",temp,0,nspecies)) {
        for (int i=0; i<nspecies; ++i) {
            molmass[i] = temp[i];
        }
    }
    if (pp.queryarr("rhobar",temp,0,nspecies)) {
        for (int i=0; i<nspecies; ++i) {
            rhobar[i] = temp[i];
        }
    }
    pp.query("rho0",rho0);
    pp.query("mach0",mach0);
    if (pp.queryarr("diameter",temp,0,nspecies)) {
        for (int i=0; i<nspecies; ++i) {
            diameter[i] = temp[i];
        }
    }
    if (pp.queryarr("dof",temp,0,nspecies)) {
        for (int i=0; i<nspecies; ++i) {
            dof[i] = temp[i];
        }
    }
    if (pp.queryarr("e0",temp,0,nspecies)) {
        for (int i=0; i<nspecies; ++i) {
            e0[i] = temp[i];
        }
    }
    if (pp.queryarr("hcv",temp,0,nspecies)) {
        for (int i=0; i<nspecies; ++i) {
            hcv[i] = temp[i];
        }
    }
    if (pp.queryarr("hcp",temp,0,nspecies)) {
        for (int i=0; i<nspecies; ++i) {
            hcp[i] = temp[i];
        }
    }
    pp.query("variance_coef_mom",variance_coef_mom);
    pp.query("variance_coef_mass",variance_coef_mass);
    pp.query("variance_coef_ener",variance_coef_ener);
    pp.query("k_B",k_B);
    pp.query("h_bar",h_bar);
    pp.query("Runiv",Runiv);
    avogadro = Runiv / k_B;
    if (pp.query("avogadro",avogadro) ) {
        Runiv = k_B * avogadro;
    }
    if (pp.query("Runiv",Runiv) && pp.query("avogadro",avogadro)) {
        Abort("Cannot specify both Runiv and avogadro");
    }
    if (pp.queryarr("T_init",temp)) {
        for (int i=0; i<nspecies; ++i) {
            T_init[i] = temp[i];
        }
    }
    pp.query("algorithm_type",algorithm_type);
    pp.query("barodiffusion_type",barodiffusion_type);
    pp.query("seed",seed);
    pp.query("visc_coef",visc_coef);
    pp.query("visc_type",visc_type);
    pp.query("advection_type",advection_type);
    pp.query("filtering_width",filtering_width);
    pp.query("stoch_stress_form",stoch_stress_form);
    pp.queryarr("u_init",u_init,0,2);
    pp.query("perturb_width",perturb_width);
    pp.query("smoothing_width",smoothing_width);
    pp.query("radius_cyl",radius_cyl);
    pp.query("radius_outer",radius_outer);
    pp.query("film_thickness",film_thickness);
    pp.query("initial_variance_mom",initial_variance_mom);
    pp.query("initial_variance_mass",initial_variance_mass);
    pp.query("domega",domega);
    if (pp.queryarr("bc_vel_lo",temp_int,0,AMREX_SPACEDIM)) {
        for (int i=0; i<AMREX_SPACEDIM; ++i) {
            bc_vel_lo[i] = temp_int[i];
        }
    }
    if (pp.queryarr("bc_vel_hi",temp_int,0,AMREX_SPACEDIM)) {
        for (int i=0; i<AMREX_SPACEDIM; ++i) {
            bc_vel_hi[i] = temp_int[i];
        }
    }
    if (pp.queryarr("bc_es_lo",temp_int,0,AMREX_SPACEDIM)) {
        for (int i=0; i<AMREX_SPACEDIM; ++i) {
            bc_es_lo[i] = temp_int[i];
        }
    }
    if (pp.queryarr("bc_es_hi",temp_int,0,AMREX_SPACEDIM)) {
        for (int i=0; i<AMREX_SPACEDIM; ++i) {
            bc_es_hi[i] = temp_int[i];
        }
    }
    if (pp.queryarr("bc_mass_lo",temp_int,0,AMREX_SPACEDIM)) {
        for (int i=0; i<AMREX_SPACEDIM; ++i) {
            bc_mass_lo[i] = temp_int[i];
        }
    }
    if (pp.queryarr("bc_mass_hi",temp_int,0,AMREX_SPACEDIM)) {
        for (int i=0; i<AMREX_SPACEDIM; ++i) {
            bc_mass_hi[i] = temp_int[i];
        }
    }
    if (pp.queryarr("bc_therm_lo",temp_int,0,AMREX_SPACEDIM)) {
        for (int i=0; i<AMREX_SPACEDIM; ++i) {
            bc_therm_lo[i] = temp_int[i];
        }
    }
    if (pp.queryarr("bc_therm_hi",temp_int,0,AMREX_SPACEDIM)) {
        for (int i=0; i<AMREX_SPACEDIM; ++i) {
            bc_therm_hi[i] = temp_int[i];
        }
    }
    if (pp.queryarr("p_lo",temp,0,AMREX_SPACEDIM)) {
        for (int i=0; i<AMREX_SPACEDIM; ++i) {
            p_lo[i] = temp[i];
        }
    }
    if (pp.queryarr("p_hi",temp,0,AMREX_SPACEDIM)) {
        for (int i=0; i<AMREX_SPACEDIM; ++i) {
            p_hi[i] = temp[i];
        }
    }
    if (pp.queryarr("t_lo",temp,0,AMREX_SPACEDIM)) {
        for (int i=0; i<AMREX_SPACEDIM; ++i) {
            t_lo[i] = temp[i];
        }
    }
    if (pp.queryarr("t_hi",temp,0,AMREX_SPACEDIM)) {
        for (int i=0; i<AMREX_SPACEDIM; ++i) {
            t_hi[i] = temp[i];
        }
    }
    if (pp.queryarr("rho_lo",temp,0,AMREX_SPACEDIM)) {
        for (int i=0; i<AMREX_SPACEDIM; ++i) {
            rho_lo[i] = temp[i];
        }
    }
    if (pp.queryarr("rho_hi",temp,0,AMREX_SPACEDIM)) {
        for (int i=0; i<AMREX_SPACEDIM; ++i) {
            rho_hi[i] = temp[i];
        }
    }
    if (pp.queryarr("n_lo",temp,0,AMREX_SPACEDIM)) {
        for (int i=0; i<AMREX_SPACEDIM; ++i) {
            n_lo[i] = temp[i];
        }
    }
    if (pp.queryarr("n_hi",temp,0,AMREX_SPACEDIM)) {
        for (int i=0; i<AMREX_SPACEDIM; ++i) {
            n_hi[i] = temp[i];
        }
    }
    if (pp.queryarr("Yk0",temp,0,nspecies)) {
        for (int i=0; i<nspecies; ++i) {
            Yk0[i] = temp[i];
        }
    }
    if (pp.queryarr("bc_Yk_x_lo",temp,0,nspecies)) {
        for (int i=0; i<nspecies; ++i) {
            bc_Yk_x_lo[i] = temp[i];
        }
    }
    if (pp.queryarr("bc_Yk_x_hi",temp,0,nspecies)) {
        for (int i=0; i<nspecies; ++i) {
            bc_Yk_x_hi[i] = temp[i];
        }
    }
    if (pp.queryarr("bc_Yk_y_lo",temp,0,nspecies)) {
        for (int i=0; i<nspecies; ++i) {
            bc_Yk_y_lo[i] = temp[i];
        }
    }
    if (pp.queryarr("bc_Yk_y_hi",temp,0,nspecies)) {
        for (int i=0; i<nspecies; ++i) {
            bc_Yk_y_hi[i] = temp[i];
        }
    }
    if (pp.queryarr("bc_Yk_z_lo",temp,0,nspecies)) {
        for (int i=0; i<nspecies; ++i) {
            bc_Yk_z_lo[i] = temp[i];
        }
    }
    if (pp.queryarr("bc_Yk_z_hi",temp,0,nspecies)) {
        for (int i=0; i<nspecies; ++i) {
            bc_Yk_z_hi[i] = temp[i];
        }
    }
    if (pp.queryarr("bc_Xk_x_lo",temp,0,nspecies)) {
        for (int i=0; i<nspecies; ++i) {
            bc_Xk_x_lo[i] = temp[i];
        }
    }
    if (pp.queryarr("bc_Xk_x_hi",temp,0,nspecies)) {
        for (int i=0; i<nspecies; ++i) {
            bc_Xk_x_hi[i] = temp[i];
        }
    }
    if (pp.queryarr("bc_Xk_y_lo",temp,0,nspecies)) {
        for (int i=0; i<nspecies; ++i) {
            bc_Xk_y_lo[i] = temp[i];
        }
    }
    if (pp.queryarr("bc_Xk_y_hi",temp,0,nspecies)) {
        for (int i=0; i<nspecies; ++i) {
            bc_Xk_y_hi[i] = temp[i];
        }
    }
    if (pp.queryarr("bc_Xk_z_lo",temp,0,nspecies)) {
        for (int i=0; i<nspecies; ++i) {
            bc_Xk_z_lo[i] = temp[i];
        }
    }
    if (pp.queryarr("bc_Xk_z_hi",temp,0,nspecies)) {
        for (int i=0; i<nspecies; ++i) {
            bc_Xk_z_hi[i] = temp[i];
        }
    }

    if (pp.queryarr("wallspeed_x_lo",temp,0,AMREX_SPACEDIM)) {
        for (int i=0; i<AMREX_SPACEDIM; ++i) {
            wallspeed_x_lo[i] = temp[i];
        }
    }
    if (pp.queryarr("wallspeed_y_lo",temp,0,AMREX_SPACEDIM)) {
        for (int i=0; i<AMREX_SPACEDIM; ++i) {
            wallspeed_y_lo[i] = temp[i];
        }
    }
    if (pp.queryarr("wallspeed_z_lo",temp,0,AMREX_SPACEDIM)) {
        for (int i=0; i<AMREX_SPACEDIM; ++i) {
            wallspeed_z_lo[i] = temp[i];
        }
    }
    if (pp.queryarr("wallspeed_x_hi",temp,0,AMREX_SPACEDIM)) {
        for (int i=0; i<AMREX_SPACEDIM; ++i) {
            wallspeed_x_hi[i] = temp[i];
        }
    }
    if (pp.queryarr("wallspeed_y_hi",temp,0,AMREX_SPACEDIM)) {
        for (int i=0; i<AMREX_SPACEDIM; ++i) {
            wallspeed_y_hi[i] = temp[i];
        }
    }
    if (pp.queryarr("wallspeed_z_hi",temp,0,AMREX_SPACEDIM)) {
        for (int i=0; i<AMREX_SPACEDIM; ++i) {
            wallspeed_z_hi[i] = temp[i];
        }
    }

    pp.query("bc_rhotot_x_lo",bc_rhotot_x_lo);
    pp.query("bc_rhotot_x_lo",bc_rhotot_x_hi);
    pp.query("bc_rhotot_y_lo",bc_rhotot_y_lo);
    pp.query("bc_rhotot_y_lo",bc_rhotot_y_hi);
    pp.query("bc_rhotot_z_lo",bc_rhotot_z_lo);
    pp.query("bc_rhotot_z_lo",bc_rhotot_z_hi);
    pp.queryarr("wallspeed_lo",wallspeed_lo,0,(AMREX_SPACEDIM-1)*AMREX_SPACEDIM);
    pp.queryarr("wallspeed_hi",wallspeed_hi,0,(AMREX_SPACEDIM-1)*AMREX_SPACEDIM);
    if (pp.queryarr("potential_lo",temp,0,AMREX_SPACEDIM)) {
        for (int i=0; i<AMREX_SPACEDIM; ++i) {
            potential_lo[i] = temp[i];
        }
    }
    if (pp.queryarr("potential_hi",temp,0,AMREX_SPACEDIM)) {
        for (int i=0; i<AMREX_SPACEDIM; ++i) {
            potential_hi[i] = temp[i];
        }
    }
    pp.query("dsmc_boundaries",dsmc_boundaries);
    pp.query("n_burn",n_burn);
    pp.query("phonon_sound_speed",phonon_sound_speed);
    pp.query("tau_i",tau_i);
    pp.query("tau_ta",tau_ta);
    pp.query("tau_la",tau_la);
    pp.query("toggleTimeFrac",toggleTimeFrac);
    pp.query("struct_fact_int",struct_fact_int);
    pp.query("radialdist_int",radialdist_int);
    pp.query("cartdist_int",cartdist_int);
    pp.query("n_steps_skip",n_steps_skip);
    pp.query("binSize",binSize);
    pp.query("searchDist",searchDist);
    pp.query("project_dir",project_dir);
    pp.query("slicepoint",slicepoint);
    pp.queryarr("max_grid_projection",max_grid_projection,0,AMREX_SPACEDIM-1);
    pp.query("histogram_unit",histogram_unit);
    pp.queryarr("density_weights",density_weights,0,nspecies);
    pp.queryarr("shift_cc_to_boundary",shift_cc_to_boundary,0,AMREX_SPACEDIM*2);
    pp.query("particle_placement",particle_placement);
    pp.query("particle_input",particle_input);
    pp.queryarr("particle_count",particle_count,0,nspecies);
    pp.queryarr("p_move_tog",p_move_tog,0,nspecies);
    pp.queryarr("p_force_tog",p_force_tog,0,nspecies);
    pp.query("particle_neff",particle_neff);
    pp.queryarr("particle_n0",particle_n0,0,nspecies);
    pp.queryarr("mass",mass,0,nspecies);
    pp.queryarr("nfrac",nfrac,0,nspecies);
    if (pp.queryarr("p_int_tog",temp_int,0,nspecies*nspecies)) {
        for (int i=0; i<nspecies*nspecies; ++i) {
            p_int_tog[i] = temp_int[i];
        }
    }
    if (pp.queryarr("eepsilon",temp,0,nspecies*nspecies)) {
        for (int i=0; i<nspecies*nspecies; ++i) {
            eepsilon[i] = temp[i];
        }
    }
    if (pp.queryarr("sigma",temp,0,nspecies*nspecies)) {
        for (int i=0; i<nspecies*nspecies; ++i) {
            sigma[i] = temp[i];
        }
    }
    if (pp.queryarr("rmin",temp,0,nspecies*nspecies)) {
        for (int i=0; i<nspecies*nspecies; ++i) {
            rmin[i] = temp[i];
        }
    }
    if (pp.queryarr("rmax",temp,0,nspecies*nspecies)) {
        for (int i=0; i<nspecies*nspecies; ++i) {
            rmax[i] = temp[i];
        }
    }
    if (pp.queryarr("alpha_pp",temp,0,nspecies*nspecies)) {
        for (int i=0; i<nspecies*nspecies; ++i) {
            alpha_pp[i] = temp[i];
        }
    }
    if (pp.queryarr("alpha_pw",temp,0,nspecies*nspecies)) {
        for (int i=0; i<nspecies*nspecies; ++i) {
            alpha_pw[i] = temp[i];
        }
    }
    if (pp.queryarr("friction_pp",temp,0,nspecies*nspecies)) {
        for (int i=0; i<nspecies*nspecies; ++i) {
            friction_pp[i] = temp[i];
        }
    }
    if (pp.queryarr("friction_pw",temp,0,nspecies*nspecies)) {
        for (int i=0; i<nspecies*nspecies; ++i) {
            friction_pw[i] = temp[i];
        }
    }
    if (pp.queryarr("phi_domain",temp,0,nspecies)) {
        for (int i=0; i<nspecies; ++i) {
            phi_domain[i] = temp[i];
        }
    }
    if (pp.queryarr("p_int_tog_wall",temp_int,0,nspecies)) {
        for (int i=0; i<nspecies; ++i) {
            p_int_tog_wall[i] = temp_int[i];
        }
    }
    if (pp.queryarr("msd_int",temp_int,0,nspecies)) {
        for (int i=0; i<nspecies; ++i) {
            msd_int[i] = temp_int[i];
        }
    }
    if (pp.queryarr("msd_len",temp_int,0,nspecies)) {
        for (int i=0; i<nspecies; ++i) {
            msd_len[i] = temp_int[i];
        }
    }
    if (pp.queryarr("eepsilon_wall",temp,0,nspecies)) {
        for (int i=0; i<nspecies; ++i) {
            eepsilon_wall[i] = temp[i];
        }
    }
    if (pp.queryarr("sigma_wall",temp,0,nspecies)) {
        for (int i=0; i<nspecies; ++i) {
            sigma_wall[i] = temp[i];
        }
    }
    if (pp.queryarr("rmin_wall",temp,0,nspecies)) {
        for (int i=0; i<nspecies; ++i) {
            rmin_wall[i] = temp[i];
        }
    }
    if (pp.queryarr("offset_wall",temp,0,nspecies)) {
        for (int i=0; i<nspecies; ++i) {
            offset_wall[i] = temp[i];
        }
    }
    if (pp.queryarr("rmax_wall",temp,0,nspecies)) {
        for (int i=0; i<nspecies; ++i) {
            rmax_wall[i] = temp[i];
        }
    }
    pp.query("poisson_verbose",poisson_verbose);
    pp.query("poisson_bottom_verbose",poisson_bottom_verbose);
    pp.query("poisson_max_iter",poisson_max_iter);
    pp.query("poisson_rel_tol",poisson_rel_tol);
    pp.query("permittivity",permittivity);
    pp.query("wall_mob",wall_mob);
    pp.query("particle_grid_refine",particle_grid_refine);
    pp.query("es_grid_refine",es_grid_refine);
    pp.queryarr("diff",diff,0,nspecies);
    pp.query("all_dry",all_dry);
    pp.query("fluid_tog",fluid_tog);
    pp.query("es_tog",es_tog);
    pp.query("drag_tog",drag_tog);
    pp.query("move_tog",move_tog);
    pp.query("rfd_tog",rfd_tog);
    pp.query("dry_move_tog",dry_move_tog);
    pp.query("sr_tog",sr_tog);
    pp.query("graphene_tog",graphene_tog);
    pp.query("thermostat_tog",thermostat_tog);
    pp.query("zero_net_force",zero_net_force);
    pp.query("crange",crange);
    pp.query("images",images);
    pp.queryarr("eamp",eamp,0,3);
    pp.queryarr("efreq",efreq,0,3);
    pp.queryarr("ephase",ephase,0,3);
    pp.queryarr("body_force_density",body_force_density,0,3);
    pp.query("plot_ascii",plot_ascii);
    pp.query("plot_means",plot_means);
    pp.query("plot_vars",plot_vars);
    pp.query("plot_mom3",plot_mom3);
    pp.query("plot_covars",plot_covars);
    pp.query("plot_cross",plot_cross);
    pp.query("plot_deltaY_dir",plot_deltaY_dir);
    pp.query("particle_motion",particle_motion);
    pp.query("turb_a",turb_a);
    pp.query("turb_b",turb_b);
    pp.query("turb_c",turb_c);
    pp.query("turb_d",turb_d);
    pp.query("turb_alpha",turb_alpha);
    pp.query("turbForcing",turbForcing);

    if (nspecies > MAX_SPECIES) {
        Abort("InitializeCommonNamespace: nspecies > MAX_SPECIES; re-compile with a larger MAX_SPEC as a compiler input");
    }

    if (wallspeed_x_lo[0] != 0.) {
        Abort("you are specifying a normal velocity on a wall; wallspeed_x_lo[0] must be 0");
    }
    if (wallspeed_x_hi[0] != 0.) {
        Abort("you are specifying a normal velocity on a wall; wallspeed_x_hi[0] must be 0");
    }
    if (wallspeed_y_lo[1] != 0.) {
        Abort("you are specifying a normal velocity on a wall; wallspeed_y_lo[1] must be 0");
    }
    if (wallspeed_y_hi[1] != 0.) {
        Abort("you are specifying a normal velocity on a wall; wallspeed_y_hi[1] must be 0");
    }
#if (AMREX_SPACEDIM == 3)
    if (wallspeed_z_lo[2] != 0.) {
        Abort("you are specifying a normal velocity on a wall; wallspeed_z_lo[2] must be 0");
    }
    if (wallspeed_z_hi[2] != 0.) {
        Abort("you are specifying a normal velocity on a wall; wallspeed_z_hi[2] must be 0");
    }
#endif


}