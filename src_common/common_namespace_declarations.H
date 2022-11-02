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
amrex::Vector<int>         common::max_grid_size_structfact;
amrex::Vector<int>         common::max_particle_tile_size;
amrex::Real                common::cell_depth;
amrex::IntVect             common::ngc;
AMREX_GPU_MANAGED int      common::nvars;
AMREX_GPU_MANAGED int      common::nprimvars;

int                        common::membrane_cell;
int                        common::cross_cell;
amrex::Real                common::transmission;

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
int                        common::particle_restart;
int                        common::print_int;
int                        common::project_eos_int;
AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> common::grav;
AMREX_GPU_MANAGED int      common::nspecies;
AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> common::molmass;
AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> common::rhobar;
AMREX_GPU_MANAGED amrex::Real common::rho0;
AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> common::diameter;
amrex::Vector<int>         common::dof;
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
amrex::Real                common::smoothing_width;
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

AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> common::bc_Yk_x_lo;
AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> common::bc_Yk_x_hi;
AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> common::bc_Yk_y_lo;
AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> common::bc_Yk_y_hi;
AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> common::bc_Yk_z_lo;
AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> common::bc_Yk_z_hi;

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

int                        common::struct_fact_int;
int                        common::radialdist_int;
int                        common::cartdist_int;
int                        common::n_steps_skip;
AMREX_GPU_MANAGED amrex::Real common::binSize;
AMREX_GPU_MANAGED amrex::Real common::searchDist;
int                        common::project_dir;
int                        common::slicepoint;
amrex::Vector<int>         common::max_grid_projection;
int                        common::histogram_unit;
amrex::Vector<amrex::Real> common::density_weights;
amrex::Vector<int>         common::shift_cc_to_boundary;

int                        common::particle_placement;
amrex::Vector<int>         common::particle_count;
amrex::Vector<int>         common::p_move_tog;
amrex::Vector<int>         common::p_force_tog;

amrex::Real                common::particle_neff;
amrex::Vector<amrex::Real> common::particle_n0;
amrex::Vector<amrex::Real> common::mass;
amrex::Vector<amrex::Real> common::nfrac;

AMREX_GPU_MANAGED amrex::GpuArray<int, MAX_SPECIES*MAX_SPECIES> common::p_int_tog;
AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES*MAX_SPECIES> common::eepsilon;
AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES*MAX_SPECIES> common::sigma;
AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES*MAX_SPECIES> common::rmin;
AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES*MAX_SPECIES> common::rmax;

AMREX_GPU_MANAGED amrex::GpuArray<int, MAX_SPECIES> common::p_int_tog_wall;
AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> common::eepsilon_wall;
AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> common::sigma_wall;
AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> common::rmin_wall;
AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> common::rmax_wall;

int                        common::poisson_verbose;
int                        common::poisson_bottom_verbose;
int                        common::poisson_max_iter;

amrex::Real                common::poisson_rel_tol;
AMREX_GPU_MANAGED amrex::Real common::permittivity;
int                        common::induced_charge_eo;
AMREX_GPU_MANAGED int         common::zero_eps_on_wall_type;
AMREX_GPU_MANAGED amrex::Real common::zero_eps_wall_left_end;
AMREX_GPU_MANAGED amrex::Real common::zero_eps_wall_right_start;
AMREX_GPU_MANAGED int         common::wall_mob;

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
AMREX_GPU_MANAGED int	   common::body_force;

int                        common::crange;

AMREX_GPU_MANAGED int      common::images;
AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> common::eamp;
AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> common::efreq;
AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> common::ephase;
AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> common::F_body;

int                        common::plot_ascii;
int                        common::plot_means;
int                        common::plot_vars;
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
