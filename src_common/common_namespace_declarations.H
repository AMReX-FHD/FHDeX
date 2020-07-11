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
amrex::Vector<int>         common::n_cells;
amrex::Vector<int>         common::max_grid_size;
amrex::Vector<int>         common::max_particle_tile_size;
amrex::Real                common::cell_depth;
amrex::IntVect             common::ngc;
int                        common::nvars;
int                        common::nprimvars;

int                        common::membrane_cell;
int                        common::cross_cell;
amrex::Real                common::transmission;

AMREX_GPU_MANAGED int      common::pkernel_fluid;
AMREX_GPU_MANAGED int      common::pkernel_es;
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
int                        common::prob_type;
int                        common::restart;
int                        common::particle_restart;
int                        common::print_int;
int                        common::project_eos_int;
amrex::Vector<amrex::Real> common::grav;
AMREX_GPU_MANAGED int      common::nspecies;
amrex::Vector<amrex::Real> common::molmass;
amrex::Vector<amrex::Real> common::rhobar;
amrex::Vector<amrex::Real> common::diameter;
amrex::Vector<int>         common::dof;
amrex::Vector<amrex::Real> common::hcv;
amrex::Vector<amrex::Real> common::hcp;


amrex::Real                common::rho0;
amrex::Real                common::variance_coef_mom;
amrex::Real                common::variance_coef_mass;
amrex::Real                common::k_B;
amrex::Real                common::Runiv;
amrex::Vector<Real>        common::T_init;
int                        common::algorithm_type;
int                        common::barodiffusion_type;
int                        common::use_bl_rng;
int                        common::seed;
int                        common::seed_momentum;
int                        common::seed_diffusion;
int                        common::seed_reaction;
int                        common::seed_init_mass;
int                        common::seed_init_momentum;
amrex::Real                common::visc_coef;
int                        common::visc_type;
int                        common::advection_type;
int                        common::filtering_width;
int                        common::stoch_stress_form;
amrex::Vector<amrex::Real> common::u_init;
amrex::Real                common::perturb_width;
amrex::Real                common::smoothing_width;
amrex::Real                common::initial_variance_mom;
amrex::Real                common::initial_variance_mass;
amrex::Real                common::domega;


amrex::Vector<int>         common::bc_vel_lo;
amrex::Vector<int>         common::bc_vel_hi;
amrex::Vector<int>         common::bc_es_lo;
amrex::Vector<int>         common::bc_es_hi;
amrex::Vector<int>         common::bc_mass_lo;
amrex::Vector<int>         common::bc_mass_hi;
amrex::Vector<int>         common::bc_therm_lo;
amrex::Vector<int>         common::bc_therm_hi;

amrex::Vector<amrex::Real> common::p_lo;
amrex::Vector<amrex::Real> common::p_hi;

amrex::Vector<amrex::Real> common::t_lo;
amrex::Vector<amrex::Real> common::t_hi;

amrex::Vector<amrex::Real> common::wallspeed_lo;
amrex::Vector<amrex::Real> common::wallspeed_hi;

amrex::Vector<amrex::Real> common::potential_lo;
amrex::Vector<amrex::Real> common::potential_hi;

int                        common::struct_fact_int;
int                        common::radialdist_int;
int                        common::cartdist_int;
int                        common::n_steps_skip;
AMREX_GPU_MANAGED amrex::Real common::binSize;
AMREX_GPU_MANAGED amrex::Real common::searchDist;
int                        common::project_dir;
amrex::Vector<int>         common::max_grid_projection;
int                        common::histogram_unit;
amrex::Vector<amrex::Real> common::density_weights;
amrex::Vector<int>         common::shift_cc_to_boundary;

int                        common::particle_placement;
amrex::Vector<int>         common::particle_count;
amrex::Vector<int>         common::p_move_tog;
amrex::Vector<int>         common::p_force_tog;
amrex::Vector<int>         common::p_int_tog;
amrex::Real                common::particle_neff;
amrex::Vector<amrex::Real> common::particle_n0;
amrex::Vector<amrex::Real> common::mass;
amrex::Vector<amrex::Real> common::nfrac;

amrex::Vector<amrex::Real> common::eepsilon;
amrex::Vector<amrex::Real> common::sigma;

int                        common::poisson_verbose;
int                        common::poisson_bottom_verbose;
int                        common::poisson_max_iter;

amrex::Real                common::poisson_rel_tol;
AMREX_GPU_MANAGED amrex::Real common::permittivity;
amrex::Real                common::cut_off;
amrex::Real                common::rmin;

amrex::Real                common::particle_grid_refine;
amrex::Real                common::es_grid_refine;
amrex::Vector<amrex::Real> common::diff;
int                        common::all_dry;

int                        common::fluid_tog;
int                        common::es_tog;
int                        common::drag_tog;
int                        common::move_tog;
int                        common::rfd_tog;
int                        common::dry_move_tog;
int                        common::sr_tog;
int                        common::graphene_tog;
int	                   common::thermostat_tog;
int	                   common::zero_net_force;

int                        common::crange;

AMREX_GPU_MANAGED int      common::images;
amrex::Vector<amrex::Real> common::eamp;
amrex::Vector<amrex::Real> common::efreq;
amrex::Vector<amrex::Real> common::ephase;

int                        common::plot_ascii;
int                        common::particle_motion;
int                        common::solve_chem;
amrex::Real                common::diffcoeff;
amrex::Real                common::scaling_factor;
amrex::Real                common::source_strength;
int                        common::regrid_int;
int                        common::do_reflux;

