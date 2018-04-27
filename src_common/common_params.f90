module common_params_module

  implicit none

  integer, parameter :: MAX_SPECIES = 6
  integer, parameter :: MAX_SPACEDIM = 3

  ! Begin the declarations of the ParmParse parameters

  double precision,  allocatable, save :: prob_lo(:)
  double precision,  allocatable, save :: prob_hi(:)
  integer,           allocatable, save :: n_cells(:)
  integer,           allocatable, save :: max_grid_size(:)
  double precision,               save :: fixed_dt
  double precision,               save :: cfl
  integer,                        save :: max_step
  integer,                        save :: plot_int
  character (len=:), allocatable, save :: plot_base_name
  integer,                        save :: chk_int
  character (len=:), allocatable, save :: chk_base_name
  integer,                        save :: prob_type
  integer,                        save :: restart
  integer,                        save :: print_int
  integer,                        save :: project_eos_int
  double precision,  allocatable, save :: grav(:)
  integer,                        save :: nspecies
  double precision,  allocatable, save :: molmass(:)
  double precision,  allocatable, save :: rhobar(:)
  double precision,               save :: rho0
  double precision,               save :: variance_coef_mom
  double precision,               save :: variance_coef_mass
  double precision,               save :: k_B
  double precision,               save :: Runiv
  integer,                        save :: algorithm_type
  integer,                        save :: barodiffusion_type
  logical,                        save :: use_bl_rng
  integer,                        save :: seed
  integer,                        save :: seed_momentum
  integer,                        save :: seed_diffusion
  integer,                        save :: seed_reaction
  integer,                        save :: seed_init_mass
  integer,                        save :: seed_init_momentum
  double precision,               save :: visc_coef
  integer,                        save :: visc_type
  integer,                        save :: advection_type
  integer,                        save :: filtering_width
  integer,                        save :: stoch_stress_form
  double precision,  allocatable, save :: u_init(:)
  double precision,               save :: perturb_width
  double precision,               save :: smoothing_width
  double precision,               save :: initial_variance_mom
  double precision,               save :: initial_variance_mass
  integer,           allocatable, save :: bc_lo(:)
  integer,           allocatable, save :: bc_hi(:)
  double precision,  allocatable, save :: wallspeed_lo_x(:)
  double precision,  allocatable, save :: wallspeed_hi_x(:)
  double precision,  allocatable, save :: wallspeed_lo_y(:)
  double precision,  allocatable, save :: wallspeed_hi_y(:)
  double precision,  allocatable, save :: wallspeed_lo_z(:)
  double precision,  allocatable, save :: wallspeed_hi_z(:)
  integer,                        save :: histogram_unit
  double precision,  allocatable, save :: density_weights(:)
  integer,           allocatable, save :: shift_cc_to_boundary_lo(:)
  integer,           allocatable, save :: shift_cc_to_boundary_hi(:)

  ! End the declarations of the ParmParse parameters

contains

  subroutine read_common_params() bind(C, name="read_common_params")

    use amrex_parmparse_module, only: amrex_parmparse_build, &
                                      amrex_parmparse_destroy, &
                                      amrex_parmparse

    type (amrex_parmparse) :: pp

    ! allocate arrays
    allocate(prob_lo(MAX_SPACEDIM))
    allocate(prob_hi(MAX_SPACEDIM))
    allocate(n_cells(MAX_SPACEDIM))
    allocate(character(len=1)::plot_base_name)
    allocate(max_grid_size(MAX_SPACEDIM))
    allocate(grav(MAX_SPACEDIM))
    allocate(molmass(MAX_SPECIES))
    allocate(rhobar(MAX_SPECIES))
    allocate(u_init(2))
    allocate(bc_lo(MAX_SPACEDIM))
    allocate(bc_hi(MAX_SPACEDIM))
    allocate(wallspeed_lo_x(2))
    allocate(wallspeed_hi_x(2))
    allocate(wallspeed_lo_y(2))
    allocate(wallspeed_hi_y(2))
    allocate(wallspeed_lo_z(2))
    allocate(wallspeed_hi_z(2))
    allocate(density_weights(MAX_SPECIES))
    allocate(shift_cc_to_boundary_lo(MAX_SPACEDIM))
    allocate(shift_cc_to_boundary_hi(MAX_SPACEDIM))

    ! default values
    prob_lo(:) = 0.d0
    plot_base_name = "plt"
    fixed_dt = 1.d0

    ! read in from inputs file
    
    call amrex_parmparse_build(pp)
    call pp%queryarr("prob_lo",prob_lo)
    call pp%query("fixed_dt", fixed_dt)
    call pp%query("plot_base_name",plot_base_name)

  end subroutine read_common_params

end module common_params_module
