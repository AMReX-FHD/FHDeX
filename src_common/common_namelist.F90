module common_namelist_module

  use iso_c_binding, only: c_char
  use amrex_string_module, only: amrex_string_c_to_f, amrex_string_f_to_c

  implicit none

  integer, parameter :: MAX_SPECIES = 10
  integer, parameter :: LOHI = 2

  double precision,   save :: prob_lo(AMREX_SPACEDIM)
  double precision,   save :: prob_hi(AMREX_SPACEDIM)
  integer,            save :: n_cells(AMREX_SPACEDIM)
  integer,            save :: max_grid_size(AMREX_SPACEDIM)
  double precision,   save :: fixed_dt
  double precision,   save :: cfl
  integer,            save :: max_step
  integer,            save :: plot_int
  character(len=128), save :: plot_base_name
  integer,            save :: chk_int
  character(len=128), save :: chk_base_name
  integer,            save :: prob_type
  integer,            save :: restart
  integer,            save :: print_int
  integer,            save :: project_eos_int
  double precision,   save :: grav(AMREX_SPACEDIM)
  integer,            save :: nspecies
  double precision,   save :: molmass(MAX_SPECIES)
  double precision,   save :: rhobar(MAX_SPECIES)
  double precision,   save :: rho0
  double precision,   save :: variance_coef_mom
  double precision,   save :: variance_coef_mass
  double precision,   save :: k_B
  double precision,   save :: Runiv
  integer,            save :: algorithm_type
  integer,            save :: advection_type
  integer,            save :: barodiffusion_type
  integer,            save :: use_bl_rng
  integer,            save :: seed
  integer,            save :: seed_momentum
  integer,            save :: seed_diffusion
  integer,            save :: seed_reaction
  integer,            save :: seed_init_mass
  integer,            save :: seed_init_momentum
  double precision,   save :: visc_coef
  integer,            save :: visc_type
  integer,            save :: filtering_width
  integer,            save :: stoch_stress_form
  double precision,   save :: u_init(2)
  double precision,   save :: perturb_width
  double precision,   save :: smoothing_width
  double precision,   save :: initial_variance_mom
  double precision,   save :: initial_variance_mass
  integer,            save :: bc_lo(AMREX_SPACEDIM)
  integer,            save :: bc_hi(AMREX_SPACEDIM)
  double precision,   save :: wallspeed_lo(AMREX_SPACEDIM-1,AMREX_SPACEDIM)
  double precision,   save :: wallspeed_hi(AMREX_SPACEDIM-1,AMREX_SPACEDIM)
  integer,            save :: histogram_unit
  double precision,   save :: density_weights(MAX_SPECIES)
  integer,            save :: shift_cc_to_boundary(AMREX_SPACEDIM,LOHI)
  
  ! Problem specification
  namelist /common/ prob_lo       ! physical lo coordinate
  namelist /common/ prob_hi       ! physical hi coordinate
  namelist /common/ n_cells       ! number of cells in domain
  namelist /common/ max_grid_size ! max number of cells in a box

  ! Time-step control
  namelist /common/ fixed_dt
  namelist /common/ cfl

  ! Controls for number of steps between actions
  namelist /common/ max_step
  namelist /common/ plot_int
  namelist /common/ plot_base_name
  namelist /common/ chk_int
  namelist /common/ chk_base_name
  namelist /common/ prob_type
  namelist /common/ restart
  namelist /common/ print_int
  namelist /common/ project_eos_int

  ! Physical parameters
  namelist /common/ grav
  namelist /common/ nspecies
  namelist /common/ molmass
  namelist /common/ rhobar
  namelist /common/ rho0

  ! stochastic forcing amplitudes (1 for physical values, 0 to run them off)
  namelist /common/ variance_coef_mom
  namelist /common/ variance_coef_mass
  namelist /common/ k_B
  namelist /common/ Runiv

  ! Algorithm control / selection
  namelist /common/ algorithm_type
  namelist /common/ advection_type
  namelist /common/ barodiffusion_type
  namelist /common/ use_bl_rng

  ! random number seed (for HydroGrid RNGs)
  ! 0        = unpredictable seed based on clock
  ! positive = fixed seed
  namelist /common/ seed

  ! Random number seeds for each physical process for use_bl_rng=T
  ! for positive value, the value is assigned as seed value
  ! for 0, a positive value is randomly chosen
  ! if -1 (only for restart), RNGs status is restored from checkpoint data
  namelist /common/ seed_momentum
  namelist /common/ seed_diffusion
  namelist /common/ seed_reaction
  namelist /common/ seed_init_mass
  namelist /common/ seed_init_momentum

  ! Viscous friction L phi operator
  ! if abs(visc_type) = 1, L = div beta grad
  ! if abs(visc_type) = 2, L = div [ beta (grad + grad^T) ]
  ! if abs(visc_type) = 3, L = div [ beta (grad + grad^T) + I (gamma - (2/3)*beta) div ]
  ! positive = assume constant coefficients
  ! negative = assume spatially-varying coefficients
  namelist /common/ visc_coef
  namelist /common/ visc_type

  ! Stochastic momentum flux controls:
  namelist /common/ filtering_width
  namelist /common/ stoch_stress_form

  ! Initial conditions
  namelist /common/ u_init
  namelist /common/ perturb_width
  namelist /common/ smoothing_width
  namelist /common/ initial_variance_mom
  namelist /common/ initial_variance_mass

  ! Boundary conditions
  ! ----------------------
  ! BC specifications:
  ! -1 = periodic
  ! 100 = no-slip wall      (Dir condition for normal vel; Dir velocity condition for trans vel)
  ! 101 = no-slip reservoir (Dir condition for normal vel; Dir velocity condition for trans vel)
  ! 200 = slip wall         (Dir condition for normal vel; Dir traction condition for trans vel)
  ! 201 = slip reservoir    (Dir condition for normal vel; Dir traction condition for trans vel)
  ! For a complete list see ???
  namelist /common/ bc_lo
  namelist /common/ bc_hi

  ! Each no-slip wall may be moving with a specified tangential 
  ! velocity along the tangential directions
  ! In 2D:
  ! wallspeed_lo/hi_x - yvel
  ! wallspeed_lo/hi_y - xvel
  ! In 3D:
  ! wallspeed_lo/hi_x - yvel,zvel
  ! wallspeed_lo/hi_y - xvel,zvel
  ! wallspeed_lo/hi_z - xvel,yvel
  namelist /common/ wallspeed_lo
  namelist /common/ wallspeed_hi

  ! These are mostly used for reaction-diffusion: 
  namelist /common/ histogram_unit
  namelist /common/ density_weights
  namelist /common/ shift_cc_to_boundary

contains

  ! read in fortran namelist into common_params_module
  subroutine read_common_namelist(inputs_file,length) bind(C, name="read_common_namelist")

    integer               , value         :: length
    character(kind=c_char), intent(in   ) :: inputs_file(length)

    ! default values
    prob_lo(:) = 0.d0
    prob_hi(:) = 1.d0
    n_cells(:) = 1
    max_grid_size(:) = 1
    fixed_dt = 1.
    cfl = 0.5
    max_step = 1
    plot_int = 0
    plot_base_name = "plt"
    chk_int = 0
    chk_base_name = "chk"
    prob_type = 1
    restart = -1
    print_int = 0
    project_eos_int = -1
    grav(:) = 0.d0
    nspecies = 2
    molmass(:) = 1.d0
    rhobar(:) = 1.d0
    rho0 = 1.
    variance_coef_mom = 1.
    variance_coef_mass = 1.
    k_B = 1.
    Runiv = 8.314462175d7
    algorithm_type = 0
    advection_type = 0
    barodiffusion_type = 0
    use_bl_rng = 0
    seed = 1
    seed_momentum = 1
    seed_diffusion = 1
    seed_reaction = 1
    seed_init_mass = 1
    seed_init_momentum = 1
    visc_coef = 1.
    visc_type = 1
    filtering_width = 0
    stoch_stress_form = 1
    u_init(:) = 0.d0
    perturb_width = 0.
    smoothing_width = 1.
    initial_variance_mom = 0.
    initial_variance_mass = 0.
    bc_lo(:) = 0
    bc_hi(:) = 0
    wallspeed_lo(:,:) = 0
    wallspeed_hi(:,:) = 0
    histogram_unit = 0
    density_weights(:) = 0.d0
    shift_cc_to_boundary(:,:) = 0

    ! read in common namelist
    open(unit=100, file=amrex_string_c_to_f(inputs_file), status='old', action='read')
    read(unit=100, nml=common)
    close(unit=100)

  end subroutine read_common_namelist

  ! copy contents of common_params_module to C++ common namespace
  subroutine initialize_common_namespace(prob_lo_in, prob_hi_in, n_cells_in, &
                                         max_grid_size_in, &
                                         fixed_dt_in, cfl_in, max_step_in, plot_int_in, &
                                         plot_base_name_in, plot_base_name_len, chk_int_in, &
                                         chk_base_name_in, chk_base_name_len, prob_type_in, &
                                         restart_in, print_int_in, project_eos_int_in, &
                                         grav_in, nspecies_in, molmass_in, rhobar_in, &
                                         rho0_in, variance_coef_mom_in, &
                                         variance_coef_mass_in, &
                                         k_B_in, Runiv_in, algorithm_type_in, & 
                                         advection_type_in, &
                                         barodiffusion_type_in, use_bl_rng_in, seed_in, &
                                         seed_momentum_in, seed_diffusion_in, &
                                         seed_reaction_in, &
                                         seed_init_mass_in, seed_init_momentum_in, &
                                         visc_coef_in, visc_type_in, &
                                         filtering_width_in, stoch_stress_form_in, &
                                         u_init_in, perturb_width_in, smoothing_width_in, &
                                         initial_variance_mom_in, initial_variance_mass_in, &
                                         bc_lo_in, bc_hi_in, &
                                         wallspeed_lo_in, wallspeed_hi_in, &
                                         histogram_unit_in, density_weights_in, &
                                         shift_cc_to_boundary_in) &
                                         bind(C, name="initialize_common_namespace")

    double precision,       intent(inout) :: prob_lo_in(AMREX_SPACEDIM)
    double precision,       intent(inout) :: prob_hi_in(AMREX_SPACEDIM)
    integer,                intent(inout) :: n_cells_in(AMREX_SPACEDIM)
    integer,                intent(inout) :: max_grid_size_in(AMREX_SPACEDIM)
    double precision,       intent(inout) :: fixed_dt_in
    double precision,       intent(inout) :: cfl_in
    integer,                intent(inout) :: max_step_in
    integer,                intent(inout) :: plot_int_in
    integer               , value         :: plot_base_name_len
    character(kind=c_char), intent(inout) :: plot_base_name_in(plot_base_name_len)
    integer,                intent(inout) :: chk_int_in
    integer               , value         :: chk_base_name_len
    character(kind=c_char), intent(inout) :: chk_base_name_in(chk_base_name_len)
    integer,                intent(inout) :: prob_type_in
    integer,                intent(inout) :: restart_in
    integer,                intent(inout) :: print_int_in
    integer,                intent(inout) :: project_eos_int_in
    double precision,       intent(inout) :: grav_in(AMREX_SPACEDIM)
    integer,                intent(inout) :: nspecies_in
    double precision,       intent(inout) :: molmass_in(MAX_SPECIES)
    double precision,       intent(inout) :: rhobar_in(MAX_SPECIES)
    double precision,       intent(inout) :: rho0_in
    double precision,       intent(inout) :: variance_coef_mom_in
    double precision,       intent(inout) :: variance_coef_mass_in
    double precision,       intent(inout) :: k_B_in
    double precision,       intent(inout) :: Runiv_in
    integer,                intent(inout) :: algorithm_type_in
    integer,                intent(inout) :: advection_type_in
    integer,                intent(inout) :: barodiffusion_type_in
    integer,                intent(inout) :: use_bl_rng_in
    integer,                intent(inout) :: seed_in
    integer,                intent(inout) :: seed_momentum_in
    integer,                intent(inout) :: seed_diffusion_in
    integer,                intent(inout) :: seed_reaction_in
    integer,                intent(inout) :: seed_init_mass_in
    integer,                intent(inout) :: seed_init_momentum_in
    double precision,       intent(inout) :: visc_coef_in
    integer,                intent(inout) :: visc_type_in
    integer,                intent(inout) :: filtering_width_in
    integer,                intent(inout) :: stoch_stress_form_in
    double precision,       intent(inout) :: u_init_in(2)
    double precision,       intent(inout) :: perturb_width_in
    double precision,       intent(inout) :: smoothing_width_in
    double precision,       intent(inout) :: initial_variance_mom_in
    double precision,       intent(inout) :: initial_variance_mass_in
    integer,                intent(inout) :: bc_lo_in(AMREX_SPACEDIM)
    integer,                intent(inout) :: bc_hi_in(AMREX_SPACEDIM)
    double precision,       intent(inout) :: wallspeed_lo_in(AMREX_SPACEDIM-1,AMREX_SPACEDIM)
    double precision,       intent(inout) :: wallspeed_hi_in(AMREX_SPACEDIM-1,AMREX_SPACEDIM)
    integer,                intent(inout) :: histogram_unit_in
    double precision,       intent(inout) :: density_weights_in(MAX_SPECIES)
    integer,                intent(inout) :: shift_cc_to_boundary_in(AMREX_SPACEDIM,LOHI)

    prob_lo_in = prob_lo
    prob_hi_in = prob_hi
    n_cells_in = n_cells
    max_grid_size_in = max_grid_size
    fixed_dt_in = fixed_dt
    cfl_in = cfl
    max_step_in = max_step
    plot_int_in = plot_int
    plot_base_name_in = amrex_string_f_to_c(plot_base_name)
    chk_int_in = chk_int
    chk_base_name_in = amrex_string_f_to_c(chk_base_name)
    prob_type_in = prob_type
    restart_in = restart
    print_int_in = print_int
    project_eos_int_in = project_eos_int
    grav_in = grav
    nspecies_in = nspecies
    molmass_in = molmass
    rhobar_in = rhobar
    rho0_in= rho0
    variance_coef_mom_in = variance_coef_mom
    variance_coef_mass_in = variance_coef_mass
    k_B_in = k_B
    Runiv_in = Runiv
    algorithm_type_in = algorithm_type
    advection_type_in = advection_type
    barodiffusion_type_in = barodiffusion_type
    use_bl_rng_in = use_bl_rng
    seed_in = seed
    seed_momentum_in = seed_momentum
    seed_diffusion_in = seed_diffusion
    seed_reaction_in = seed_reaction
    seed_init_mass_in = seed_init_mass
    seed_init_momentum_in = seed_init_momentum
    visc_coef_in = visc_coef
    visc_type_in = visc_type
    filtering_width_in = filtering_width
    stoch_stress_form_in = stoch_stress_form
    u_init_in = u_init
    perturb_width_in = perturb_width
    smoothing_width_in = smoothing_width
    initial_variance_mom_in = initial_variance_mom
    initial_variance_mass_in = initial_variance_mass
    bc_lo_in = bc_lo
    bc_hi_in = bc_hi
    wallspeed_lo_in = wallspeed_lo
    wallspeed_hi_in = wallspeed_hi
    histogram_unit_in = histogram_unit
    density_weights_in = density_weights
    shift_cc_to_boundary_in = shift_cc_to_boundary

  end subroutine initialize_common_namespace

end module common_namelist_module
