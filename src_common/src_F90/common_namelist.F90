#include <AMReX_Config.H>


module common_namelist_module

  use iso_c_binding, only: c_char
  use amrex_string_module, only: amrex_string_c_to_f, amrex_string_f_to_c
  use amrex_error_module
  
  implicit none

  integer, parameter :: MAX_SPECIES = 4
  integer, parameter :: LOHI = 2

  double precision,   save :: prob_lo(AMREX_SPACEDIM)
  double precision,   save :: prob_hi(AMREX_SPACEDIM)
  integer,            save :: n_cells(AMREX_SPACEDIM)
  integer,            save :: max_grid_size(AMREX_SPACEDIM)
  integer,            save :: max_grid_size_structfact(AMREX_SPACEDIM)
  integer,            save :: max_particle_tile_size(AMREX_SPACEDIM)
  double precision,   save :: cell_depth

  integer,            save :: ngc(AMREX_SPACEDIM)
  integer,            save :: nvars
  integer,            save :: nprimvars
  integer,            save :: membrane_cell
  integer,            save :: cross_cell
  double precision,   save :: transmission

  double precision,   save :: qval(MAX_SPECIES)
  integer,            save :: pkernel_fluid(MAX_SPECIES) ! GALEN - FLUID KERNEL
  integer,            save :: pkernel_es(MAX_SPECIES)
  integer,            save :: eskernel_fluid(MAX_SPECIES) ! ES KERNEL
  double precision,   save :: eskernel_beta(MAX_SPECIES) ! ES KERNEL: beta

  double precision,   save :: mass(MAX_SPECIES)
  double precision,   save :: nfrac(MAX_SPECIES)
  
  integer,            save :: particle_placement
  integer,            save :: particle_count(MAX_SPECIES)
  integer,            save :: p_move_tog(MAX_SPECIES)
  integer,            save :: p_force_tog(MAX_SPECIES)
  integer,            save :: p_int_tog(MAX_SPECIES*MAX_SPECIES)
  integer,            save :: p_int_tog_wall(MAX_SPECIES)
  double precision,   save :: particle_n0(MAX_SPECIES)
  double precision,   save :: particle_neff
  
  double precision,   save :: fixed_dt
  double precision,   save :: cfl

  double precision,   save :: rfd_delta

  integer,            save :: max_step
  integer,            save :: plot_int
  integer,            save :: plot_stag
  character(len=128), save :: plot_base_name
  integer,            save :: chk_int
  character(len=128), save :: chk_base_name
  integer,            save :: prob_type
  integer,            save :: restart
  integer,            save :: particle_restart
  integer,            save :: print_int
  integer,            save :: project_eos_int
  
  double precision,   save :: grav(AMREX_SPACEDIM)
  double precision,   save :: rhobar(MAX_SPECIES)
  double precision,   save :: rho0

  integer,            save :: nspecies
  double precision,   save :: molmass(MAX_SPECIES)
  double precision,   save :: diameter(MAX_SPECIES)

  integer,            save :: dof(MAX_SPECIES)
  double precision,   save :: hcv(MAX_SPECIES)
  double precision,   save :: hcp(MAX_SPECIES)

  double precision,   save :: variance_coef_mom
  double precision,   save :: variance_coef_mass
  double precision,   save :: k_B
  double precision,   save :: Runiv
  double precision,   save :: T_init(2)

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
  double precision,   save :: domega

  integer,            save :: bc_vel_lo(AMREX_SPACEDIM)
  integer,            save :: bc_vel_hi(AMREX_SPACEDIM)
  integer,            save :: bc_es_lo(AMREX_SPACEDIM)
  integer,            save :: bc_es_hi(AMREX_SPACEDIM)
  integer,            save :: bc_mass_lo(AMREX_SPACEDIM)
  integer,            save :: bc_mass_hi(AMREX_SPACEDIM)
  integer,            save :: bc_therm_lo(AMREX_SPACEDIM)
  integer,            save :: bc_therm_hi(AMREX_SPACEDIM)

  double precision,   save :: p_lo(AMREX_SPACEDIM)
  double precision,   save :: p_hi(AMREX_SPACEDIM)
  
  double precision,   save :: t_lo(AMREX_SPACEDIM)
  double precision,   save :: t_hi(AMREX_SPACEDIM)
  
  double precision,   save :: bc_Yk_x_lo(MAX_SPECIES)
  double precision,   save :: bc_Yk_x_hi(MAX_SPECIES)
  double precision,   save :: bc_Yk_y_lo(MAX_SPECIES)
  double precision,   save :: bc_Yk_y_hi(MAX_SPECIES)
  double precision,   save :: bc_Yk_z_lo(MAX_SPECIES)
  double precision,   save :: bc_Yk_z_hi(MAX_SPECIES)
  
  double precision,   save :: bc_Xk_x_lo(MAX_SPECIES)
  double precision,   save :: bc_Xk_x_hi(MAX_SPECIES)
  double precision,   save :: bc_Xk_y_lo(MAX_SPECIES)
  double precision,   save :: bc_Xk_y_hi(MAX_SPECIES)
  double precision,   save :: bc_Xk_z_lo(MAX_SPECIES)
  double precision,   save :: bc_Xk_z_hi(MAX_SPECIES)

  double precision,   save :: wallspeed_lo(AMREX_SPACEDIM-1,AMREX_SPACEDIM)
  double precision,   save :: wallspeed_hi(AMREX_SPACEDIM-1,AMREX_SPACEDIM)

  double precision,   save :: potential_lo(AMREX_SPACEDIM)
  double precision,   save :: potential_hi(AMREX_SPACEDIM)

  integer,            save :: struct_fact_int
  integer,            save :: radialdist_int
  integer,            save :: cartdist_int
  integer,            save :: n_steps_skip
  double precision,   save :: binsize
  double precision,   save :: searchdist

  integer,            save :: project_dir
  integer,            save :: slicepoint
  integer,            save :: max_grid_projection(AMREX_SPACEDIM-1)

  integer,            save :: histogram_unit
  double precision,   save :: density_weights(MAX_SPECIES)
  integer,            save :: shift_cc_to_boundary(AMREX_SPACEDIM,LOHI)

  double precision,   save :: permittivity
  integer,            save :: induced_charge_eo
  integer,            save :: zero_eps_on_wall_type
  double precision,   save :: zero_eps_wall_left_end
  double precision,   save :: zero_eps_wall_right_start
  integer,            save :: wall_mob
  double precision,   save :: rmin(MAX_SPECIES*MAX_SPECIES)
  double precision,   save :: rmax(MAX_SPECIES*MAX_SPECIES)
  double precision,   save :: eepsilon(MAX_SPECIES*MAX_SPECIES)
  double precision,   save :: sigma(MAX_SPECIES*MAX_SPECIES)

  double precision,   save :: rmin_wall(MAX_SPECIES)
  double precision,   save :: rmax_wall(MAX_SPECIES)
  double precision,   save :: eepsilon_wall(MAX_SPECIES)
  double precision,   save :: sigma_wall(MAX_SPECIES)
  
  integer,            save :: poisson_verbose
  integer,            save :: poisson_bottom_verbose
  integer,            save :: poisson_max_iter
  double precision,   save :: poisson_rel_tol

  double precision,   save :: particle_grid_refine
  double precision,   save :: es_grid_refine
  double precision,   save :: diff(MAX_SPECIES)
  integer,            save :: all_dry

  integer,            save :: fluid_tog
  integer,            save :: es_tog
  integer,            save :: drag_tog
  integer,            save :: move_tog
  integer,            save :: rfd_tog
  integer,            save :: dry_move_tog
  integer,            save :: sr_tog
  integer,            save :: graphene_tog
  integer,            save :: crange
  integer,            save :: thermostat_tog
  integer,            save :: zero_net_force
  integer,            save :: body_force

  integer,            save :: images
  double precision,   save :: eamp(3)
  double precision,   save :: efreq(3)
  double precision,   save :: ephase(3)
  double precision,   save :: F_body(3)

  integer,            save :: plot_ascii
  integer,            save :: plot_means
  integer,            save :: plot_vars
  integer,            save :: particle_motion

  integer,            save :: solve_chem
  double precision,   save :: diffcoeff
  double precision,   save :: scaling_factor
  double precision,   save :: source_strength
  integer,            save :: regrid_int
  integer,            save :: do_reflux

  double precision,   save :: turb_a
  double precision,   save :: turb_b
  integer,            save :: turbForcing

  ! Problem specification
  namelist /common/ prob_lo       ! physical lo coordinate
  namelist /common/ prob_hi       ! physical hi coordinate
  namelist /common/ n_cells       ! number of cells in domain
  namelist /common/ max_grid_size ! max number of cells in a box
  namelist /common/ max_grid_size_structfact ! max number of cells in a box for structure factor
  namelist /common/ max_particle_tile_size ! max number of cells in a box
  namelist /common/ cell_depth

  namelist /common/ ngc           ! number of ghost cells
  namelist /common/ nvars         ! number of conserved variables
  namelist /common/ nprimvars     ! number of primative variables
  namelist /common/ membrane_cell ! location of membrane
  namelist /common/ cross_cell    ! cell to compute spatial correlation
  namelist /common/ transmission

  namelist /common/ qval                ! charge on an ion
  namelist /common/ pkernel_fluid       ! peskin kernel for fluid
  namelist /common/ pkernel_es          ! peskin kernel for es
  namelist /common/ eskernel_fluid      ! ES kernel for fluid
  namelist /common/ eskernel_beta       ! ES kernel for fluid: beta

  namelist /common/ mass
  namelist /common/ nfrac

  namelist /common/ particle_placement
  namelist /common/ particle_count
  namelist /common/ p_move_tog
  namelist /common/ p_force_tog
  namelist /common/ p_int_tog
  namelist /common/ p_int_tog_wall
  namelist /common/ particle_n0
  namelist /common/ particle_neff

  ! Time-step control
  namelist /common/ fixed_dt
  namelist /common/ cfl

  !random finite difference size, fraction of cell size
  namelist /common/ rfd_delta

  ! Controls for number of steps between actions
  namelist /common/ max_step
  namelist /common/ plot_int
  namelist /common/ plot_stag
  namelist /common/ plot_base_name
  namelist /common/ chk_int
  namelist /common/ chk_base_name
  namelist /common/ prob_type
  namelist /common/ restart
  namelist /common/ particle_restart
  namelist /common/ print_int
  namelist /common/ project_eos_int

  ! Physical parameters
  namelist /common/ grav
  namelist /common/ rhobar
  namelist /common/ rho0

  ! Kinetic parameters
  namelist /common/ nspecies
  namelist /common/ molmass
  namelist /common/ diameter

  namelist /common/ dof
  namelist /common/ hcv
  namelist /common/ hcp

  ! stochastic forcing amplitudes (1 for physical values, 0 to run them off)
  namelist /common/ variance_coef_mom
  namelist /common/ variance_coef_mass
  namelist /common/ k_B
  namelist /common/ Runiv
  namelist /common/ T_init

  ! Algorithm control / selection
  namelist /common/ algorithm_type
  namelist /common/ advection_type
  namelist /common/ barodiffusion_type
  namelist /common/ use_bl_rng

  ! random number seed
  ! 0        = unpredictable seed based on clock
  ! positive = fixed seed
  namelist /common/ seed

  ! as assortment of other seeds in case one needs different engines
  ! implementation is problem-dependent
  ! 0        = unpredictable seed based on clock
  ! positive = fixed seed
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
  namelist /common/ domega

  ! Boundary conditions
  namelist /common/ bc_vel_lo
  namelist /common/ bc_vel_hi
  namelist /common/ bc_es_lo
  namelist /common/ bc_es_hi
  namelist /common/ bc_mass_lo
  namelist /common/ bc_mass_hi
  namelist /common/ bc_therm_lo
  namelist /common/ bc_therm_hi

  ! Pressure drop are periodic inflow/outflow walls (bc_[hi,lo]=-2).
  namelist /common/ p_lo
  namelist /common/ p_hi

  namelist /common/ t_lo
  namelist /common/ t_hi

  ! c_i boundary conditions
  namelist /common/ bc_Yk_x_lo
  namelist /common/ bc_Yk_x_hi
  namelist /common/ bc_Yk_y_lo
  namelist /common/ bc_Yk_y_hi
  namelist /common/ bc_Yk_z_lo
  namelist /common/ bc_Yk_z_hi

  namelist /common/ bc_Xk_x_lo
  namelist /common/ bc_Xk_x_hi
  namelist /common/ bc_Xk_y_lo
  namelist /common/ bc_Xk_y_hi
  namelist /common/ bc_Xk_z_lo
  namelist /common/ bc_Xk_z_hi

  ! Each no-slip wall may be moving with a specified tangential

  namelist /common/ wallspeed_lo
  namelist /common/ wallspeed_hi

  namelist /common/ potential_lo
  namelist /common/ potential_hi

  ! structure factor and radial/cartesian pair correlation function analysis
  namelist /common/ struct_fact_int
  namelist /common/ radialdist_int
  namelist /common/ cartdist_int
  namelist /common/ n_steps_skip
  namelist /common/ binsize
  namelist /common/ searchdist

  ! projection
  namelist /common/ project_dir
  namelist /common/ slicepoint
  namelist /common/ max_grid_projection

  ! These are mostly used for reaction-diffusion:
  namelist /common/ histogram_unit
  namelist /common/ density_weights
  namelist /common/ shift_cc_to_boundary

  namelist /common/ permittivity
  namelist /common/ induced_charge_eo
  namelist /common/ zero_eps_on_wall_type
  namelist /common/ zero_eps_wall_left_end
  namelist /common/ zero_eps_wall_right_start
  namelist /common/ wall_mob
  namelist /common/ rmin
  namelist /common/ rmax
  namelist /common/ eepsilon
  namelist /common/ sigma

  namelist /common/ rmin_wall
  namelist /common/ rmax_wall
  namelist /common/ eepsilon_wall
  namelist /common/ sigma_wall
  
  namelist /common/ poisson_verbose
  namelist /common/ poisson_bottom_verbose
  namelist /common/ poisson_max_iter
  namelist /common/ poisson_rel_tol

  namelist /common/ particle_grid_refine
  namelist /common/ es_grid_refine
  namelist /common/ diff
  namelist /common/ all_dry

  namelist /common/ fluid_tog
  namelist /common/ es_tog
  namelist /common/ drag_tog
  namelist /common/ move_tog
  namelist /common/ rfd_tog
  namelist /common/ dry_move_tog
  namelist /common/ sr_tog
  namelist /common/ graphene_tog
  namelist /common/ crange
  namelist /common/ thermostat_tog
  namelist /common/ zero_net_force
  namelist /common/ body_force


  namelist /common/ images
  namelist /common/ eamp
  namelist /common/ efreq
  namelist /common/ ephase
  namelist /common/ F_body

  namelist /common/ plot_ascii
  namelist /common/ plot_means
  namelist /common/ plot_vars
  namelist /common/ particle_motion

  ! chemistry
  namelist /common/ solve_chem
  namelist /common/ diffcoeff
  namelist /common/ source_strength
  namelist /common/ scaling_factor
  namelist /common/ regrid_int
  namelist /common/ do_reflux

  ! turblent forcing parameters
  namelist /common/ turb_a
  namelist /common/ turb_b
  namelist /common/ turbForcing

contains

  ! read in fortran namelist into common_params_module
  subroutine read_common_namelist(inputs_file,length) bind(C, name="read_common_namelist")

    integer               , value         :: length
    integer                               :: narg, farg
    character(kind=c_char), intent(in   ) :: inputs_file(length)
    character(len=128) :: fname
    
    ! default values
    prob_lo(:) = 0.d0
    prob_hi(:) = 1.d0
    n_cells(:) = 1
    max_grid_size(:) = 1
    max_grid_size_structfact(:) = 1
    max_particle_tile_size(:) = 0
    cell_depth = 1.d0

    ngc(:) = 1
    ! nvars (no default)
    ! nprimvars (no default)

    membrane_cell = -1
    cross_cell = 0
    ! transmission (no default)
    
    fixed_dt = 1.
    cfl = 0.5
    rfd_delta = 1e-5
    max_step = 1
    plot_int = 0
    plot_stag = 0
    plot_base_name = "plt"
    chk_int = 0
    chk_base_name = "chk"
    prob_type = 1
    restart = -1
    particle_restart = -1
    print_int = 0
    project_eos_int = -1
    grav(:) = 0.d0
    nspecies = 2
    molmass(:) = 1.d0
    diameter(:) = 1.d0
    rhobar(:) = 1.d0
    rho0 = 1.
    variance_coef_mom = 1.
    variance_coef_mass = 1.
    k_B = 1.38064852e-16
    Runiv = 8.314462175e7
    T_init(:) = 1.d0
    algorithm_type = 0
    advection_type = 0
    barodiffusion_type = 0
    use_bl_rng = 0
    seed = 0
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
    domega=0.d0
    bc_vel_lo(:) = 0
    bc_vel_hi(:) = 0
    bc_es_lo(:) = 0
    bc_es_hi(:) = 0
    bc_mass_lo(:) = 0
    bc_mass_hi(:) = 0
    bc_therm_lo(:) = 0
    bc_therm_hi(:) = 0

    drag_tog = 0
    particle_grid_refine = 1
    es_grid_refine = 1

    t_lo(:) = 0
    t_hi(:) = 0
    bc_Yk_x_lo(:) = 0.d0
    bc_Yk_x_hi(:) = 0.d0
    bc_Yk_y_lo(:) = 0.d0
    bc_Yk_y_hi(:) = 0.d0
    bc_Yk_z_lo(:) = 0.d0
    bc_Yk_z_hi(:) = 0.d0
    bc_Xk_x_lo(:) = 0.d0
    bc_Xk_x_hi(:) = 0.d0
    bc_Xk_y_lo(:) = 0.d0
    bc_Xk_y_hi(:) = 0.d0
    bc_Xk_z_lo(:) = 0.d0
    bc_Xk_z_hi(:) = 0.d0
    p_lo(:) = 0
    p_hi(:) = 0
    wallspeed_lo(:,:) = 0
    wallspeed_hi(:,:) = 0
    potential_lo(:) = 0
    potential_hi(:) = 0
    struct_fact_int = 0
    radialdist_int = 0
    cartdist_int = 0
    n_steps_skip = 0
    binsize = 0.
    searchDist = 0.
    project_dir = -1
    slicepoint = -1
    max_grid_projection(:) = 1
    histogram_unit = 0
    density_weights(:) = 0.d0
    shift_cc_to_boundary(:,:) = 0

    poisson_verbose = 1
    poisson_bottom_verbose = 0
    poisson_max_iter = 100
    poisson_rel_tol = 1.d-10

    p_move_tog(:) = 1
    p_force_tog(:) = 1
    p_int_tog(:) = 1

    pkernel_fluid(:) = 4
    pkernel_es(:) = 4
    eskernel_fluid(:) = -1
    eskernel_beta(:) = -1
    solve_chem = 0
    diffcoeff  = 0.001
    scaling_factor = 0.1
    source_strength = 0.1
    regrid_int = 25
    do_reflux  = 0

    turb_a = 1.d0
    turb_b = 1.d0
    turbForcing = 0

    plot_means = 0
    plot_vars = 0    
    particle_motion = 0

    graphene_tog = 0
    thermostat_tog = 0
    zero_net_force = 0
    body_force = 0

    all_dry = 0
    particle_neff = 1

    wall_mob = 1

    eamp(:) =  0
    efreq(:) = 0
    ephase(:) = 0
    F_body(:) = 0
    induced_charge_eo = 0
    zero_eps_on_wall_type = 0

    qval(:) = 0

    crange = maxval(pkernel_es) + 1

    ! read in common namelist
    open(unit=100, file=amrex_string_c_to_f(inputs_file), status='old', action='read')
    read(unit=100, nml=common)
    close(unit=100)

    if (nspecies > MAX_SPECIES) then
       call amrex_abort("nspecies > MAX_SPECIES; change in common_functions.H and common_namelist.F90")
    end if
    

  end subroutine read_common_namelist

  ! copy contents of common_params_module to C++ common namespace
  subroutine initialize_common_namespace(prob_lo_in, prob_hi_in, n_cells_in, &
                                         max_grid_size_in, max_grid_size_structfact_in, &
                                         max_particle_tile_size_in, cell_depth_in, ngc_in, &
                                         nvars_in, nprimvars_in, &
                                         membrane_cell_in, cross_cell_in, transmission_in, &
                                         qval_in, pkernel_fluid_in, pkernel_es_in,&
                                         eskernel_fluid_in, eskernel_beta_in,&
                                         fixed_dt_in, cfl_in, rfd_delta_in, max_step_in, plot_int_in, plot_stag_in, &
                                         plot_base_name_in, plot_base_name_len, chk_int_in, &
                                         chk_base_name_in, chk_base_name_len, prob_type_in, &
                                         restart_in, particle_restart_in, &
                                         print_int_in, project_eos_int_in, &
                                         grav_in, nspecies_in, molmass_in, diameter_in, &
                                         dof_in, hcv_in, hcp_in, rhobar_in, &
                                         rho0_in, variance_coef_mom_in, &
                                         variance_coef_mass_in, &
                                         k_B_in, Runiv_in, T_init_in, algorithm_type_in, &
                                         advection_type_in, &
                                         barodiffusion_type_in, use_bl_rng_in, seed_in, &
                                         seed_momentum_in, seed_diffusion_in, &
                                         seed_reaction_in, &
                                         seed_init_mass_in, seed_init_momentum_in, &
                                         visc_coef_in, visc_type_in, &
                                         filtering_width_in, stoch_stress_form_in, &
                                         u_init_in, perturb_width_in, smoothing_width_in, &
                                         initial_variance_mom_in, initial_variance_mass_in, &
                                         domega_in, bc_vel_lo_in, bc_vel_hi_in, &
                                         bc_es_lo_in, bc_es_hi_in,  &
                                         bc_mass_lo_in, bc_mass_hi_in,  &
                                         bc_therm_lo_in, bc_therm_hi_in,  &
                                         p_lo_in, p_hi_in, &
                                         t_lo_in, t_hi_in, &
                                         bc_Yk_x_lo_in, bc_Yk_x_hi_in, &
                                         bc_Yk_y_lo_in, bc_Yk_y_hi_in, &
                                         bc_Yk_z_lo_in, bc_Yk_z_hi_in, &
                                         bc_Xk_x_lo_in, bc_Xk_x_hi_in, &
                                         bc_Xk_y_lo_in, bc_Xk_y_hi_in, &
                                         bc_Xk_z_lo_in, bc_Xk_z_hi_in, &
                                         wallspeed_lo_in, wallspeed_hi_in, &
                                         potential_lo_in, potential_hi_in, &
                                         struct_fact_int_in, radialdist_int_in, &
                                         cartdist_int_in, n_steps_skip_in, &
                                         binsize_in, searchdist_in, &
                                         project_dir_in, slicepoint_in, max_grid_projection_in, &
                                         histogram_unit_in, density_weights_in, &
                                         shift_cc_to_boundary_in, &
                                         particle_placement_in, particle_count_in, p_move_tog_in, &
                                         p_force_tog_in, p_int_tog_in, p_int_tog_wall_in, particle_neff_in,&
                                         particle_n0_in, mass_in, nfrac_in, permittivity_in, &
                                         induced_charge_eo_in, zero_eps_on_wall_type_in, zero_eps_wall_left_end_in, &
                                         zero_eps_wall_right_start_in, &
                                         wall_mob_in, rmin_in, rmax_in, eepsilon_in, sigma_in, rmin_wall_in, rmax_wall_in, eepsilon_wall_in, sigma_wall_in, poisson_verbose_in, &
                                         poisson_bottom_verbose_in, poisson_max_iter_in, poisson_rel_tol_in, &
                                         particle_grid_refine_in, es_grid_refine_in, diff_in, all_dry_in, &
                                         fluid_tog_in, es_tog_in, drag_tog_in, move_tog_in, rfd_tog_in, &
                                         dry_move_tog_in, sr_tog_in, graphene_tog_in, crange_in, &
                                         thermostat_tog_in, zero_net_force_in, body_force_in, images_in, eamp_in, efreq_in, ephase_in, F_body_in, &
                                         plot_ascii_in, plot_means_in, plot_vars_in, &
                                         solve_chem_in, diffcoeff_in, scaling_factor_in, &
                                         source_strength_in, regrid_int_in, do_reflux_in, particle_motion_in, &
                                         turb_a_in, turb_b_in, turbForcing_in) &
                                         bind(C, name="initialize_common_namespace")

    double precision,       intent(inout) :: prob_lo_in(AMREX_SPACEDIM)
    double precision,       intent(inout) :: prob_hi_in(AMREX_SPACEDIM)
    integer,                intent(inout) :: n_cells_in(AMREX_SPACEDIM)
    integer,                intent(inout) :: max_grid_size_in(AMREX_SPACEDIM)
    integer,                intent(inout) :: max_grid_size_structfact_in(AMREX_SPACEDIM)
    integer,                intent(inout) :: max_particle_tile_size_in(AMREX_SPACEDIM)
    double precision,       intent(inout) :: cell_depth_in

    double precision,       intent(inout) :: mass_in(MAX_SPECIES)
    double precision,       intent(inout) :: nfrac_in(MAX_SPECIES)
    double precision,       intent(inout) :: particle_n0_in(MAX_SPECIES)
    double precision,       intent(inout) :: particle_neff_in
    integer,                intent(inout) :: particle_count_in(MAX_SPECIES)
    integer,                intent(inout) :: p_move_tog_in(MAX_SPECIES)
    integer,                intent(inout) :: p_force_tog_in(MAX_SPECIES)
    integer,                intent(inout) :: p_int_tog_in(MAX_SPECIES*MAX_SPECIES)
    integer,                intent(inout) :: p_int_tog_wall_in(MAX_SPECIES)
    integer,                intent(inout) :: particle_placement_in
    
    double precision,       intent(inout) :: fixed_dt_in
    double precision,       intent(inout) :: cfl_in
    double precision,       intent(inout) :: rfd_delta_in

    integer,                intent(inout) :: ngc_in(AMREX_SPACEDIM)
    integer,                intent(inout) :: nvars_in
    integer,                intent(inout) :: nprimvars_in

    integer,                intent(inout) :: membrane_cell_in
    integer,                intent(inout) :: cross_cell_in
    double precision,       intent(inout) :: transmission_in

    double precision,       intent(inout) :: qval_in(MAX_SPECIES)
    integer,                intent(inout) :: pkernel_fluid_in(MAX_SPECIES)
    integer,                intent(inout) :: pkernel_es_in(MAX_SPECIES)
    integer,                intent(inout) :: eskernel_fluid_in(MAX_SPECIES)
    double precision,       intent(inout) :: eskernel_beta_in(MAX_SPECIES)

    integer,                intent(inout) :: max_step_in
    integer,                intent(inout) :: plot_int_in
    integer,                intent(inout) :: plot_stag_in
    integer               , value         :: plot_base_name_len
    character(kind=c_char), intent(inout) :: plot_base_name_in(plot_base_name_len)
    integer,                intent(inout) :: chk_int_in
    integer               , value         :: chk_base_name_len
    character(kind=c_char), intent(inout) :: chk_base_name_in(chk_base_name_len)
    integer,                intent(inout) :: prob_type_in
    integer,                intent(inout) :: restart_in
    integer,                intent(inout) :: particle_restart_in
    integer,                intent(inout) :: print_int_in
    integer,                intent(inout) :: project_eos_int_in
    double precision,       intent(inout) :: grav_in(AMREX_SPACEDIM)
    integer,                intent(inout) :: nspecies_in
    double precision,       intent(inout) :: molmass_in(MAX_SPECIES)
    double precision,       intent(inout) :: diameter_in(MAX_SPECIES)
    integer,                intent(inout) :: dof_in(MAX_SPECIES)
    double precision,       intent(inout) :: hcv_in(MAX_SPECIES)
    double precision,       intent(inout) :: hcp_in(MAX_SPECIES)

    double precision,       intent(inout) :: rhobar_in(MAX_SPECIES)
    double precision,       intent(inout) :: rho0_in
    double precision,       intent(inout) :: variance_coef_mom_in
    double precision,       intent(inout) :: variance_coef_mass_in
    double precision,       intent(inout) :: k_B_in
    double precision,       intent(inout) :: Runiv_in
    double precision,       intent(inout) :: T_init_in(2)
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
    double precision,       intent(inout) :: domega_in
    integer,                intent(inout) :: bc_vel_lo_in(AMREX_SPACEDIM)
    integer,                intent(inout) :: bc_vel_hi_in(AMREX_SPACEDIM)
    integer,                intent(inout) :: bc_es_lo_in(AMREX_SPACEDIM)
    integer,                intent(inout) :: bc_es_hi_in(AMREX_SPACEDIM)
    integer,                intent(inout) :: bc_mass_lo_in(AMREX_SPACEDIM)
    integer,                intent(inout) :: bc_mass_hi_in(AMREX_SPACEDIM)
    integer,                intent(inout) :: bc_therm_lo_in(AMREX_SPACEDIM)
    integer,                intent(inout) :: bc_therm_hi_in(AMREX_SPACEDIM)

    double precision,       intent(inout) :: p_lo_in(AMREX_SPACEDIM)
    double precision,       intent(inout) :: p_hi_in(AMREX_SPACEDIM)
    double precision,       intent(inout) :: t_lo_in(AMREX_SPACEDIM)
    double precision,       intent(inout) :: t_hi_in(AMREX_SPACEDIM)
    double precision,       intent(inout) :: bc_Yk_x_lo_in(MAX_SPECIES)
    double precision,       intent(inout) :: bc_Yk_x_hi_in(MAX_SPECIES)
    double precision,       intent(inout) :: bc_Yk_y_lo_in(MAX_SPECIES)
    double precision,       intent(inout) :: bc_Yk_y_hi_in(MAX_SPECIES)
    double precision,       intent(inout) :: bc_Yk_z_lo_in(MAX_SPECIES)
    double precision,       intent(inout) :: bc_Yk_z_hi_in(MAX_SPECIES)
    double precision,       intent(inout) :: bc_Xk_x_lo_in(MAX_SPECIES)
    double precision,       intent(inout) :: bc_Xk_x_hi_in(MAX_SPECIES)
    double precision,       intent(inout) :: bc_Xk_y_lo_in(MAX_SPECIES)
    double precision,       intent(inout) :: bc_Xk_y_hi_in(MAX_SPECIES)
    double precision,       intent(inout) :: bc_Xk_z_lo_in(MAX_SPECIES)
    double precision,       intent(inout) :: bc_Xk_z_hi_in(MAX_SPECIES)
    double precision,       intent(inout) :: wallspeed_lo_in(AMREX_SPACEDIM-1,AMREX_SPACEDIM)
    double precision,       intent(inout) :: wallspeed_hi_in(AMREX_SPACEDIM-1,AMREX_SPACEDIM)

    double precision,       intent(inout) :: potential_lo_in(AMREX_SPACEDIM)
    double precision,       intent(inout) :: potential_hi_in(AMREX_SPACEDIM)

    integer,                intent(inout) :: struct_fact_int_in
    integer,                intent(inout) :: radialdist_int_in
    integer,                intent(inout) :: cartdist_int_in
    integer,                intent(inout) :: n_steps_skip_in
    double precision,       intent(inout) :: binsize_in
    double precision,       intent(inout) :: searchdist_in
    integer,                intent(inout) :: project_dir_in
    integer,                intent(inout) :: slicepoint_in
    integer,                intent(inout) :: max_grid_projection_in(AMREX_SPACEDIM-1)
    integer,                intent(inout) :: histogram_unit_in
    double precision,       intent(inout) :: density_weights_in(MAX_SPECIES)
    integer,                intent(inout) :: shift_cc_to_boundary_in(AMREX_SPACEDIM,LOHI)

    double precision,       intent(inout) :: eepsilon_in(MAX_SPECIES*MAX_SPECIES)
    double precision,       intent(inout) :: sigma_in(MAX_SPECIES*MAX_SPECIES)
    double precision,       intent(inout) :: eepsilon_wall_in(MAX_SPECIES)
    double precision,       intent(inout) :: sigma_wall_in(MAX_SPECIES)

    double precision,       intent(inout) :: permittivity_in
    integer,                intent(inout) :: induced_charge_eo_in
    integer,                intent(inout) :: zero_eps_on_wall_type_in
    double precision,       intent(inout) :: zero_eps_wall_left_end_in
    double precision,       intent(inout) :: zero_eps_wall_right_start_in
    integer,                intent(inout) :: wall_mob_in
    double precision,       intent(inout) :: rmin_in(MAX_SPECIES*MAX_SPECIES)
    double precision,       intent(inout) :: rmax_in(MAX_SPECIES*MAX_SPECIES)
    double precision,       intent(inout) :: rmin_wall_in(MAX_SPECIES)
    double precision,       intent(inout) :: rmax_wall_in(MAX_SPECIES)
    double precision,       intent(inout) :: poisson_rel_tol_in

    integer,                intent(inout) :: poisson_max_iter_in
    integer,                intent(inout) :: poisson_verbose_in
    integer,                intent(inout) :: poisson_bottom_verbose_in

    double precision,       intent(inout) :: particle_grid_refine_in
    double precision,       intent(inout) :: es_grid_refine_in
    double precision,       intent(inout) :: diff_in(MAX_SPECIES)
    integer,                intent(inout) :: all_dry_in
    
    integer,                intent(inout) :: fluid_tog_in
    integer,                intent(inout) :: es_tog_in
    integer,                intent(inout) :: drag_tog_in
    integer,                intent(inout) :: move_tog_in
    integer,                intent(inout) :: rfd_tog_in
    integer,                intent(inout) :: dry_move_tog_in
    integer,                intent(inout) :: sr_tog_in
    integer,                intent(inout) :: crange_in
    integer,                intent(inout) :: graphene_tog_in
    integer,                intent(inout) :: thermostat_tog_in
    integer,                intent(inout) :: zero_net_force_in
    integer,                intent(inout) :: body_force_in

    integer,                intent(inout) :: images_in
    double precision,       intent(inout) :: eamp_in(3)
    double precision,       intent(inout) :: efreq_in(3)
    double precision,       intent(inout) :: ephase_in(3)
    double precision,       intent(inout) :: F_body_in(3)

    integer,                intent(inout) :: plot_ascii_in
    integer,                intent(inout) :: plot_means_in
    integer,                intent(inout) :: plot_vars_in
    
    integer,                intent(inout) :: solve_chem_in
    double precision,       intent(inout) :: diffcoeff_in
    double precision,       intent(inout) :: scaling_factor_in
    double precision,       intent(inout) :: source_strength_in
    integer,                intent(inout) :: regrid_int_in
    integer,                intent(inout) :: do_reflux_in
    integer,                intent(inout) :: particle_motion_in

    double precision,       intent(inout) :: turb_a_in
    double precision,       intent(inout) :: turb_b_in
    integer,                intent(inout) :: turbForcing_in

    prob_lo_in = prob_lo
    prob_hi_in = prob_hi
    n_cells_in = n_cells
    max_grid_size_in = max_grid_size
    max_grid_size_structfact_in = max_grid_size_structfact
    max_particle_tile_size_in = max_particle_tile_size
    cell_depth_in = cell_depth
    ngc_in = ngc
    nvars_in = nvars
    nprimvars_in = nprimvars
    membrane_cell_in = membrane_cell
    cross_cell_in = cross_cell
    transmission_in = transmission

    qval_in = qval
    pkernel_fluid_in = pkernel_fluid
    pkernel_es_in = pkernel_es
    eskernel_fluid_in = eskernel_fluid
    eskernel_beta_in = eskernel_beta

    fixed_dt_in = fixed_dt
    cfl_in = cfl
    rfd_delta_in = rfd_delta
    max_step_in = max_step
    plot_int_in = plot_int
    plot_stag_in = plot_stag
    plot_base_name_in = amrex_string_f_to_c(plot_base_name)
    chk_int_in = chk_int
    chk_base_name_in = amrex_string_f_to_c(chk_base_name)
    prob_type_in = prob_type
    restart_in = restart
    particle_restart_in = particle_restart
    print_int_in = print_int
    project_eos_int_in = project_eos_int
    grav_in = grav
    nspecies_in = nspecies
    molmass_in = molmass
    rhobar_in = rhobar
    diameter_in = diameter
    dof_in = dof
    hcv_in = hcv
    hcp_in = hcp
    rho0_in= rho0
    variance_coef_mom_in = variance_coef_mom
    variance_coef_mass_in = variance_coef_mass
    k_B_in = k_B
    Runiv_in = Runiv
    T_init_in = T_init
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
    domega_in=domega
    bc_vel_lo_in = bc_vel_lo
    bc_vel_hi_in = bc_vel_hi
    bc_es_lo_in = bc_es_lo
    bc_es_hi_in = bc_es_hi
    bc_mass_lo_in = bc_mass_lo
    bc_mass_hi_in = bc_mass_hi
    bc_therm_lo_in = bc_therm_lo
    bc_therm_hi_in = bc_therm_hi

    p_lo_in = p_lo
    p_hi_in = p_hi
    t_lo_in = t_lo
    t_hi_in = t_hi
    bc_Yk_x_lo_in = bc_Yk_x_lo
    bc_Yk_x_hi_in = bc_Yk_x_hi
    bc_Yk_y_lo_in = bc_Yk_y_lo
    bc_Yk_y_hi_in = bc_Yk_y_hi
    bc_Yk_z_lo_in = bc_Yk_z_lo
    bc_Yk_z_hi_in = bc_Yk_z_hi
    bc_Xk_x_lo_in = bc_Xk_x_lo
    bc_Xk_x_hi_in = bc_Xk_x_hi
    bc_Xk_y_lo_in = bc_Xk_y_lo
    bc_Xk_y_hi_in = bc_Xk_y_hi
    bc_Xk_z_lo_in = bc_Xk_z_lo
    bc_Xk_z_hi_in = bc_Xk_z_hi
    wallspeed_lo_in = wallspeed_lo
    wallspeed_hi_in = wallspeed_hi

    potential_lo_in = potential_lo
    potential_hi_in = potential_hi

    struct_fact_int_in = struct_fact_int
    radialdist_int_in = radialdist_int
    cartdist_int_in = cartdist_int
    n_steps_skip_in = n_steps_skip
    binsize_in = binsize
    searchdist_in = searchdist
    project_dir_in = project_dir
    slicepoint_in = slicepoint
    max_grid_projection_in = max_grid_projection
    histogram_unit_in = histogram_unit
    density_weights_in = density_weights
    shift_cc_to_boundary_in = shift_cc_to_boundary

    mass_in = mass
    nfrac_in = nfrac
    particle_n0_in = particle_n0
    particle_neff_in = particle_neff
    particle_count_in = particle_count
    p_move_tog_in = p_move_tog
    p_force_tog_in = p_force_tog
    p_int_tog_in = p_int_tog
    p_int_tog_wall_in = p_int_tog_wall
    particle_placement_in = particle_placement

    poisson_verbose_in = poisson_verbose
    poisson_bottom_verbose_in = poisson_bottom_verbose
    poisson_max_iter_in = poisson_max_iter
    poisson_rel_tol_in = poisson_rel_tol
    permittivity_in = permittivity
    induced_charge_eo_in = induced_charge_eo
    zero_eps_on_wall_type_in = zero_eps_on_wall_type
    zero_eps_wall_left_end_in = zero_eps_wall_left_end
    zero_eps_wall_right_start_in = zero_eps_wall_right_start
    wall_mob_in = wall_mob
    rmin_in = rmin
    rmax_in = rmax
    eepsilon_in = eepsilon
    sigma_in = sigma

    rmin_wall_in = rmin_wall
    rmax_wall_in = rmax_wall
    eepsilon_wall_in = eepsilon_wall
    sigma_wall_in = sigma_wall

    particle_grid_refine_in = particle_grid_refine
    es_grid_refine_in = es_grid_refine
    diff_in = diff
    all_dry_in = all_dry
    
    fluid_tog_in = fluid_tog
    es_tog_in = es_tog
    drag_tog_in = drag_tog
    move_tog_in = move_tog
    rfd_tog_in = rfd_tog
    dry_move_tog_in = dry_move_tog
    sr_tog_in = sr_tog
    crange_in = crange
    graphene_tog_in = graphene_tog
    thermostat_tog_in = thermostat_tog
    zero_net_force_in = zero_net_force
    body_force_in = body_force

    images_in = images
    eamp_in = eamp
    efreq_in = efreq
    ephase_in = ephase
    F_body_in = F_body

    plot_ascii_in = plot_ascii
    plot_means_in = plot_means
    plot_vars_in = plot_vars
    
    solve_chem_in = solve_chem
    diffcoeff_in  = diffcoeff
    scaling_factor_in  = scaling_factor
    source_strength_in  = source_strength

    regrid_int_in = regrid_int
    do_reflux_in  = do_reflux
    particle_motion_in = particle_motion

    turb_a_in = turb_a
    turb_b_in = turb_b
    turbForcing_in = turbForcing

  end subroutine initialize_common_namespace

  subroutine set_max_step(max_step_in) bind(C, name="set_max_step")

    integer, intent(in   ) :: max_step_in

    max_step = max_step_in
    
  end subroutine set_max_step

  subroutine set_domega(domega_in) bind(C, name="set_domega")

    double precision, intent(in   ) :: domega_in

    domega = domega_in
    
  end subroutine set_domega  

end module common_namelist_module
