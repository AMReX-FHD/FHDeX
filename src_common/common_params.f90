module common_params_module

  use amrex_error_module

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

  integer, allocatable, save :: test_array(:,:,:)
  namelist /probin_common/ test_array

  ! End the declarations of the ParmParse parameters

contains

  subroutine read_common_params() bind(C, name="read_common_params")

    use amrex_parmparse_module, only: amrex_parmparse_build, &
                                      amrex_parmparse_destroy, &
                                      amrex_parmparse

    type (amrex_parmparse) :: pp

    ! for reading in namelist from probin
    integer :: un
    logical :: lexist

    integer :: i,j,k

    ! allocate arrays
    allocate(prob_lo(MAX_SPACEDIM))
    allocate(prob_hi(MAX_SPACEDIM))
    allocate(n_cells(MAX_SPACEDIM))
    allocate(max_grid_size(MAX_SPACEDIM))
    allocate(character(len=1)::plot_base_name)
    allocate(character(len=1)::chk_base_name)
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

    allocate(test_array(3,3,3))

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
    barodiffusion_type = 0
    use_bl_rng = .false.
    seed = 1
    seed_momentum = 1
    seed_diffusion = 1
    seed_reaction = 1
    seed_init_mass = 1
    seed_init_momentum = 1
    visc_coef = 1.
    visc_type = 1
    advection_type = 0
    filtering_width = 0
    stoch_stress_form = 1
    u_init(:) = 0.d0
    perturb_width = 0.
    smoothing_width = 1.
    initial_variance_mom = 0.
    initial_variance_mass = 0.
    bc_lo(:) = 0
    bc_hi(:) = 0
    wallspeed_lo_x(:) = 0
    wallspeed_hi_x(:) = 0
    wallspeed_lo_y(:) = 0
    wallspeed_hi_y(:) = 0
    wallspeed_lo_z(:) = 0
    wallspeed_hi_z(:) = 0
    histogram_unit = 0
    density_weights(:) = 0.d0
    shift_cc_to_boundary_lo(:) = 0
    shift_cc_to_boundary_hi(:) = 0

    test_array(:,:,:) = 0

    ! read in from inputs file    
    call amrex_parmparse_build(pp)

    call pp%queryarr("prob_lo",prob_lo);
    call pp%queryarr("prob_hi",prob_hi);
    call pp%queryarr("n_cells",n_cells);
    call pp%queryarr("max_grid_size",max_grid_size);
    call pp%query("fixed_dt",fixed_dt);
    call pp%query("cfl",cfl);
    call pp%query("max_step",max_step);
    call pp%query("plot_int",plot_int);
    call pp%query("plot_base_name",plot_base_name);
    call pp%query("chk_int",chk_int);
    call pp%query("chk_base_name",chk_base_name);
    call pp%query("prob_type",prob_type);
    call pp%query("restart",restart);
    call pp%query("print_int",print_int);
    call pp%query("project_eos_int",project_eos_int);
    call pp%queryarr("grav",grav);
    call pp%query("nspecies",nspecies);
    call pp%queryarr("molmass",molmass);
    call pp%queryarr("rhobar",rhobar);
    call pp%query("rho0",rho0);
    call pp%query("variance_coef_mom",variance_coef_mom);
    call pp%query("variance_coef_mass",variance_coef_mass);
    call pp%query("k_B",k_B);
    call pp%query("Runiv",Runiv);
    call pp%query("algorithm_type",algorithm_type);
    call pp%query("barodiffusion_type",barodiffusion_type);
    call pp%query("use_bl_rng",use_bl_rng);
    call pp%query("seed",seed);
    call pp%query("seed_momentum",seed_momentum);
    call pp%query("seed_diffusion",seed_diffusion);
    call pp%query("seed_reaction",seed_reaction);
    call pp%query("seed_init_mass",seed_init_mass);
    call pp%query("seed_init_momentum",seed_init_momentum);
    call pp%query("visc_coef",visc_coef);
    call pp%query("visc_type",visc_type);
    call pp%query("advection_type",advection_type);
    call pp%query("filtering_width",filtering_width);
    call pp%query("stoch_stress_form",stoch_stress_form);
    call pp%queryarr("u_init",u_init);
    call pp%query("perturb_width",perturb_width);
    call pp%query("smoothing_width",smoothing_width);
    call pp%query("initial_variance_mom",initial_variance_mom);
    call pp%query("initial_variance_mass",initial_variance_mass);
    call pp%queryarr("bc_lo",bc_lo);
    call pp%queryarr("bc_hi",bc_hi);
    call pp%queryarr("wallspeed_lo_x",wallspeed_lo_x);
    call pp%queryarr("wallspeed_hi_x",wallspeed_hi_x);
    call pp%queryarr("wallspeed_lo_y",wallspeed_lo_y);
    call pp%queryarr("wallspeed_hi_y",wallspeed_hi_y);
    call pp%queryarr("wallspeed_lo_z",wallspeed_lo_z);
    call pp%queryarr("wallspeed_hi_z",wallspeed_hi_z);
    call pp%query("histogram_unit",histogram_unit);
    call pp%queryarr("density_weights",density_weights);
    call pp%queryarr("shift_cc_to_boundary_lo",shift_cc_to_boundary_lo);
    call pp%queryarr("shift_cc_to_boundary_hi",shift_cc_to_boundary_hi);

    call amrex_parmparse_destroy(pp)

    ! now the multidimensional arrays in the namelist
    ! these exist in common_params_module (fortran) only
    ! if we end up needing these in C++ we'll have to make accessor functions
    inquire(file = "probin_2d", exist = lexist )
    if ( lexist ) then
       un = 100 ! each namelist needs a different unit number
       open(unit = un, file = "probin_2d", status = 'old', action = 'read')
       read(unit = un, nml = probin_common)
       close(unit = un)
    else
       call amrex_error("invalid probin file")
    end if

  end subroutine read_common_params

  subroutine finalize_common_params() bind(C, name="finalize_common_params")

    ! deallocate parameters to quiet valgrind
    deallocate(prob_lo)
    deallocate(prob_hi)
    deallocate(n_cells)
    deallocate(max_grid_size)
    deallocate(plot_base_name)
    deallocate(chk_base_name)
    deallocate(grav)
    deallocate(molmass)
    deallocate(rhobar)
    deallocate(u_init)
    deallocate(bc_lo)
    deallocate(bc_hi)
    deallocate(wallspeed_lo_x)
    deallocate(wallspeed_hi_x)
    deallocate(wallspeed_lo_y)
    deallocate(wallspeed_hi_y)
    deallocate(wallspeed_lo_z)
    deallocate(wallspeed_hi_z)
    deallocate(density_weights)
    deallocate(shift_cc_to_boundary_lo)
    deallocate(shift_cc_to_boundary_hi)

  end subroutine finalize_common_params

end module common_params_module
