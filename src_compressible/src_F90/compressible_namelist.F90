module compressible_namelist_module

  use iso_c_binding, only: c_char
  use amrex_string_module, only: amrex_string_c_to_f, amrex_string_f_to_c

  implicit none

  integer, parameter :: MAX_SPECIES = 10
  integer, parameter :: LOHI = 2

  integer,            save :: mass_bc_lo(AMREX_SPACEDIM)
  integer,            save :: mass_bc_hi(AMREX_SPACEDIM)
  integer,            save :: therm_bc_lo(AMREX_SPACEDIM)
  integer,            save :: therm_bc_hi(AMREX_SPACEDIM)
  integer,            save :: vel_bc_lo(AMREX_SPACEDIM)
  integer,            save :: vel_bc_hi(AMREX_SPACEDIM)
  double precision,   save :: Yk_bc(AMREX_SPACEDIM,LOHI,MAX_SPECIES)
  double precision,   save :: Xk_bc(AMREX_SPACEDIM,LOHI,MAX_SPECIES)

  ! Boundary conditions
  namelist /compressible/ mass_bc_lo  ! mass boundary conditions
  namelist /compressible/ mass_bc_hi
  namelist /compressible/ therm_bc_lo ! temperature boundary conditions
  namelist /compressible/ therm_bc_hi
  namelist /compressible/ vel_bc_lo   ! velocity field boundary conditions
  namelist /compressible/ vel_bc_hi
  namelist /compressible/ Yk_bc       ! mass fraction wall boundary value
  namelist /compressible/ Xk_bc       ! mole fraction wall boundary value

contains

  ! read in fortran namelist into compressible_params_module
  subroutine read_compressible_namelist(inputs_file,length) bind(C, name="read_compressible_namelist")

    integer               , value         :: length
    character(kind=c_char), intent(in   ) :: inputs_file(length)

    ! default values
    mass_bc_lo(:) = 0
    mass_bc_hi(:) = 0
    therm_bc_lo(:) = 0
    therm_bc_hi(:) = 0
    vel_bc_lo(:) = 0
    vel_bc_hi(:) = 0
    Yk_bc(:,:,:) = 0.d0
    Xk_bc(:,:,:) = 0.d0

    ! read in compressible namelist
    open(unit=100, file=amrex_string_c_to_f(inputs_file), status='old', action='read')
    read(unit=100, nml=compressible)
    close(unit=100)

  end subroutine read_compressible_namelist

  ! copy contents of compressible_params_module to C++ compressible namespace
  subroutine initialize_compressible_namespace( mass_bc_lo_in, mass_bc_hi_in, &
                                                therm_bc_lo_in, therm_bc_hi_in, &
                                                vel_bc_lo_in, vel_bc_hi_in, &
                                                Yk_bc_in, Xk_bc_in) &
                                                bind(C, name="initialize_compressible_namespace")

    double precision,       intent(inout) :: mass_bc_lo_in(AMREX_SPACEDIM), mass_bc_hi_in(AMREX_SPACEDIM), therm_bc_lo_in(AMREX_SPACEDIM), therm_bc_hi_in(AMREX_SPACEDIM), vel_bc_lo_in(AMREX_SPACEDIM), vel_bc_hi_in(AMREX_SPACEDIM),Yk_bc_in(AMREX_SPACEDIM,LOHI,MAX_SPECIES), Xk_bc_in(AMREX_SPACEDIM,LOHI,MAX_SPECIES)

    mass_bc_lo_in = mass_bc_lo
    mass_bc_hi_in = mass_bc_hi
    therm_bc_lo_in = therm_bc_lo
    therm_bc_hi_in = therm_bc_hi
    vel_bc_lo_in = vel_bc_lo
    vel_bc_hi_in = vel_bc_hi
    Yk_bc_in = Yk_bc
    Xk_bc_in = Xk_bc

  end subroutine initialize_compressible_namespace

end module compressible_namelist_module
