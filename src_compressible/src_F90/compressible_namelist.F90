module compressible_namelist_module

  use iso_c_binding, only: c_char
  use amrex_string_module, only: amrex_string_c_to_f, amrex_string_f_to_c

  implicit none

  integer, parameter :: MAX_SPECIES = 10
  integer, parameter :: LOHI = 2

  double precision,   save :: Yk_lo(AMREX_SPACEDIM,MAX_SPECIES)
  double precision,   save :: Yk_hi(AMREX_SPACEDIM,MAX_SPECIES)
  double precision,   save :: Xk_lo(AMREX_SPACEDIM,MAX_SPECIES)
  double precision,   save :: Xk_hi(AMREX_SPACEDIM,MAX_SPECIES)
  
  ! Multispecies parameters
  namelist /compressible/ Yk_lo       ! lo mass fraction wall boundary value
  namelist /compressible/ Yk_hi       ! hi mass fraction wall boundary value
  namelist /compressible/ Xk_lo       ! lo mole fraction wall boundary value
  namelist /compressible/ Xk_hi       ! hi mole fraction wall boundary value

contains

  ! read in fortran namelist into compressible_params_module
  subroutine read_compressible_namelist(inputs_file,length) bind(C, name="read_compressible_namelist")

    integer               , value         :: length
    character(kind=c_char), intent(in   ) :: inputs_file(length)

    ! default values
    Yk_lo(:,:) = -1.d0
    Yk_hi(:,:) = -1.d0
    Xk_lo(:,:) = -1.d0
    Xk_hi(:,:) = -1.d0

    ! read in compressible namelist
    open(unit=100, file=amrex_string_c_to_f(inputs_file), status='old', action='read')
    ! print*, inputs_file
    ! stop
    read(unit=100, nml=compressible)
    close(unit=100)

  end subroutine read_compressible_namelist

  ! copy contents of compressible_params_module to C++ compressible namespace
  subroutine initialize_compressible_namespace( Yk_lo_in, Yk_hi_in, &
                                                Xk_lo_in, Xk_hi_in) &
                                               bind(C, name="initialize_compressible_namespace")

    double precision,       intent(inout) :: Yk_lo_in(AMREX_SPACEDIM,MAX_SPECIES), Yk_hi_in(AMREX_SPACEDIM,MAX_SPECIES) 
    double precision,       intent(inout) :: Xk_lo_in(AMREX_SPACEDIM,MAX_SPECIES), Xk_hi_in(AMREX_SPACEDIM,MAX_SPECIES)

    Yk_lo_in = Yk_lo
    Yk_hi_in = Yk_hi
    Xk_lo_in = Xk_lo
    Xk_hi_in = Xk_hi

  end subroutine initialize_compressible_namespace

end module compressible_namelist_module
