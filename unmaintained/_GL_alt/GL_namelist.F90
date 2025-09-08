module GL_namelist_module

  use iso_c_binding, only: c_char
  use amrex_string_module, only: amrex_string_c_to_f, amrex_string_f_to_c

  implicit none


  double precision,   save :: acoef,bcoef,ccoef,dcoef
  double precision,   save :: umbrella
  double precision,   save :: diff_coef
  double precision,   save :: noise_coef
  double precision,   save :: phi0
  double precision,   save :: phi_inc
  double precision,   save :: rad
  integer,            save :: n_inc_phi
  integer,            save :: Plot_Skip
  integer,            save :: Number_of_Samples
  integer,            save :: Equil
  double precision,   save :: alpha
  double precision,   save :: r1
  double precision,   save :: r2
  double precision,   save :: umbrella_max
  double precision,   save :: umbrella_min
  integer,            save :: adaptive
  integer,            save :: Reverse








  ! Problem specification
  namelist /GL_params/ acoef
  namelist /GL_params/ bcoef
  namelist /GL_params/ ccoef
  namelist /GL_params/ dcoef
  namelist /GL_params/ umbrella
  namelist /GL_params/ diff_coef
  namelist /GL_params/ noise_coef
  namelist /GL_params/ phi0
  namelist /GL_params/ phi_inc
  namelist /GL_params/ n_inc_phi
  namelist /GL_params/ rad
  namelist /GL_params/ Plot_Skip
  namelist /GL_params/ Number_of_Samples
  namelist /GL_params/ Equil
  namelist /GL_params/ alpha
  namelist /GL_params/ r1
  namelist /GL_params/ r2
  namelist /GL_params/ umbrella_max
  namelist /GL_params/ umbrella_min
  namelist /GL_params/ adaptive
  namelist /GL_params/ Reverse









contains

  ! read in fortran namelist into common_params_module
  subroutine read_GL_namelist(inputs_file,length) bind(C, name="read_GL_namelist")

    integer               , value         :: length
    character(kind=c_char), intent(in   ) :: inputs_file(length)

    phi_inc = 0.d0
    n_inc_phi = 10000000

    ! read in common namelist
    open(unit=100, file=amrex_string_c_to_f(inputs_file), status='old', action='read')
    read(unit=100, nml=GL_params)
    close(unit=100)

  end subroutine read_GL_namelist

  ! copy contents of common_params_module to C++ common namespace
  subroutine initialize_GL_namespace(acoef_in,bcoef_in,ccoef_in,dcoef_in,umbrella_in,diff_coef_in,noise_coef_in) &
                                         bind(C, name="initialize_GL_namespace")


  double precision :: acoef_in,bcoef_in,ccoef_in,dcoef_in
  double precision :: umbrella_in
  double precision :: diff_coef_in
  double precision :: noise_coef_in
  double precision :: phi0_in
  double precision :: rad_in


  acoef_in = acoef
  bcoef_in = bcoef
  ccoef_in = ccoef
  dcoef_in = dcoef

  umbrella_in = umbrella
  diff_coef_in = diff_coef
  noise_coef_in = noise_coef

  phi0_in = phi0

  rad_in = rad

  end subroutine initialize_GL_namespace

end module GL_namelist_module
