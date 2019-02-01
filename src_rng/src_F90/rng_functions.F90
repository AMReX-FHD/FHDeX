module rng_functions_module

  use amrex_error_module
  use amrex_fort_module
  use amrex_constants_module
  use common_namelist_module

  use bl_types
  use bl_random_module
  use parallel
  use bl_error_module


  implicit none
  private

  public :: rng_initialize, &
       &    get_particle_normal, get_particle_normal_func, get_selector, get_uniform, get_uniform_func, &
       &    get_angles, get_half_angles, get_fhd_normal_func

  type(bl_rng_engine)      , save :: rng_eng_fhd
  type(bl_rng_engine)      , save :: rng_eng_particle
  type(bl_rng_engine)      , save :: rng_eng_select
  type(bl_rng_engine)      , save :: rng_eng_scatter_theta
  type(bl_rng_engine)      , save :: rng_eng_scatter_phi
  type(bl_rng_engine)      , save :: rng_eng_general

  type(bl_rng_normal)      , save :: nm_fhd
  type(bl_rng_normal)      , save :: nm_particle
  type(bl_rng_uniform_real), save :: un_select
  type(bl_rng_uniform_real), save :: un_costheta
  type(bl_rng_uniform_real), save :: un_phi
  type(bl_rng_uniform_real), save :: un_general

contains

  subroutine rng_initialize(fhd, particle, sel, theta, phi, general) bind(c,name='rng_initialize')

    implicit none

    integer, intent(in) :: fhd, particle, sel, theta, phi, general

    call bl_rng_build_engine(rng_eng_fhd, fhd)
    call bl_rng_build_engine(rng_eng_particle, particle)
    call bl_rng_build_engine(rng_eng_select, sel)
    call bl_rng_build_engine(rng_eng_scatter_theta, theta)
    call bl_rng_build_engine(rng_eng_scatter_phi, phi)
    call bl_rng_build_engine(rng_eng_general, general)

    call bl_rng_build_distro(nm_fhd, 0.0d0, 1.0d0)
    call bl_rng_build_distro(nm_particle, 0.0d0, 1.0d0)
    call bl_rng_build_distro(un_select, 0.0d0, 1.0d0)
    call bl_rng_build_distro(un_costheta, 0.0d0, 1.0d0)
    call bl_rng_build_distro(un_general, 0.0d0, 1.0d0)
    call bl_rng_build_distro(un_phi, 0.0d0, 2d0*3.14159265359)



  end subroutine rng_initialize


  subroutine get_particle_normal(test)

      double precision, intent(inout) :: test

      test = bl_rng_get(nm_particle, rng_eng_particle)

  end subroutine get_particle_normal

  function get_particle_normal_func() result(test) bind(c,name='get_particle_normal_func')

      double precision test

      test = bl_rng_get(nm_particle, rng_eng_particle)

  end function get_particle_normal_func

  subroutine get_selector(test, length)

      integer, intent(inout) :: test
      integer, intent(in)    :: length

      test = ceiling(bl_rng_get(un_select, rng_eng_select)*length)

  end subroutine get_selector

  subroutine get_angles(costheta, sintheta, cosphi, sinphi) bind(c,name="get_angles")

      double precision, intent(inout) :: costheta, sintheta, cosphi, sinphi
      double precision phi

      costheta = 2d0*bl_rng_get(un_costheta, rng_eng_scatter_theta) - 1d0
      sintheta = sqrt(1d0 - costheta**2)

      phi = bl_rng_get(un_phi, rng_eng_scatter_phi)
      cosphi = cos(phi)
      sinphi = sin(phi)

  end subroutine get_angles

  subroutine get_half_angles(costheta, sintheta, cosphi, sinphi)

      double precision, intent(inout) :: costheta, sintheta, cosphi, sinphi
      double precision phi

      costheta = bl_rng_get(un_costheta, rng_eng_scatter_theta)
      sintheta = sqrt(1d0 - costheta**2)

      phi = bl_rng_get(un_phi, rng_eng_scatter_phi)
      cosphi = cos(phi)
      sinphi = sin(phi)

  end subroutine get_half_angles

  subroutine get_uniform(test)

      double precision, intent(inout) :: test

      test = bl_rng_get(un_general, rng_eng_general)

  end subroutine get_uniform

  function get_uniform_func() result(test) bind(c,name='get_uniform_func')

      double precision test

      test = bl_rng_get(un_general, rng_eng_general)

  end function get_uniform_func

  function get_fhd_normal_func() result(test) bind(c,name='get_fhd_normal_func')

      double precision test

      test = bl_rng_get(nm_fhd, rng_eng_fhd)

  end function get_fhd_normal_func

end module rng_functions_module
