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
            get_particle_normal, &
            rng_eng_fhd, &
            rng_eng_particle

  type(bl_rng_engine)      , save :: rng_eng_fhd
  type(bl_rng_engine)      , save :: rng_eng_particle

  type(bl_rng_normal)      , save :: nm_fhd
  type(bl_rng_normal)      , save :: nm_particle

contains

  subroutine rng_initialize() bind(c,name='rng_initialize')

    implicit none

    double precision :: test

    call bl_rng_build_engine(rng_eng_fhd, 2)
    call bl_rng_build_engine(rng_eng_particle, 2)       

    call bl_rng_build_distro(nm_fhd, 0.0d0, 1.0d0)
    call bl_rng_build_distro(nm_particle, 0.0d0, 1.0d0)

  end subroutine rng_initialize


  subroutine get_particle_normal(test)

      double precision, intent(inout) :: test

      test = bl_rng_get(nm_particle, rng_eng_particle)

  end subroutine get_particle_normal
  
end module rng_functions_module































