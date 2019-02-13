subroutine set_particle_properties(pstate, pradius, pdensity, pvol, pmass, omoi, omega) &
     bind(C, name="set_particle_properties")

  use amrex_fort_module, only: rt => amrex_real
  use iso_c_binding ,    only: c_int
  use constant,          only: pi
  use param,            only: zero

  implicit none

  integer(c_int), intent(in)  :: pstate
  real(rt),   intent(in)  :: pradius, pdensity
  real(rt),   intent(out) :: pvol, pmass, omoi, omega

  pvol  = (4.0d0/3.0d0)*pi*pradius**3
  pmass = pvol * pdensity
  omoi  = 2.5d0/(pmass * pradius**2)
  omega = 0.d0

end subroutine set_particle_properties
