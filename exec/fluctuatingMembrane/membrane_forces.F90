module membrane_forces_module
  use amrex_fort_module, only: amrex_real
  use iso_c_binding
   
contains

! Daniel: Please confirm the following things:
! 1) The particles are in their original orders
! 2) Are the particles ordered by species also? Perhaps it is better to pass a species array here as input as well?
! 3) Is spec3xForce really intent(inout) or intent(in), that is, does this routine increment the fources or compute forces
! 4) Can we get time as an input argument of this routine also please? This should be based on the temporal integrator used, i.e., at which position in time the particle positions are.

subroutine user_force_calc(spec3xPos, spec3yPos, spec3zPos, spec3xForce, spec3yForce, spec3zForce, length) bind(c,name="user_force_calc")
  ! This routine should increment the forces supplied here

  
  implicit none

  integer(c_int),          intent(in   )         :: length
  real(c_double),          intent(in   )         :: spec3xPos(length), spec3yPos(length), spec3zPos(length)
  real(c_double),          intent(inout)         :: spec3xForce(length), spec3yForce(length), spec3zForce(length)

  integer i

  do i = 1, length

    write(*,*) "Setting force on particle=", i, " pos: ", spec3xPos(i), spec3yPos(i), spec3zPos(i)
    spec3xForce(i)=1.0d0

  enddo
  
end subroutine user_force_calc

end module membrane_forces_module
