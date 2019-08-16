module membrane_forces_module
  use amrex_fort_module, only: amrex_real
  use iso_c_binding
   
contains

! For Daniel and Andy:
! 0) I still don't get how particle initialization (positions) works but will try to get answers from Daniel in person.
! 1) Please confirm the particles are in their original order
! 2) Are the particles ordered by species also? Perhaps it is better to pass a species array here as input as well?
! 3) Is spec3xForce really intent(inout) or intent(in), that is, does this routine increment the fources or compute forces
! 4) Can we get time as an input argument of this routine also please? This should be based on the temporal integrator used, i.e., at which position in time the particle positions are.
! 5) I have added routines user_force_calc_init and user_force_calc_destroy here, which are to be called to initialize the user force code (e.g., read namelist with parameters for membrane) and to open files and then close them. We can adjust the interfaces as needed.
! 6) When the particle positions are output from inside the C++ code, are they sorted in original order?
! 7) How are namelists read in this C++/Fortran code -- see .

subroutine user_force_calc_init(unit) bind(c,name="user_force_calc_init")
   ! Read namelists, open files, etc.
   integer, intent(in) :: unit ! Unit to use to read namelist from, or maybe we want to pass file name, or what?
end subroutine user_force_calc_init

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

subroutine user_force_calc_destroy() bind(c,name="user_force_calc_destroy")
   ! Free files etc
end subroutine user_force_calc_destroy

end module membrane_forces_module
