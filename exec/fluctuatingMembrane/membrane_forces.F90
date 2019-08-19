module membrane_forces_module
  use iso_c_binding
  use amrex_fort_module, only: amrex_real
  use amrex_string_module, only: amrex_string_c_to_f, amrex_string_f_to_c
  use common_namelist_module
  use HydroGridModule
  public
  
  integer, save :: nmem=100 ! Make the membrane out of nmem^2 markers
  integer, save :: call_hydroGrid=0 ! If positive, how many steps between calls to HydroGrid

  ! For statistical analysis of fluctuating fields
  type (HydroGrid), target, save :: hydro_grid  
   
contains

! For Daniel and Andy:
! 0) Seems not to like boxes with no particles --- should be fixed. -Fixed DRL
! 1) Please confirm the particles are in their original order: CHANGE: Please pass all particles here ordered
! 2) Are the particles ordered by species also? Perhaps it is better to pass a species array here as input as well?
! 4) Can we get time AND time step as an input argument of this routine also please? This should be based on the temporal integrator used, i.e., at which position in time the particle positions are.
! 5) I have added routines user_force_calc_init and user_force_calc_destroy here, which are to be called to initialize the user force code (e.g., read namelist with parameters for membrane) and to open files and then close them. We can adjust the interfaces as needed.
! 6) When the particle positions are output from inside the C++ code, are they sorted in original order? It is only in binary in paraview

subroutine user_force_calc_init(inputs_file,length) bind(c,name="user_force_calc_init")
   ! Read namelists, open files, etc.
   integer(c_int), intent(in) :: length
   character(kind=c_char), intent(in   ) :: inputs_file(length)
   ! CHANGE: Call this routine after all the other stuff has been initialized but before time loop
   ! CHANGE: How do I get access to dx/dy/dz and dt?
   
   namelist /FluctuatingMembrane/ Nmem, call_hydroGrid
   
   open(unit=100, file=amrex_string_c_to_f(inputs_file), status='old', action='read')
   read(unit=100, nml=FluctuatingMembrane)
   close(unit=100)
   
   if(call_hydroGrid>0) then
      call createHydroAnalysis (hydro_grid, &
         nCells=(/nmem,nmem,1/), nSpecies = 1, &
         isSingleFluid = .true., nVelocityDimensions = 2, nPassiveScalars = 0, &
         systemLength = (/prob_hi(1:2)-prob_lo(1:2), 1.0d0/), timestep = call_hydroGrid*fixed_dt)   
   end if

end subroutine user_force_calc_init

subroutine user_force_calc(spec3xPos, spec3yPos, spec3zPos, spec3xForce, spec3yForce, spec3zForce, length) bind(c,name="user_force_calc")
  ! This routine should increment the forces supplied here
  ! CHANGE: Make length an array of length nspecies
  ! CHANGE: Add time and time_step as argument
  
  implicit none

  integer(c_int),          intent(in   )         :: length
  real(c_double),          intent(in   )         :: spec3xPos(length), spec3yPos(length), spec3zPos(length)
  real(c_double),          intent(out  )         :: spec3xForce(length), spec3yForce(length), spec3zForce(length)

  integer :: i


  spec3xForce=0.0d0
  spec3yForce=0.0d0
  spec3zForce=0.0d0

  do i = 1, length

    write(*,*) "Setting force on particle=", i, " pos: ", spec3xPos(i), spec3yPos(i), spec3zPos(i)

  enddo
  spec3xForce(1)=1.0d-20
  if(length>1) spec3xForce(2)=-spec3xForce(1) ! Equal and opposite to ensure solvability
  
end subroutine user_force_calc

subroutine user_force_calc_destroy() bind(c,name="user_force_calc_destroy")
   ! Free files etc
   if(call_hydroGrid>0) then
      call writeToFiles (hydro_grid)
      call destroyHydroAnalysis (hydro_grid)
   end if
end subroutine user_force_calc_destroy

end module membrane_forces_module
