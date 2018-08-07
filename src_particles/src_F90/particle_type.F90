module cell_sorted_particle_module
  use amrex_fort_module, only: amrex_real, amrex_particle_real
  use iso_c_binding ,    only: c_int
  
  implicit none
  private
  
  public particle_t, remove_particle_from_cell
  
  type, bind(C) :: particle_t
#if (BL_SPACEDIM == 3)
     real(amrex_particle_real) :: pos(3)
#endif
#if (BL_SPACEDIM == 2)
     real(amrex_particle_real) :: pos(2)
#endif


     real(amrex_particle_real) :: vel(3)     
     real(amrex_particle_real) :: mass
     real(amrex_particle_real) :: R
     real(amrex_particle_real) :: radius
     real(amrex_particle_real) :: accel_factor
     real(amrex_particle_real) :: drag_factor
     real(amrex_particle_real) :: angular_vel(3)

     real(amrex_particle_real) :: dir(3)
     real(amrex_particle_real) :: propulsion


     integer(c_int)            :: id         
     integer(c_int)            :: cpu        
     integer(c_int)            :: sorted     
     integer(c_int)            :: i
     integer(c_int)            :: j
     integer(c_int)            :: k
  end type particle_t

contains
  
  subroutine remove_particle_from_cell(cell_parts, cell_np, new_np, i)
    
    use iso_c_binding, only: c_int
    
    implicit none
    
    integer(c_int), intent(inout) :: cell_parts(cell_np)
    integer(c_int), intent(in   ) :: cell_np
    integer(c_int), intent(inout) :: new_np
    integer(c_int), intent(in   ) :: i 

    cell_parts(i) = cell_parts(new_np)
    new_np = new_np - 1
        
  end subroutine remove_particle_from_cell
  
end module cell_sorted_particle_module
