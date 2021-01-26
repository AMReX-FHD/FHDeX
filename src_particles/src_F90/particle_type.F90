module cell_sorted_particle_module
  use amrex_fort_module, only: amrex_real, amrex_particle_real
  use iso_c_binding ,    only: c_int
  
  implicit none
  private
  
  public particle_t, dsmc_particle_t, remove_particle_from_cell
  
  type, bind(C) :: particle_t
#if (BL_SPACEDIM == 3)
     real(amrex_particle_real) :: pos(3)
#endif
#if (BL_SPACEDIM == 2)
     real(amrex_particle_real) :: pos(2)
#endif
     real(amrex_particle_real) :: radius
     real(amrex_particle_real) :: vel(3)
     real(amrex_particle_real) :: force(3)
     real(amrex_particle_real) :: posP(3)
     real(amrex_particle_real) :: velP(3)
     real(amrex_particle_real) :: forceP(3)

     real(amrex_particle_real) :: velO(3)     
     real(amrex_particle_real) :: forceO(3)
     real(amrex_particle_real) :: localvel(3)
     real(amrex_particle_real) :: mass
     real(amrex_particle_real) :: R
     real(amrex_particle_real) :: q
     real(amrex_particle_real) :: accel_factor
     real(amrex_particle_real) :: drag_factor

     real(amrex_particle_real) :: origin(3)
     real(amrex_particle_real) :: abspos(3)
     real(amrex_particle_real) :: travel_time
     real(amrex_particle_real) :: diff_av
     real(amrex_particle_real) :: step_count
     real(amrex_particle_real) :: multi
     real(amrex_particle_real) :: dry_diff
     real(amrex_particle_real) :: wet_diff
     real(amrex_particle_real) :: total_diff
     real(amrex_particle_real) :: sigma
     real(amrex_particle_real) :: eepsilon
     real(amrex_particle_real) :: potential
     real(amrex_particle_real) :: p3m_radius
     real(amrex_particle_real) :: spring

     integer(c_int)            :: id         
     integer(c_int)            :: cpu        
     integer(c_int)            :: sorted     
     integer(c_int)            :: i
     integer(c_int)            :: j
     integer(c_int)            :: k
     integer(c_int)            :: species
     integer(c_int)            :: visible
     integer(c_int)            :: pinned
  end type particle_t

  type, bind(C) :: dsmc_particle_t
#if (BL_SPACEDIM == 3)
     real(amrex_particle_real) :: pos(3)
#endif
#if (BL_SPACEDIM == 2)
     real(amrex_particle_real) :: pos(2)
#endif
     real(amrex_particle_real) :: radius
     real(amrex_particle_real) :: velJ(3)
     real(amrex_particle_real) :: forceJ(3)
     real(amrex_particle_real) :: posP(3)
     real(amrex_particle_real) :: velP(3)
     real(amrex_particle_real) :: forceP(3)

     real(amrex_particle_real) :: vel(3)     
     real(amrex_particle_real) :: force(3)
     real(amrex_particle_real) :: localvel(3)
     real(amrex_particle_real) :: mass
     real(amrex_particle_real) :: R
     real(amrex_particle_real) :: q
     real(amrex_particle_real) :: accel_factor
     real(amrex_particle_real) :: drag_factor

     real(amrex_particle_real) :: origin(3)
     real(amrex_particle_real) :: abspos(3)
     real(amrex_particle_real) :: travel_time
     real(amrex_particle_real) :: diff_av
     real(amrex_particle_real) :: step_count
     real(amrex_particle_real) :: multi
     real(amrex_particle_real) :: dry_diff
     real(amrex_particle_real) :: wet_diff
     real(amrex_particle_real) :: total_diff
     real(amrex_particle_real) :: sigma
     real(amrex_particle_real) :: eepsilon
     real(amrex_particle_real) :: potential
     real(amrex_particle_real) :: p3m_radius

     integer(c_int)            :: id         
     integer(c_int)            :: cpu        
     integer(c_int)            :: sorted     
     integer(c_int)            :: i
     integer(c_int)            :: j
     integer(c_int)            :: k
     integer(c_int)            :: species
  end type dsmc_particle_t

contains
  
  subroutine remove_particle_from_cell(cell_parts, cell_np, new_np, i) &
                             bind(c,name="remove_particle_from_cell")
    
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
