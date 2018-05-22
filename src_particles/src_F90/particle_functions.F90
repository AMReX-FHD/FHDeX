module particle_functions_module

  use amrex_error_module
  use amrex_fort_module
  use common_namelist_module

  implicit none
  private

  public  f_particle
  
  type, bind(C)  :: f_particle
     real(amrex_particle_real)   :: pos(3) 
     real(amrex_particle_real)   :: vel(3)
     integer(c_int)              :: id
     integer(c_int)              :: cpu
  end type f_particle

contains

  subroutine move_particles(particles, np, dt, prob_lo, prob_hi) &
       bind(c,name='move_particles')
    
    implicit none

    integer,          intent(in   )         :: np
    type(f_particle), intent(inout), target :: particles(np)
    double precision, intent(in   )         :: dt
    double precision, intent(in   )         :: prob_lo(3), prob_hi(3)

    integer i
    type(f_particle), pointer :: p

    do i = 1, np
       
       p => particles(i)

      print *, p%id
       print *, p%pos(1), p%vel(1)

       p%pos(1) = p%pos(1) + p%vel(1) * dt 
       p%pos(2) = p%pos(2) + p%vel(2) * dt 
       p%pos(3) = p%pos(3) + p%vel(3) * dt 

    end do

  end subroutine move_particles

end module particle_functions_module

