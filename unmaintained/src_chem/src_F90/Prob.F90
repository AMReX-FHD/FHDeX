subroutine initdata(level, time, lo, hi, &
     con, con_lo, con_hi, &
     dx, prob_lo) bind(C, name="initdata")

  use amrex_fort_module, only : amrex_spacedim

  implicit none
  integer, intent(in) :: level, lo(3), hi(3), con_lo(3), con_hi(3)
  double precision, intent(in) :: time
  double precision, intent(inout) :: con(con_lo(1):con_hi(1), &
       &                                 con_lo(2):con_hi(2), &
       &                                 con_lo(3):con_hi(3))
  double precision, intent(in) :: dx(3), prob_lo(3)

  integer          :: i,j,k
  double precision :: x,y,z,r2
 ! Use for initial condition for concentration
  do k=lo(3),hi(3)
     do j=lo(2),hi(2)
        do i=lo(1),hi(1)
              con(i,j,k) = 0.d0! + exp(-r2)
        end do
     end do
  end do

end subroutine initdata
