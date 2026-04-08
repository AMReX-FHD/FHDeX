subroutine init_vel(lo, hi, vel, vello, velhi, dx, prob_lo, prob_hi, di, reallo, realhi) bind(C, name="init_vel")

  use amrex_fort_module, only : amrex_real

  implicit none

  integer         , intent(in   ) :: lo(3), hi(3), vello(3), velhi(3), di
  real(amrex_real), intent(inout) :: vel(vello(1):velhi(1),vello(2):velhi(2),vello(3):velhi(3))
  real(amrex_real), intent(in   ) :: reallo(3), realhi(3)
  real(amrex_real), intent(in   ) :: prob_lo(3)
  real(amrex_real), intent(in   ) :: prob_hi(3)
  real(amrex_real), intent(in   ) :: dx(3)

  integer          :: i,j,k
  double precision :: pos(3),center(3),partdom,itVec(3),relpos(3),rad

  center = (realhi - reallo)/2d0;
  partdom = ((realhi(1) - reallo(1))/2d0)**2;

  print *, partdom, realhi(1), reallo(1)

  if (di .EQ. 0) then
     do k = lo(3), hi(3)
        do j = lo(2), hi(2)
           do i = lo(1), hi(1) + 1

              itVec(1) = dble(i)*dx(1)
              itVec(2) = dble(j)*dx(2)
              itVec(3) = dble(k)*dx(3)

              pos = reallo + itVec
              relpos = pos - center
              rad = DOT_PRODUCT(relpos,relpos)

              !print *, "rad: ", rad, "partdom: ", partdom

              if (rad .LT. partdom) then
                 vel(i,j,k) = -200*exp(-rad/(10*partdom*partdom))*relpos(2)
              else
                 vel(i,j,k) = 0d0
              endif

           end do
        end do
     end do
  endif

  if (di .EQ. 1) then
     do k = lo(3), hi(3)
        do j = lo(2), hi(2) + 1
           do i = lo(1), hi(1)

              itVec(1) = dble(i)*dx(1)
              itVec(2) = dble(j)*dx(2)
              itVec(3) = dble(k)*dx(3)

              pos = reallo + itVec
              relpos = pos - center
              rad = DOT_PRODUCT(relpos,relpos)

              if (rad .LT. partdom) then
                 vel(i,j,k) = 100*exp(-rad/(10*partdom*partdom))*relpos(1)
              else
                 vel(i,j,k) = 0
              endif

           end do
        end do
     end do
  endif

  if (di .EQ. 2) then
     do k = lo(3), hi(3) + 1
        do j = lo(2), hi(2)
           do i = lo(1), hi(1)

              itVec(1) = dble(i)*dx(1)
              itVec(2) = dble(j)*dx(2)
              itVec(3) = dble(k)*dx(3)

              pos = reallo + itVec
              relpos = pos - center
              rad = DOT_PRODUCT(relpos,relpos)

              if (rad .LT. partdom) then
                 vel(i,j,k) = 0
              else
                 vel(i,j,k) = 0
              endif

           end do
        end do
     end do
  endif


end subroutine init_vel


