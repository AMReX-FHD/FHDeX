subroutine init_vel(lo, hi, vel, vello, velhi, dx, prob_lo, prob_hi, di, &
                    reallo, realhi) bind(C, name="init_vel")

  use amrex_fort_module, only : amrex_real

  implicit none

  integer         , intent(in   ) :: lo(3), hi(3), vello(3), velhi(3), di
  real(amrex_real), intent(inout) :: vel(vello(1):velhi(1),vello(2):velhi(2),vello(3):velhi(3))
  real(amrex_real), intent(in   ) :: reallo(3), realhi(3)
  real(amrex_real), intent(in   ) :: prob_lo(3)
  real(amrex_real), intent(in   ) :: prob_hi(3)
  real(amrex_real), intent(in   ) :: dx(3)

  integer          :: i,j,k
  double precision :: pos(3),center(3),partdom,itVec(3),relpos(3),rad,zshft

#if (AMREX_SPACEDIM == 2)
  zshft = 0.0d0
#elif (AMREX_SPACEDIM == 3)
  zshft = 0.5d0
#endif

  center = (realhi - reallo)/2d0;
  partdom = ((realhi(1) - reallo(1))/4d0)**2;

  if (di .EQ. 0) then
     do k = lo(3), hi(3)
        do j = lo(2), hi(2)
           do i = lo(1), hi(1) + 1

              itVec(1) = dble(i)*dx(1)
              itVec(2) = (dble(j)+0.5d0)*dx(2)
              itVec(3) = (dble(k)+zshft)*dx(3)

              pos = reallo + itVec
              relpos = pos - center
              rad = DOT_PRODUCT(relpos,relpos)

              ! if (rad .LT. partdom) then
              !    vel(i,j,k) = -200*exp(-rad/(10*partdom*partdom))*relpos(2)
              ! else
              !    vel(i,j,k) = 0d0
              ! endif

              vel(i,j,k) = 1.d0 + exp(-rad**2/(1.0d0*partdom*partdom))

           end do
        end do
     end do
  endif

  if (di .EQ. 1) then
     do k = lo(3), hi(3)
        do j = lo(2), hi(2) + 1
           do i = lo(1), hi(1)

              itVec(1) = (dble(i)+0.5d0)*dx(1)
              itVec(2) = dble(j)*dx(2)
              itVec(3) = (dble(k)+zshft)*dx(3)

              pos = reallo + itVec
              relpos = pos - center
              rad = DOT_PRODUCT(relpos,relpos)

              ! if (rad .LT. partdom) then
              !    vel(i,j,k) = -200*exp(-rad/(10*partdom*partdom))*relpos(1)
              ! else
              !    vel(i,j,k) = 0
              ! endif

              vel(i,j,k) = 1.d0 + exp(-rad**2/(1.0d0*partdom*partdom))

           end do
        end do
     end do
  endif

  if (di .EQ. 2) then
     do k = lo(3), hi(3) + 1
        do j = lo(2), hi(2)
           do i = lo(1), hi(1)

              itVec(1) = (dble(i)+0.5d0)*dx(1)
              itVec(2) = (dble(j)+0.5d0)*dx(2)
              itVec(3) = dble(k)*dx(3)

              pos = reallo + itVec
              relpos = pos - center
              rad = DOT_PRODUCT(relpos,relpos)

              ! if (rad .LT. partdom) then
              !    vel(i,j,k) = 0
              ! else
              !    vel(i,j,k) = 0
              ! endif

              vel(i,j,k) = 1.d0 + exp(-rad**2/(1.0d0*partdom*partdom))

           end do
        end do
     end do
  endif


end subroutine init_vel


