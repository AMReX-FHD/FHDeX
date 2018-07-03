module compute_div_and_grad_module

  use amrex_fort_module, only : amrex_real
  implicit none

  private

  public :: compute_div, compute_grad

contains

  ! Computes divergence at cell centres from face values
#if (AMREX_SPACEDIM == 2)
subroutine compute_div(lo, hi, phix, phixlo, phixhi, phiy, phiylo, phiyhi, &
                       div, divlo, divhi, dx, increment) bind(C, name="compute_div")

  integer         , intent(in   ) :: lo(2), hi(2), phixlo(2), phixhi(2), phiylo(2), phiyhi(2)
  integer         , intent(in   ) :: divlo(2), divhi(2), increment
  real(amrex_real), intent(in   ) :: phix(phixlo(1):phixhi(1),phixlo(2):phixhi(2))
  real(amrex_real), intent(in   ) :: phiy(phiylo(1):phiyhi(1),phiylo(2):phiyhi(2))
  real(amrex_real), intent(inout) :: div(divlo(1):divhi(1),divlo(2):divhi(2))
  real(amrex_real), intent(in   ) :: dx(2)

  ! local variables
  integer i,j
  double precision dxInv, dyInv

  dxInv = 1.d0/dx(1)
  dyInv = 1.d0/dx(2)

  if (increment .eq. 0) then
     ! compute divergence
     do j = lo(2), hi(2)
     do i = lo(1), hi(1)
        div(i,j) = (-phix(i,j) + phix(i+1,j))*dxInv + &
                   (-phiy(i,j) + phiy(i,j+1))*dyInv
     end do
     end do
  else
     ! increment divergence
     do j = lo(2), hi(2)
     do i = lo(1), hi(1)
        div(i,j) = div(i,j) + (-phix(i,j) + phix(i+1,j))*dxInv + &
                              (-phiy(i,j) + phiy(i,j+1))*dyInv
     end do
     end do
  endif

end subroutine compute_div
#endif

#if (AMREX_SPACEDIM==3)
subroutine compute_div(lo, hi, phix, phixlo, phixhi, phiy, phiylo, phiyhi, phiz, phizlo, phizhi, &
                       div, divlo, divhi, dx, increment) bind(C, name="compute_div")

  integer         , intent(in   ) :: lo(3), hi(3)
  integer         , intent(in   ) :: phixlo(3), phixhi(3), phiylo(3), phiyhi(3), phizlo(3), phizhi(3)
  integer         , intent(in   ) :: divlo(3), divhi(3), increment
  real(amrex_real), intent(inout) :: phix(phixlo(1):phixhi(1),phixlo(2):phixhi(2),phixlo(3):phixhi(3))
  real(amrex_real), intent(inout) :: phiy(phiylo(1):phiyhi(1),phiylo(2):phiyhi(2),phiylo(3):phiyhi(3))
  real(amrex_real), intent(inout) :: phiz(phizlo(1):phizhi(1),phizlo(2):phizhi(2),phizlo(3):phizhi(3))
  real(amrex_real), intent(inout) :: div(divlo(1):divhi(1),divlo(2):divhi(2),divlo(3):divhi(3))
  real(amrex_real), intent(in   ) :: dx(3) 

  ! local variables
  integer i,j,k
  double precision dxInv, dyInv, dzInv

  dxInv = 1.d0/dx(1)
  dyInv = 1.d0/dx(2)
  dzInv = 1.d0/dx(3)

  if (increment .eq. 0) then
     ! compute divergence
     do k = lo(3), hi(3)
     do j = lo(2), hi(2)
     do i = lo(1), hi(1)
        div(i,j,k) = (-phix(i,j,k) + phix(i+1,j,k))*dxInv + &
                     (-phiy(i,j,k) + phiy(i,j+1,k))*dyInv + &
                     (-phiz(i,j,k) + phiy(i,j,k+1))*dzInv
     end do
     end do
     end do
  else
     ! increment divergence
     do k = lo(3), hi(3)
     do j = lo(2), hi(2)
     do i = lo(1), hi(1)
        div(i,j,k) = div(i,j,k) + (-phix(i,j,k) + phix(i+1,j,k))*dxInv + &
                                  (-phiy(i,j,k) + phiy(i,j+1,k))*dyInv + &
                                  (-phiz(i,j,k) + phiy(i,j,k+1))*dzInv
     end do
     end do
     end do
  endif

end subroutine compute_div
#endif

! Computes gradient of cell centred scalar field at cell faces
#if (AMREX_SPACEDIM==2)
subroutine compute_grad(lo, hi, gphix, gphixlo, gphixhi, gphiy, gphiylo, gphiyhi, &
                        phi, philo, phihi, dx) bind(C, name="compute_grad")

  integer         , intent(in   ) :: lo(2), hi(2), gphixlo(2), gphixhi(2), gphiylo(2), gphiyhi(2)
  integer         , intent(in   ) :: philo(2), phihi(2)
  real(amrex_real), intent(inout) :: gphix(gphixlo(1):gphixhi(1),gphixlo(2):gphixhi(2))
  real(amrex_real), intent(inout) :: gphiy(gphiylo(1):gphiyhi(1),gphiylo(2):gphiyhi(2))
  real(amrex_real), intent(in   ) :: phi(philo(1):phihi(1),philo(2):phihi(2))
  real(amrex_real), intent(in   ) :: dx(2) 

  ! local variables
  integer i,j
  double precision dxInv, dyInv;

  dxInv = 1.d0/dx(1)
  dyInv = 1.d0/dx(2)

  ! Currently assuming periodic

  ! x faces
  do j = lo(2), hi(2)
  do i = lo(1), hi(1)+1
     gphix(i,j) = (phi(i,j) - phi(i-1,j))*dxInv
  end do
  end do

  ! y faces
  do j = lo(2), hi(2)+1
  do i = lo(1), hi(1)
     gphiy(i,j) = (phi(i,j) - phi(i,j-1))*dyInv
  end do
  end do

end subroutine compute_grad
#endif

#if (AMREX_SPACEDIM==3)
subroutine compute_grad(lo, hi, gphix, gphixlo, gphixhi, gphiy, gphiylo, gphiyhi, gphiz, gphizlo, gphizhi, &
                        phi, philo, phihi, dx) bind(C, name="compute_grad")

  integer         , intent(in   ) :: lo(3), hi(3)
  integer         , intent(in   ) :: gphixlo(3), gphixhi(3), gphiylo(3), gphiyhi(3), gphizlo(3), gphizhi(3)
  integer         , intent(in   ) :: philo(3), phihi(3)
  real(amrex_real), intent(inout) :: gphix(gphixlo(1):gphixhi(1),gphixlo(2):gphixhi(2),gphixlo(3):gphixhi(3))
  real(amrex_real), intent(inout) :: gphiy(gphiylo(1):gphiyhi(1),gphiylo(2):gphiyhi(2),gphiylo(3):gphiyhi(3))
  real(amrex_real), intent(inout) :: gphiz(gphizlo(1):gphizhi(1),gphizlo(2):gphizhi(2),gphizlo(3):gphizhi(3))
  real(amrex_real), intent(in   ) :: phi(philo(1):phihi(1),philo(2):phihi(2),philo(3):phihi(3))
  real(amrex_real), intent(in   ) :: dx(3) 

  ! local variables
  integer i,j,k
  double precision dxInv, dyInv, dzInv

  dxInv = 1.d0/dx(1)
  dyInv = 1.d0/dx(2)
  dzInv = 1.d0/dx(3)

  ! Currently assuming periodic

  ! x faces
  do k = lo(3), hi(3)
  do j = lo(2), hi(2)
  do i = lo(1), hi(1)+1
     gphix(i,j,k) = (phi(i,j,k) - phi(i-1,j,k))*dxInv
  end do
  end do
  end do

  ! y faces
  do k = lo(3), hi(3)
  do j = lo(2), hi(2)+1
  do i = lo(1), hi(1)
     gphiy(i,j,k) = (phi(i,j,k) - phi(i,j-1,k))*dyInv
  end do
  end do
  end do

  ! z faces
  do k = lo(3), hi(3)+1
  do j = lo(2), hi(2)
  do i = lo(1), hi(1)
     gphiy(i,j,k) = (phi(i,j,k) - phi(i,j,k-1))*dzInv
  end do
  end do
  end do

end subroutine compute_grad
#endif

end module compute_div_and_grad_module
