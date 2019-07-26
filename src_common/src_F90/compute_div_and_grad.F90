module compute_div_and_grad_module

  use amrex_fort_module, only : amrex_real
  implicit none

  private

contains

! Computes gradient of cell centred scalar field at cell centre
#if (AMREX_SPACEDIM==2)
subroutine compute_grad_cc(lo, hi, gphix, gphixlo, gphixhi, gphiy, gphiylo, gphiyhi, &
                        phi, philo, phihi, dx) bind(C, name="compute_grad_cc")

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

end subroutine compute_grad_cc
#endif

#if (AMREX_SPACEDIM==3)
subroutine compute_grad_cc(lo, hi, gphix, gphixlo, gphixhi, gphiy, gphiylo, gphiyhi, gphiz, gphizlo, gphizhi, &
                        phi, philo, phihi, dx) bind(C, name="compute_grad_cc")

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
  double precision dxInv, dyInv, dzInv, dxInv2, dyInv2, dzInv2

  dxInv = 1.d0/dx(1)
  dyInv = 1.d0/dx(2)
  dzInv = 1.d0/dx(3)

  ! x faces
  do k = lo(3), hi(3)
  do j = lo(2), hi(2)
  do i = lo(1), hi(1)
     gphix(i,j,k) = (phi(i+1,j,k) - phi(i-1,j,k))*dxInv
  end do
  end do
  end do

  ! y faces
  do k = lo(3), hi(3)
  do j = lo(2), hi(2)
  do i = lo(1), hi(1)
     gphiy(i,j,k) = (phi(i,j+1,k) - phi(i,j-1,k))*dyInv

  end do
  end do
  end do

  ! z faces
  do k = lo(3), hi(3)
  do j = lo(2), hi(2)
  do i = lo(1), hi(1)
     gphiz(i,j,k) = (phi(i,j,k+1) - phi(i,j,k-1))*dzInv
  end do
  end do
  end do

end subroutine compute_grad_cc
#endif

end module compute_div_and_grad_module
