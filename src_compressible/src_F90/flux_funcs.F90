module flux_module

  use amrex_fort_module, only : amrex_real
  implicit none

  private

  public :: diff_flux

contains

  subroutine diff_flux(lo,hi, cons, prim, eta, zeta, kappa, xflux, yflux, zflux, dx, nvars, nprimvars) bind(C,name="diff_flux")

      integer         , intent(in   ) :: lo(3),hi(3), nvars, nprimvars
      real(amrex_real), intent(in   ) :: dx(3)
      real(amrex_real), intent(inout) :: xflux(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3), nvars)
      real(amrex_real), intent(inout) :: yflux(lo(1):hi(1),lo(2):hi(2)+1,lo(3):hi(3), nvars)
#if (AMREX_SPACEDIM == 3)
      real(amrex_real), intent(inout) :: zflux(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)+1, nvars)
#endif

      real(amrex_real), intent(in   ) :: cons(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3), nvars)
      real(amrex_real), intent(in   ) :: prim(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3), nprimvars)

      real(amrex_real), intent(in   ) :: eta(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
      real(amrex_real), intent(in   ) :: zeta(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
      real(amrex_real), intent(in   ) :: kappa(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))

      integer :: i,j,k
      
      do k = lo(3),hi(3)
        do j = lo(2),hi(2)
          do i = lo(1),hi(1)

          end do
        end do
      end do

  end subroutine diff_flux
  
  subroutine hyp_flux(lo,hi, cons, prim, eta, zeta, kappa, xflux, yflux, zflux, dx, nvars, nprimvars) bind(C,name="hyp_flux")

      integer         , intent(in   ) :: lo(3),hi(3), nvars, nprimvars
      real(amrex_real), intent(in   ) :: dx(3)
      real(amrex_real), intent(inout) :: xflux(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3), nvars)
      real(amrex_real), intent(inout) :: yflux(lo(1):hi(1),lo(2):hi(2)+1,lo(3):hi(3), nvars)
#if (AMREX_SPACEDIM == 3)
      real(amrex_real), intent(inout) :: zflux(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)+1, nvars)
#endif

      real(amrex_real), intent(in   ) :: cons(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3), nvars)
      real(amrex_real), intent(in   ) :: prim(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3), nprimvars)

      real(amrex_real), intent(in   ) :: eta(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
      real(amrex_real), intent(in   ) :: zeta(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
      real(amrex_real), intent(in   ) :: kappa(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))

      integer :: i,j,k
      
      do k = lo(3),hi(3)
        do j = lo(2),hi(2)
          do i = lo(1),hi(1)

          end do
        end do
      end do

  end subroutine hyp_flux


end module flux_module
