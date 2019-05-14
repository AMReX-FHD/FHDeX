module compute_basic_stats_module

  use amrex_fort_module, only : amrex_real
  implicit none

  private

  public :: compute_means, compute_vars

contains

subroutine compute_means(instfab, inlo, inhi, insize, meanfab, mlo, mhi, msize, incomp, mcomp, steps) bind(C, name="compute_means")

  integer         , intent(in   ) :: inlo(3), inhi(3), mlo(3), mhi(3), insize, msize, incomp, mcomp, steps
  real(amrex_real), intent(in   ) :: instfab(inlo(1):inhi(1),inlo(2):inhi(2),inlo(3):inhi(3),1:insize)
  real(amrex_real), intent(inout) :: meanfab(mlo(1):mhi(1),mlo(2):mhi(2),mlo(3):mhi(3),1:msize)

  ! local variables
  integer i,j, k
  double precision stepsinv, stepsminusone

  stepsinv = 1/steps
  stepsminusone = steps - 1

  do k = mlo(2), mhi(2)
    do j = mlo(2), mhi(2)
      do i = mlo(1), mhi(1)

        meanfab(i,j,k,mcomp) = (meanfab(i,j,k,mcomp)*stepsminusone + instfab(i,j,k,incomp))*stepsinv

      end do
    end do
  end do

end subroutine compute_means

subroutine compute_vars(instfab, inlo, inhi, insize, meanfab, mlo, mhi, msize, varfab, vlo, vhi, vsize, incomp, mcomp, vcomp, steps) bind(C, name="compute_vars")

  integer         , intent(in   ) :: inlo(3), inhi(3), mlo(3), mhi(3), vlo(3), vhi(3), insize, msize, vsize, incomp, mcomp, vcomp, steps
  real(amrex_real), intent(in   ) :: instfab(inlo(1):inhi(1),inlo(2):inhi(2),inlo(3):inhi(3),1:insize)
  real(amrex_real), intent(in   ) :: meanfab(mlo(1):mhi(1),mlo(2):mhi(2),mlo(3):mhi(3),1:msize)
  real(amrex_real), intent(inout) :: varfab(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3),1:vsize)

  ! local variables
  integer i,j, k
  double precision stepsinv, stepsminusone, del

  stepsinv = 1/steps
  stepsminusone = steps - 1

  do k = vlo(2), vhi(2)
    do j = vlo(2), vhi(2)
      do i = vlo(1), vhi(1)

        del = instfab(i,j,k,incomp) - meanfab(i,j,k,7)

        varfab(i,j,k,vcomp) = (varfab(i,j,k,vcomp)*stepsminusone + del**2)*stepsinv

      end do
    end do
  end do

end subroutine compute_vars

end module compute_basic_stats_module
