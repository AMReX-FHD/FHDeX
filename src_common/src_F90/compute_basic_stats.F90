module compute_basic_stats_module

  use amrex_fort_module, only : amrex_real
  implicit none

  private

  public :: compute_means, compute_vars, sum_fab

contains

subroutine compute_means(instfab, inlo, inhi, insize, meanfab, mlo, mhi, msize, incomp, mcomp, steps) bind(C, name="compute_means")

  integer         , intent(in   ) :: inlo(3), inhi(3), mlo(3), mhi(3), insize, msize, incomp, mcomp, steps
  real(amrex_real), intent(in   ) :: instfab(inlo(1):inhi(1),inlo(2):inhi(2),inlo(3):inhi(3),1:insize)
  real(amrex_real), intent(inout) :: meanfab(mlo(1):mhi(1),mlo(2):mhi(2),mlo(3):mhi(3),1:msize)

  ! local variables
  integer i,j, k
  double precision stepsinv, stepsminusone

  stepsinv = 1d0/steps
  stepsminusone = steps - 1

  do k = mlo(3), mhi(3)
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

  stepsinv = 1d0/steps
  stepsminusone = steps - 1

  do k = vlo(3), vhi(3)
    do j = vlo(2), vhi(2)
      do i = vlo(1), vhi(1)

        del = instfab(i,j,k,incomp) - meanfab(i,j,k,mcomp)

        varfab(i,j,k,vcomp) = (varfab(i,j,k,vcomp)*stepsminusone + del**2)*stepsinv

      end do
    end do
  end do

end subroutine compute_vars

subroutine sum_fab(lo, hi, infab, inlo, inhi, insize, gs, total, comp) bind(C, name="sum_fab")

  integer         , intent(in   ) :: inlo(3), inhi(3), lo(3), hi(3), insize, gs, comp
  real(amrex_real), intent(inout) :: total

  double precision, intent(in   ) :: infab(inlo(1):inhi(1),inlo(2):inhi(2),inlo(3):inhi(3), 1:insize)

  ! local variables
  integer i,j, k
  
  total = 0

  !print *, lo, hi, gs

  do k = lo(3) - gs, hi(3) + gs
    do j = lo(2) - gs, hi(2) + gs
      do i = lo(1) - gs, hi(1) + gs

        total = total + infab(i,j,k,comp)
       
      end do
    end do
  end do

end subroutine sum_fab

subroutine x_mean_fab(lo, hi, infab, inlo, inhi, insize, outfab, outlo, outhi, outsize, gs) bind(C, name="x_mean_fab")

  integer         , intent(in   ) :: outlo(3), outhi(3), inlo(3), inhi(3), lo(3), hi(3), insize, outsize, gs

  double precision, intent(in   ) :: infab(inlo(1):inhi(1),inlo(2):inhi(2),inlo(3):inhi(3), 1:insize)
  double precision, intent(inout) :: outfab(outlo(1):outhi(1),outlo(2):outhi(2),outlo(3):outhi(3), 1:outsize)

  ! local variables
  integer i,j, k, l, cellcount
  double precision xmean
 
  !print *, lo, hi, gs
  do l = 1, insize
    do k = lo(3) - gs, hi(3) + gs
      do j = lo(2) - gs, hi(2) + gs

        xmean = 0
        cellcount = 0
        do i = lo(1) - gs, hi(1) + gs

          xmean = xmean + infab(i,j,k,l)
          cellcount = cellcount + 1
         
        end do

        do i = lo(1) - gs, hi(1) + gs

          outfab(i,j,k,l) = xmean/cellcount
         
        end do
      end do
    end do
  end do

end subroutine x_mean_fab

end module compute_basic_stats_module
