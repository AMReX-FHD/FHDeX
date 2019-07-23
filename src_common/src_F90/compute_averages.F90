module compute_averages_module
  
  use amrex_error_module
  use amrex_fort_module, only : amrex_real
  implicit none

  private

  public :: compute_vert_average

contains

  subroutine compute_vert_average(inputfab, inlo, inhi, insize, meanfab, mlo, mhi, msize, dir, incomp, mcomp, ncomp) bind(C, name="compute_vert_average")

    integer         , intent(in   ) :: inlo(3),inhi(3), mlo(3),mhi(3), insize,msize, incomp,mcomp,ncomp, dir
    real(amrex_real), intent(in   ) :: inputfab(inlo(1):inhi(1),inlo(2):inhi(2),inlo(3):inhi(3),1:insize)
    real(amrex_real), intent(inout) :: meanfab(mlo(1):mhi(1),mlo(2):mhi(2),mlo(3):mhi(3),1:msize)

    ! local variables
    integer i,j,k,n
    double precision ninv

    ninv = 1.d0/(inhi(dir) - inlo(dir))

    SELECT CASE (dir)
    CASE (0)
       do n = 0, ncomp-1
          do k = inlo(3), inhi(3)
             do j = inlo(2), inhi(2)
                do i = inlo(1), inhi(1)
                   meanfab(mlo(1),j,k,mcomp+n) = meanfab(mlo(1),j,k,mcomp+n) + ninv*inputfab(i,j,k,incomp+n)
                end do
             end do
          end do
       end do
    CASE (1)
       do n = 0, ncomp-1
          do k = inlo(3), inhi(3)
             do j = inlo(2), inhi(2)
                do i = inlo(1), inhi(1)
                   meanfab(i,mlo(2),k,mcomp+n) = meanfab(i,mlo(2),k,mcomp+n) + ninv*inputfab(i,j,k,incomp+n)
                end do
             end do
          end do
       end do
    CASE (2)
       do n = 0, ncomp-1
          do k = inlo(3), inhi(3)
             do j = inlo(2), inhi(2)
                do i = inlo(1), inhi(1)
                   meanfab(i,j,mlo(3),mcomp+n) = meanfab(i,j,mlo(3),mcomp+n) + ninv*inputfab(i,j,k,incomp+n)
                end do
             end do
          end do
       end do
    CASE DEFAULT
       call bl_error('Invalid average direction')
    END SELECT

  end subroutine compute_vert_average

end module compute_averages_module
