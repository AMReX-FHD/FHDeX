module sqrt_mf_module

  implicit none

  private

contains

#if (AMREX_SPACEDIM == 1)

  subroutine sqrt_mf(lo,hi, &
                       mf, mflo, mfhi, ncomp) bind(C,name="sqrt_mf")
      
      integer         , intent(in   ) :: ncomp
      integer         , intent(in   ) :: lo(1),hi(1)
      integer         , intent(in   ) :: mflo(2),mfhi(1)
      double precision, intent(inout) :: mf(mflo(1):mfhi(1),ncomp)

      ! local
      integer :: i,j,n
      
      do n = 1,ncomp
      do i = lo(1),hi(1)
         mf(i,n) = sqrt(mf(i,n))
      end do
      end do

    end subroutine sqrt_mf

#endif

#if (AMREX_SPACEDIM == 2)

  subroutine sqrt_mf(lo,hi, &
                       mf, mflo, mfhi, ncomp) bind(C,name="sqrt_mf")
      
      integer         , intent(in   ) :: ncomp
      integer         , intent(in   ) :: lo(2),hi(2)
      integer         , intent(in   ) :: mflo(2),mfhi(2)
      double precision, intent(inout) :: mf(mflo(1):mfhi(1),mflo(2):mfhi(2),ncomp)

      ! local
      integer :: i,j,n
      
      do n = 1,ncomp
      do j = lo(2),hi(2)
      do i = lo(1),hi(1)
         mf(i,j,n) = sqrt(mf(i,j,n))
      end do
      end do
      end do

    end subroutine sqrt_mf

#endif

#if (AMREX_SPACEDIM == 3)

    subroutine sqrt_mf(lo,hi, &
                       mf, mflo, mfhi, ncomp) bind(C,name="sqrt_mf")
      
      integer         , intent(in   ) :: ncomp
      integer         , intent(in   ) :: lo(3),hi(3)
      integer         , intent(in   ) :: mflo(3),mfhi(3)
      double precision, intent(inout) :: mf(mflo(1):mfhi(1),mflo(2):mfhi(2),mflo(3):mfhi(3),ncomp)

      ! local
      integer :: i,j,k,n
     
      do n = 1,ncomp
      do k = lo(3),hi(3)
      do j = lo(2),hi(2)
      do i = lo(1),hi(1)
         mf(i,j,k,n) = sqrt(mf(i,j,k,n))
      end do
      end do
      end do
      end do
      
    end subroutine sqrt_mf
    
#endif

  end module sqrt_mf_module
