module multifab_fill_random_module

  use amrex_error_module
  use rng_functions_module

  implicit none

  private

contains

#if (AMREX_SPACEDIM == 2)

    subroutine multifab_fill_random(lo,hi, &
                                    mf,mflo,mfhi, &
                                    ncomp,comp) bind(C,name="multifab_fill_random")

      integer         , intent(in   ) :: lo(2),hi(2)
      integer         , intent(in   ) :: mflo(2),mfhi(2)
      double precision, intent(inout) :: mf(mflo(1):mfhi(1),mflo(2):mfhi(2),ncomp)
      integer         , intent(in   ) :: ncomp,comp
      
      ! local
      integer :: i,j

      !=============================
      ! fill elements of array with random numbers
      !=============================

      ! print*, "F90 hack: comp = ", comp, "/", ncomp

      do j=lo(2),hi(2)
      do i=lo(1),hi(1)
         mf(i,j,comp+1) = get_fhd_normal_func()
      end do
      end do

      !! Hack:
      ! print*, "box lo & hi: ", lo, hi
      ! print*, "mf lo & hi: ", mflo, mfhi
      ! print*, "mf x-lo bound: ", mf(mflo(1),mflo(2):mfhi(2))
      ! print*, "mf x-hi bound: ", mf(mfhi(1),mflo(2):mfhi(2))
      ! print*, "mf y-lo bound: ", mf(mflo(1):mfhi(1),mflo(2))
      ! print*, "mf y-hi bound: ", mf(mflo(1):mfhi(1),mfhi(2))

    end subroutine multifab_fill_random

#endif

#if (AMREX_SPACEDIM == 3)

    subroutine multifab_fill_random(lo,hi, &
                                    mf,mflo,mfhi, &
                                    ncomp,comp) bind(C,name="multifab_fill_random")

      integer         , intent(in   ) :: lo(3),hi(3)
      integer         , intent(in   ) :: mflo(3),mfhi(3)
      double precision, intent(inout) :: mf(mflo(1):mfhi(1),mflo(2):mfhi(2),mflo(3):mfhi(3),ncomp)
      integer         , intent(in   ) :: ncomp,comp

      ! local
      integer :: i,j,k
      
      !=============================
      ! fill elements of array with random numbers
      !=============================
      do k=lo(3),hi(3)
      do j=lo(2),hi(2)
      do i=lo(1),hi(1)
         mf(i,j,k,comp+1) = get_fhd_normal_func()
      end do
      end do
      end do

    end subroutine multifab_fill_random

#endif

    subroutine multifab_fill_random_hack(lo,hi, &
                                    mf,mflo,mfhi, &
                                    ncomp,comp, mytop, mybottom) bind(C,name="multifab_fill_random_hack")

      integer         , intent(in   ) :: lo(3),hi(3)
      double precision, intent(in   ) :: mytop,mybottom
      integer         , intent(in   ) :: mflo(3),mfhi(3)
      double precision, intent(inout) :: mf(mflo(1):mfhi(1),mflo(2):mfhi(2),mflo(3):mfhi(3),ncomp)
      integer         , intent(in   ) :: ncomp,comp

      ! local
      integer :: i,j,k
      
      !=============================
      ! fill elements of array with random numbers
      !=============================
      do k=lo(3),hi(3)
      do j=lo(2),hi(2)
      do i=lo(1)+1,hi(1)-1
         mf(i,j,k,comp+1) = get_fhd_normal_func()
      end do
      end do
      end do

      do k=lo(3),hi(3)
      do j=lo(2),hi(2)

         mf(lo(1),j,k,comp+1) = mybottom
         mf(hi(1),j,k,comp+1) = mytop

      end do
      end do

    end subroutine multifab_fill_random_hack
   

end module multifab_fill_random_module
