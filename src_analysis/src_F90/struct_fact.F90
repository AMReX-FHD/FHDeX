module struct_fact_module

  use amrex_error_module

  implicit none

  private

contains

#if (AMREX_SPACEDIM == 2)

    subroutine fft_shift(lo,hi, &
                         dft, dftlo, dfthi, ncomp, &
                         zero_avg) bind(C,name="fft_shift")

      integer         , intent(in   ) :: ncomp, zero_avg
      integer         , intent(in   ) :: lo(2),hi(2)
      integer         , intent(in   ) :: dftlo(2),dfthi(2)
      double precision, intent(inout) :: dft(dftlo(1):dfthi(1),dftlo(2):dfthi(2),ncomp)

      ! local
      integer :: i,j,n
      integer :: ip,jp
      integer :: nx,ny, nxh,nyh
      double precision :: dft_temp(dftlo(1):dfthi(1),dftlo(2):dfthi(2))
      
      nx = hi(1) - lo(1) + 1
      ny = hi(2) - lo(2) + 1

      ! print*, "nx, ny", nx, ny

      ! if ((mod(nx,2)==1).OR.(mod(ny,2)==1)) then
      !    print*, "ERROR: Odd dimensions, fftshift will not work"
      !    stop
      ! end if

      nxh = (nx+1)/2
      nyh = (ny+1)/2

      do n = 1,ncomp
         
         if (zero_avg.eq.1) then
            ! set k=0 cell to zero
            dft(lo(1),lo(2),n) = 0.0d0
         endif

         dft_temp(dftlo(1):dfthi(1),dftlo(2):dfthi(2)) = &
              dft(dftlo(1):dfthi(1),dftlo(2):dfthi(2),n)

         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               ! Find shifted indices
               ip = MOD(i+nxh,nx)
               jp = MOD(j+nyh,ny)

               ! Switch values
               dft(ip,jp,n) = dft_temp(i,j)
            end do
         end do

      end do

    end subroutine fft_shift

#endif

#if (AMREX_SPACEDIM == 3)

    subroutine fft_shift(lo,hi, &
                         dft, dftlo, dfthi, ncomp, &
                         zero_avg) bind(C,name="fft_shift")
      
      integer         , intent(in   ) :: ncomp, zero_avg
      integer         , intent(in   ) :: lo(3),hi(3)
      integer         , intent(in   ) :: dftlo(3),dfthi(3)
      double precision, intent(inout) :: dft(dftlo(1):dfthi(1),dftlo(2):dfthi(2),dftlo(3):dfthi(3),ncomp)

      ! local
      integer :: i,j,k,n
      integer :: ip,jp,kp
      integer :: nx,ny,nz, nxh,nyh,nzh
      double precision :: dft_temp(dftlo(1):dfthi(1),dftlo(2):dfthi(2),dftlo(3):dfthi(3))
      
      nx = hi(1) - lo(1) + 1
      ny = hi(2) - lo(2) + 1
      nz = hi(3) - lo(3) + 1
      
      ! if ((mod(nx,2)==1).OR.(mod(ny,2)==1).OR.(mod(nz,2)==1)) then
      !    print*, "ERROR: Odd dimensions, fftshift will not work"
      !    stop
      ! end if
      
      ! Take ceiling
      nxh = (nx+1)/2
      nyh = (ny+1)/2
      nzh = (nz+1)/2

      do n = 1,ncomp
         
         if (zero_avg.eq.1) then
            ! set k=0 cell to zero
            dft(lo(1),lo(2),lo(3),n) = 0.0d0
         endif

         dft_temp(dftlo(1):dfthi(1),dftlo(2):dfthi(2),dftlo(3):dfthi(3)) = &
              dft(dftlo(1):dfthi(1),dftlo(2):dfthi(2),dftlo(3):dfthi(3),n)

         do k = lo(3),hi(3)
            do j = lo(2),hi(2)
               do i = lo(1),hi(1)
                  ! Find shifted indices
                  ip = MOD(i+nxh,nx)
                  jp = MOD(j+nyh,ny)
                  kp = MOD(k+nzh,nz)

                  ! Switch values
                  dft(ip,jp,kp,n) = dft_temp(i,j,k)
               end do
            end do
         end do

      end do
      
    end subroutine fft_shift
    
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

end module struct_fact_module
