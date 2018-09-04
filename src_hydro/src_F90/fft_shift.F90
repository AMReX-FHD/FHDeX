module fft_shift_module

  use amrex_error_module
  ! use common_namelist_module, only: visc_type

  implicit none

  private

contains

#if (AMREX_SPACEDIM == 2)

    subroutine fft_shift(lo,hi, &
                         dft, dftlo, dfthi) bind(C,name="fft_shift")

      integer         , intent(in   ) :: lo(2),hi(2)
      integer         , intent(in   ) :: dftlo(2),dfthi(2)
      double precision, intent(inout) :: dft(dftlo(1):dfthi(1),dftlo(2):dfthi(2))

      ! local
      integer :: i,j
      integer :: ip,jp
      integer :: nx,ny, nxh,nyh
      double precision :: dft_temp(dftlo(1):dfthi(1),dftlo(2):dfthi(2))
      
      nx = hi(1) - lo(1) + 1
      ny = hi(2) - lo(2) + 1

      ! print*, "nx, ny", nx, ny
      
      if ((mod(nx,2)==1).OR.(mod(ny,2)==1)) then
         print*, "ERROR: Odd dimensions, fftshift will not work"
         stop
      end if

      nxh = nx/2
      nyh = ny/2

      dft_temp(dftlo(1):dfthi(1),dftlo(2):dfthi(2)) = &
           dft(dftlo(1):dfthi(1),dftlo(2):dfthi(2))

      do j = lo(2),hi(2)
      do i = lo(1),hi(1)
         ! Find shifted indices
         if (i.gt.nxh) then
            ip = i - (nxh + 1)
         else 
            ip = i + (nxh - 1)
         end if

         if (j.gt.nxh) then
            jp = j - (nxh + 1)
         else 
            jp = j + (nxh - 1)
         end if
         
         ! Switch values
         dft(ip,jp) = dft_temp(i,j)
      end do
      end do

    end subroutine fft_shift

#endif

#if (AMREX_SPACEDIM == 3)

    subroutine fft_shift(lo,hi, &
                         dft, dftlo, dfthi) bind(C,name="fft_shift")

      integer         , intent(in   ) :: lo(3),hi(3)
      integer         , intent(in   ) :: dftlo(3),dfthi(3)
      double precision, intent(inout) :: dft(dftlo(1):dfthi(1),dftlo(2):dfthi(2),dftlo(3):dfthi(3))

      ! local
      integer :: i,j,k
      integer :: ip,jp,kp
      integer :: nx,ny,nz, nxh,nyh,nzh
      double precision :: dft_temp(dftlo(1):dfthi(1),dftlo(2):dfthi(2),dftlo(3):dfthi(3))
      
      nx = hi(1) - lo(1) + 1
      ny = hi(2) - lo(2) + 1
      nz = hi(3) - lo(3) + 1
      
      if ((mod(nx,2)==1).OR.(mod(ny,2)==1).OR.(mod(nz,2)==1)) then
         print*, "ERROR: Odd dimensions, fftshift will not work"
         stop
      end if

      nxh = nx/2
      nyh = ny/2
      nzh = nz/2

      dft_temp(dftlo(1):dfthi(1),dftlo(2):dfthi(2),dftlo(3):dfthi(3)) = &
           dft(dftlo(1):dfthi(1),dftlo(2):dfthi(2),dftlo(3):dfthi(3))

      do k = lo(3),hi(3)
      do j = lo(2),hi(2)
      do i = lo(1),hi(1)
         ! Find shifted indices
         if (i.gt.nxh) then
            ip = i - (nxh + 1)
         else 
            ip = i + (nxh - 1)
         end if

         if (j.gt.nxh) then
            jp = j - (nxh + 1)
         else 
            jp = j + (nxh - 1)
         end if

         if (k.gt.nxh) then
            kp = k - (nxh + 1)
         else 
            kp = k + (nxh - 1)
         end if
         
         ! Switch values
         dft(ip,jp,kp) = dft_temp(i,j,k)
      end do
      end do
      end do
      
    end subroutine fft_shift
    
#endif

end module fft_shift_module
