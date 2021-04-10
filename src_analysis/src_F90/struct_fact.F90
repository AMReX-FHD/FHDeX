module struct_fact_module

  use amrex_error_module

  implicit none

  private

contains

#if (AMREX_SPACEDIM == 2)

    subroutine fft_shift(lo,hi, &
                         dft, dftlo, dfthi, &
                         dft_temp, dft_templo, dft_temphi, &
                         zero_avg) bind(C,name="fft_shift")

      integer         , intent(in   ) :: zero_avg
      integer         , intent(in   ) :: lo(2),hi(2)
      integer         , intent(in   ) :: dftlo(2),dfthi(2)
      double precision, intent(inout) :: dft(dftlo(1):dfthi(1),dftlo(2):dfthi(2))
      integer         , intent(in   ) :: dft_templo(2),dft_temphi(2)
      double precision, intent(inout) :: dft_temp(dft_templo(1):dft_temphi(1),dft_templo(2):dft_temphi(2))

      ! local
      integer :: i,j,n
      integer :: ip,jp
      integer :: nx,ny, nxh,nyh
      
      nx = hi(1) - lo(1) + 1
      ny = hi(2) - lo(2) + 1

      ! print*, "nx, ny", nx, ny

      ! if ((mod(nx,2)==1).OR.(mod(ny,2)==1)) then
      !    print*, "ERROR: Odd dimensions, fftshift will not work"
      !    stop
      ! end if

      nxh = (nx+1)/2
      nyh = (ny+1)/2

      if (zero_avg.eq.1) then
         ! set k=0 cell to zero
         dft(lo(1),lo(2)) = 0.0d0
      endif

      dft_temp(dftlo(1):dfthi(1),dftlo(2):dfthi(2)) = &
           dft(dftlo(1):dfthi(1),dftlo(2):dfthi(2))

      do j = lo(2),hi(2)
         do i = lo(1),hi(1)
            ! Find shifted indices
            ip = MOD(i+nxh,nx)
            jp = MOD(j+nyh,ny)

            ! Switch values
            dft(ip,jp) = dft_temp(i,j)
         end do
      end do

    end subroutine fft_shift

#endif

#if (AMREX_SPACEDIM == 3)

    subroutine fft_shift(lo,hi, &
                         dft, dftlo, dfthi, &
                         dft_temp, dft_templo, dft_temphi, &
                         zero_avg) bind(C,name="fft_shift")
      
      integer         , intent(in   ) :: zero_avg
      integer         , intent(in   ) :: lo(3),hi(3)
      integer         , intent(in   ) :: dftlo(3),dfthi(3)
      double precision, intent(inout) :: dft(dftlo(1):dfthi(1),dftlo(2):dfthi(2),dftlo(3):dfthi(3))
      integer         , intent(in   ) :: dft_templo(3),dft_temphi(3)
      double precision, intent(inout) :: dft_temp(dft_templo(1):dft_temphi(1),dft_templo(2):dft_temphi(2),dft_templo(3):dft_temphi(3))

      ! local
      integer :: i,j,k,n
      integer :: ip,jp,kp
      integer :: nx,ny,nz, nxh,nyh,nzh

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

      if (zero_avg.eq.1) then
         ! set k=0 cell to zero
         dft(lo(1),lo(2),lo(3)) = 0.0d0
      endif

      dft_temp(dftlo(1):dfthi(1),dftlo(2):dfthi(2),dftlo(3):dfthi(3)) = &
           dft(dftlo(1):dfthi(1),dftlo(2):dfthi(2),dftlo(3):dfthi(3))

      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               ! Find shifted indices
               ip = MOD(i+nxh,nx)
               jp = MOD(j+nyh,ny)
               kp = MOD(k+nzh,nz)

               ! Switch values
               dft(ip,jp,kp) = dft_temp(i,j,k)
            end do
         end do
      end do

    end subroutine fft_shift
    
#endif

end module struct_fact_module
