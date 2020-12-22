#include <AMReX_Config.H>


subroutine find_face_coords(real_lo, real_hi, xface, xfacelo, xfacehi, yface, yfacelo, yfacehi, &
#if (BL_SPACEDIM == 3)
                             zface, zfacelo, zfacehi, &
#endif
                             dx) bind(C, name="find_face_coords")

  use amrex_fort_module, only : amrex_real
  implicit none

  integer, intent(in) :: xfacelo(3), xfacehi(3), yfacelo(3), yfacehi(3)
#if (BL_SPACEDIM == 3)
  integer, intent(in) :: zfacelo(3), zfacehi(3)
#endif
  real(amrex_real), intent(in   ) :: real_lo(3), real_hi(3), dx(3)

  real(amrex_real), intent(inout) :: xface(xfacelo(1):xfacehi(1),xfacelo(2):xfacehi(2),xfacelo(3):xfacehi(3),1:AMREX_SPACEDIM)
  real(amrex_real), intent(inout) :: yface(yfacelo(1):yfacehi(1),yfacelo(2):yfacehi(2),yfacelo(3):yfacehi(3),1:AMREX_SPACEDIM)
#if (BL_SPACEDIM == 3)
  real(amrex_real), intent(inout) :: zface(zfacelo(1):zfacehi(1),zfacelo(2):zfacehi(2),zfacelo(3):zfacehi(3),1:AMREX_SPACEDIM)
#endif

  ! local variables
  integer i,j,k

  !x face coordinates
  do k = xfacelo(3), xfacehi(3)
    do j = xfacelo(2), xfacehi(2)
      do i = xfacelo(1), xfacehi(1) !iterate into ghost cells
         xface(i,j,k,1) = i*dx(1)+real_lo(1)
      end do
    end do
  end do

  do k = xfacelo(3), xfacehi(3)
    do j = xfacelo(2), xfacehi(2)
      do i = xfacelo(1), xfacehi(1)
         xface(i,j,k,2) = (j+0.5)*dx(2)+real_lo(2)
      end do
    end do
  end do

#if (BL_SPACEDIM == 3)
  do k = xfacelo(3), xfacehi(3)
    do j = xfacelo(2), xfacehi(2)
      do i = xfacelo(1), xfacehi(1)
         xface(i,j,k,3) = (k+0.5)*dx(3)+real_lo(3)
      end do
    end do
  end do
#endif

  do k = yfacelo(3), yfacehi(3)
    do j = yfacelo(2), yfacehi(2)
      do i = yfacelo(1), yfacehi(1)
         yface(i,j,k,1) = (i+0.5)*dx(1)+real_lo(1)
      end do
    end do
  end do

  do k = yfacelo(3), yfacehi(3)
    do j = yfacelo(2), yfacehi(2)
      do i = yfacelo(1), yfacehi(1)
         yface(i,j,k,2) = (j)*dx(2)+real_lo(2)
      end do
    end do
  end do

#if (BL_SPACEDIM == 3)
  do k = yfacelo(3), yfacehi(3)
    do j = yfacelo(2), yfacehi(2)
      do i = yfacelo(1), yfacehi(1)
         yface(i,j,k,3) = (k+0.5)*dx(3)+real_lo(3)
      end do
    end do
  end do
#endif

#if (BL_SPACEDIM == 3)
  do k = zfacelo(3), zfacehi(3)
    do j = zfacelo(2), zfacehi(2)
      do i = zfacelo(1), zfacehi(1) !iterate into ghost cells
         zface(i,j,k,1) = (i+0.5)*dx(1)+real_lo(1)
      end do
    end do
  end do

  do k = zfacelo(3), zfacehi(3)
    do j = zfacelo(2), zfacehi(2)
      do i = zfacelo(1), zfacehi(1)
         zface(i,j,k,2) = (j+0.5)*dx(2)+real_lo(2)
      end do
    end do
  end do

  do k = zfacelo(3), zfacehi(3)
    do j = zfacelo(2), zfacehi(2)
      do i = zfacelo(1), zfacehi(1)
         zface(i,j,k,3) = (k)*dx(3)+real_lo(3)
      end do
    end do
  end do
#endif

end subroutine find_face_coords

subroutine find_center_coords(real_lo, real_hi, centers, lo, hi, dx) bind(C, name="find_center_coords")

  use amrex_fort_module, only : amrex_real
  implicit none

  integer,          intent(in   ) :: lo(3), hi(3)
  real(amrex_real), intent(in   ) :: real_lo(3), real_hi(3), dx(3)

#if (BL_SPACEDIM == 3)
  real(amrex_real), intent(inout) :: centers(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1:AMREX_SPACEDIM)
#endif

#if (BL_SPACEDIM == 2)
  real(amrex_real), intent(inout) :: centers(lo(1):hi(1),lo(2):hi(2),1:AMREX_SPACEDIM)
#endif

#if (BL_SPACEDIM == 1)
  real(amrex_real), intent(inout) :: centers(lo(1):hi(1),1:AMREX_SPACEDIM)
#endif

#if (BL_SPACEDIM == 3)
  integer i,j,k

  do k = lo(3), hi(3)
    do j = lo(2), hi(2)
      do i = lo(1), hi(1) !iterate into ghost cells
         centers(i,j,k,1) = (i+0.5)*dx(1)+real_lo(1)
         centers(i,j,k,2) = (j+0.5)*dx(2)+real_lo(2)
         centers(i,j,k,3) = (k+0.5)*dx(3)+real_lo(3)
      end do
    end do
  end do
#endif

#if (BL_SPACEDIM == 2)
  integer i,j

    do j = lo(2), hi(2)
      do i = lo(1), hi(1) !iterate into ghost cells
         centers(i,j,1) = (i+0.5)*dx(1)+real_lo(1)
         centers(i,j,2) = (j+0.5)*dx(2)+real_lo(2)
      end do
    end do
#endif


#if (BL_SPACEDIM == 1)
  integer i

      do i = lo(1), hi(1) !iterate into ghost cells
         centers(i,1) = (i+0.5)*dx(1)+real_lo(1)
    end do
#endif

end subroutine find_center_coords


