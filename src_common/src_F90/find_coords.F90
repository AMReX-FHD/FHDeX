#include <AMReX_Config.H>

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


