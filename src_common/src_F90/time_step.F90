
module time_step_module

  use amrex_error_module

  implicit none

  private

contains

#if AMREX_SPACEDIM == 2
  subroutine euler_step_stag(lo, hi, old, oldlo, oldhi, stagop, stagoplo, stagophi, newdata, newdatalo, newdatahi, offset, dt) bind(C, name="euler_step_stag")

    use amrex_fort_module, only : amrex_real
    implicit none

	  integer, intent(in) :: lo(2), hi(2), oldlo(2), oldhi(2), newdatalo(2), newdatahi(2), stagoplo(2), stagophi(2), offset(2)
	  real(amrex_real), intent(in   ) :: dt
    real(amrex_real), intent(in) :: stagop(stagoplo(1):stagophi(1),stagoplo(2):stagophi(2))
    real(amrex_real), intent(inout) :: old(oldlo(1):oldhi(1),oldlo(2):oldhi(2))
    real(amrex_real), intent(inout) :: newdata(newdatalo(1):newdatahi(1),newdatalo(2):newdatahi(2))

    ! local variables
    integer i,j

    do j = lo(2), hi(2)+offset(2)
    	do i = lo(1), hi(1)+offset(1)
        newdata(i,j) = old(i,j)-stagop(i,j)*dt
        !old(i,j) = newdata(i,j)
    	end do
    end do

  end subroutine euler_step_stag
#endif

#if AMREX_SPACEDIM == 3
  subroutine euler_step_stag(lo, hi, old, oldlo, oldhi, stagop, stagoplo, stagophi, newdata, newdatalo, newdatahi, offset, dt) bind(C, name="euler_step_stag")

    use amrex_fort_module, only : amrex_real
    implicit none

	  integer, intent(in) :: lo(3), hi(3), oldlo(3), oldhi(3), newdatalo(3), newdatahi(3), stagoplo(3), stagophi(3), offset(3)
	  real(amrex_real), intent(in   ) :: dt
    real(amrex_real), intent(in) :: stagop(stagoplo(1):stagophi(1),stagoplo(2):stagophi(2),stagoplo(3):stagophi(3))
    real(amrex_real), intent(inout) :: old(oldlo(1):oldhi(1),oldlo(2):oldhi(2),oldlo(3):oldhi(3))
    real(amrex_real), intent(inout) :: newdata(newdatalo(1):newdatahi(1),newdatalo(2):newdatahi(2),newdatalo(3):newdatahi(3))

    ! local variables
    integer i,j,k

    do k = lo(3), hi(3)+offset(3)
      do j = lo(2), hi(2)+offset(2)
      	do i = lo(1), hi(1)+offset(1)
          !minus sign due to operator construction
          newdata(i,j,k) = old(i,j,k)-stagop(i,j,k)*dt
          !newdata(i,j,k) = old(i,j,k)
      	end do
      end do
    end do
  end subroutine euler_step_stag
#endif


end module time_step_module

