module convert_stag_module

  use amrex_error_module

  implicit none

  private

contains

  subroutine cc_restriction(lo_c,hi_c, &
                            phi_c, c_lo, c_hi, &
                            phi_f, f_lo, f_hi) &
                            bind (C,name="cc_restriction")

    integer         , intent(in   ) :: lo_c(3), hi_c(3)
    integer         , intent(in   ) :: c_lo(3), c_hi(3)
    integer         , intent(in   ) :: f_lo(3), f_hi(3)
    double precision, intent(inout) :: phi_c(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3))
    double precision, intent(in   ) :: phi_f(c_lo(1):c_hi(1),c_lo(2):c_hi(2),c_lo(3):c_hi(3))
    ! local
    integer :: i,j,k

#if (AMREX_SPACEDIM == 2)
    
    k = 0
    do j=lo_c(2),hi_c(2)
    do i=lo_c(1),hi_c(1)
       phi_c(i,j,k) = 0.25d0*(  phi_f(2*i,2*j  ,k) + phi_f(2*i+1,2*j  ,k) &
                              + phi_f(2*i,2*j+1,k) + phi_f(2*i+1,2*j+1,k) )
    end do
    end do

#elif (AMREX_SPACEDIM == 3)

    do k=lo_c(3),hi_c(3)
    do j=lo_c(2),hi_c(2)
    do i=lo_c(1),hi_c(1)
       phi_c(i,j,k) = 0.125d0*(  phi_f(2*i,2*j  ,2*k  ) + phi_f(2*i+1,2*j  ,2*k  ) &
                               + phi_f(2*i,2*j+1,2*k  ) + phi_f(2*i+1,2*j+1,2*k  ) &
                               + phi_f(2*i,2*j  ,2*k+1) + phi_f(2*i+1,2*j  ,2*k+1) &
                               + phi_f(2*i,2*j+1,2*k+1) + phi_f(2*i+1,2*j+1,2*k+1) )
    end do
    end do
    end do

#endif

  end subroutine cc_restriction

end module convert_stag_module
