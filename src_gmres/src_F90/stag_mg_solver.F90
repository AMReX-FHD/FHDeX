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

  subroutine stag_restriction(lo_c,hi_c, &
                              phix_c, cx_lo, cx_hi, &
                              phix_f, fx_lo, fx_hi, &
                              phiy_c, cy_lo, cy_hi, &
                              phiy_f, fy_lo, fy_hi, &
#if (AMREX_SPACEDIM == 3)
                              phiz_c, cz_lo, cz_hi, &
                              phiz_f, fz_lo, fz_hi, &
#endif
                              simple_stencil) &
                              bind (C,name="stag_restriction")

    integer         , intent(in   ) :: lo_c(3), hi_c(3)
    integer         , intent(in   ) :: cx_lo(3), cx_hi(3)
    integer         , intent(in   ) :: fx_lo(3), fx_hi(3)
    double precision, intent(inout) :: phix_c(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2),fx_lo(3):fx_hi(3))
    double precision, intent(in   ) :: phix_f(cx_lo(1):cx_hi(1),cx_lo(2):cx_hi(2),cx_lo(3):cx_hi(3))
    integer         , intent(in   ) :: cy_lo(3), cy_hi(3)
    integer         , intent(in   ) :: fy_lo(3), fy_hi(3)
    double precision, intent(inout) :: phiy_c(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2),fy_lo(3):fy_hi(3))
    double precision, intent(in   ) :: phiy_f(cy_lo(1):cy_hi(1),cy_lo(2):cy_hi(2),cy_lo(3):cy_hi(3))
#if (AMREX_SPACEDIM == 3)
    integer         , intent(in   ) :: cz_lo(3), cz_hi(3)
    integer         , intent(in   ) :: fz_lo(3), fz_hi(3)
    double precision, intent(inout) :: phiz_c(fz_lo(1):fz_hi(1),fz_lo(2):fz_hi(2),fz_lo(3):fz_hi(3))
    double precision, intent(in   ) :: phiz_f(cz_lo(1):cz_hi(1),cz_lo(2):cz_hi(2),cz_lo(3):cz_hi(3))
#endif
    integer         , intent(in   ) :: simple_stencil

    ! local
    integer :: i,j,k

#if (AMREX_SPACEDIM == 2)


#elif (AMREX_SPACEDIM == 3)


#endif

  end subroutine stag_restriction

end module convert_stag_module
