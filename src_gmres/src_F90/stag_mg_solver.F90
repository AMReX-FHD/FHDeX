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

  subroutine face_restriction(lo_c,hi_c, &
                              phix_c, cx_lo, cx_hi, &
                              phix_f, fx_lo, fx_hi, &
                              phiy_c, cy_lo, cy_hi, &
                              phiy_f, fy_lo, fy_hi, &
#if (AMREX_SPACEDIM == 3)
                              phiz_c, cz_lo, cz_hi, &
                              phiz_f, fz_lo, fz_hi, &
#endif
                              simple_stencil) &
                              bind (C,name="face_restriction")

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

      if (simple_stencil .eq. 1) then

         ! 2 point stencils
         do j=lo_c(2),hi_c(2)
         do i=lo_c(1),hi_c(1)+1
            phix_c(i,j,k) = 0.5d0*(phix_f(2*i,2*j,k) + phix_f(2*i,2*j+1,k))
         end do
         end do

         do j=lo_c(2),hi_c(2)+1
         do i=lo_c(1),hi_c(1)
            phiy_c(i,j,k) = 0.5d0*(phiy_f(2*i,2*j,k) + phiy_f(2*i+1,2*j,k))
         end do
         end do

      else

         ! 6 point stencils
         do j=lo_c(2),hi_c(2)
         do i=lo_c(1),hi_c(1)+1
            phix_c(i,j,k) = 0.25d0*(phix_f(2*i,2*j,k) + phix_f(2*i,2*j+1,k)) &
                 + 0.125d0*( phix_f(2*i+1,2*j,k) + phix_f(2*i+1,2*j+1,k) &
                 +phix_f(2*i-1,2*j,k) + phix_f(2*i-1,2*j+1,k))
         end do
         end do

         do j=lo_c(2),hi_c(2)+1
         do i=lo_c(1),hi_c(1)
            phiy_c(i,j,k) = 0.25d0*(phiy_f(2*i,2*j,k) + phiy_f(2*i+1,2*j,k)) &
                 + 0.125d0*( phiy_f(2*i,2*j+1,k) + phiy_f(2*i+1,2*j+1,k) &
                 +phiy_f(2*i,2*j-1,k) + phiy_f(2*i+1,2*j-1,k))
         end do
         end do

      end if

#elif (AMREX_SPACEDIM == 3)

      if (simple_stencil .eq. 1) then

         ! 4 point stencils
         do k=lo_c(3),hi_c(3)
         do j=lo_c(2),hi_c(2)
         do i=lo_c(1),hi_c(1)+1
            phix_c(i,j,k) = 0.25d0* ( phix_f(2*i,2*j,2*k  ) + phix_f(2*i,2*j+1,2*k  ) &
                 +phix_f(2*i,2*j,2*k+1) + phix_f(2*i,2*j+1,2*k+1) )
         end do
         end do
         end do

         do k=lo_c(3),hi_c(3)
         do j=lo_c(2),hi_c(2)+1
         do i=lo_c(1),hi_c(1)
            phiy_c(i,j,k) = 0.25d0* ( phiy_f(2*i,2*j,2*k  ) + phiy_f(2*i+1,2*j,2*k  ) &
                 +phiy_f(2*i,2*j,2*k+1) + phiy_f(2*i+1,2*j,2*k+1) )
         end do
         end do
         end do

         do k=lo_c(3),hi_c(3)+1
         do j=lo_c(2),hi_c(2)
         do i=lo_c(1),hi_c(1)
            phiz_c(i,j,k) = 0.25d0* ( phiz_f(2*i,2*j  ,2*k) + phiz_f(2*i+1,2*j  ,2*k) &
                 +phiz_f(2*i,2*j+1,2*k) + phiz_f(2*i+1,2*j+1,2*k) )
         end do
         end do
         end do

      else

         ! 12 point stencils
         do k=lo_c(3),hi_c(3)
         do j=lo_c(2),hi_c(2)
         do i=lo_c(1),hi_c(1)+1
            phix_c(i,j,k) = 0.125d0* ( phix_f(2*i,2*j,2*k  ) + phix_f(2*i,2*j+1,2*k  ) &
                 +phix_f(2*i,2*j,2*k+1) + phix_f(2*i,2*j+1,2*k+1) ) &
                 + 0.0625* ( phix_f(2*i+1,2*j,2*k  ) + phix_f(2*i+1,2*j+1,2*k  ) &
                 +phix_f(2*i+1,2*j,2*k+1) + phix_f(2*i+1,2*j+1,2*k+1) ) &
                 + 0.0625* ( phix_f(2*i-1,2*j,2*k  ) + phix_f(2*i-1,2*j+1,2*k  ) &
                 +phix_f(2*i-1,2*j,2*k+1) + phix_f(2*i-1,2*j+1,2*k+1) )
         end do
         end do
         end do

         do k=lo_c(3),hi_c(3)
         do j=lo_c(2),hi_c(2)+1
         do i=lo_c(1),hi_c(1)
            phiy_c(i,j,k) = 0.125d0* ( phiy_f(2*i,2*j,2*k  ) + phiy_f(2*i+1,2*j,2*k  ) &
                 +phiy_f(2*i,2*j,2*k+1) + phiy_f(2*i+1,2*j,2*k+1) ) &
                 + 0.0625* ( phiy_f(2*i,2*j+1,2*k  ) + phiy_f(2*i+1,2*j+1,2*k  ) &
                 +phiy_f(2*i,2*j+1,2*k+1) + phiy_f(2*i+1,2*j+1,2*k+1) ) &
                 + 0.0625* ( phiy_f(2*i,2*j-1,2*k  ) + phiy_f(2*i+1,2*j-1,2*k  ) &
                 +phiy_f(2*i,2*j-1,2*k+1) + phiy_f(2*i+1,2*j-1,2*k+1) )
         end do
         end do
         end do

         do k=lo_c(3),hi_c(3)+1
         do j=lo_c(2),hi_c(2)
         do i=lo_c(1),hi_c(1)
            phiz_c(i,j,k) = 0.125d0* ( phiz_f(2*i,2*j  ,2*k) + phiz_f(2*i+1,2*j  ,2*k) &
                 +phiz_f(2*i,2*j+1,2*k) + phiz_f(2*i+1,2*j+1,2*k) ) &
                 + 0.0625d0* ( phiz_f(2*i,2*j  ,2*k+1) + phiz_f(2*i+1,2*j  ,2*k+1) &
                 +phiz_f(2*i,2*j+1,2*k+1) + phiz_f(2*i+1,2*j+1,2*k+1) ) &
                 + 0.0625d0* ( phiz_f(2*i,2*j  ,2*k-1) + phiz_f(2*i+1,2*j  ,2*k-1) &
                 +phiz_f(2*i,2*j+1,2*k-1) + phiz_f(2*i+1,2*j+1,2*k-1) )
         end do
         end do
         end do

      end if

#endif

  end subroutine face_restriction

  subroutine edge_restriction(lo_c,hi_c, &
                              phixy_c, cxy_lo, cxy_hi, &
                              phixy_f, fxy_lo, fxy_hi, &
                              phixz_c, cxz_lo, cxz_hi, &
                              phixz_f, fxz_lo, fxz_hi, &
                              phiyz_c, cyz_lo, cyz_hi, &
                              phiyz_f, fyz_lo, fyz_hi) &
                              bind (C,name="edge_restriction")

    integer         , intent(in   ) :: lo_c(3), hi_c(3)
    integer         , intent(in   ) :: cxy_lo(3), cxy_hi(3)
    integer         , intent(in   ) :: fxy_lo(3), fxy_hi(3)
    double precision, intent(inout) :: phixy_c(fxy_lo(1):fxy_hi(1),fxy_lo(2):fxy_hi(2),fxy_lo(3):fxy_hi(3))
    double precision, intent(in   ) :: phixy_f(cxy_lo(1):cxy_hi(1),cxy_lo(2):cxy_hi(2),cxy_lo(3):cxy_hi(3))
    integer         , intent(in   ) :: cxz_lo(3), cxz_hi(3)
    integer         , intent(in   ) :: fxz_lo(3), fxz_hi(3)
    double precision, intent(inout) :: phixz_c(fxz_lo(1):fxz_hi(1),fxz_lo(2):fxz_hi(2),fxz_lo(3):fxz_hi(3))
    double precision, intent(in   ) :: phixz_f(cxz_lo(1):cxz_hi(1),cxz_lo(2):cxz_hi(2),cxz_lo(3):cxz_hi(3))
    integer         , intent(in   ) :: cyz_lo(3), cyz_hi(3)
    integer         , intent(in   ) :: fyz_lo(3), fyz_hi(3)
    double precision, intent(inout) :: phiyz_c(fyz_lo(1):fyz_hi(1),fyz_lo(2):fyz_hi(2),fyz_lo(3):fyz_hi(3))
    double precision, intent(in   ) :: phiyz_f(cyz_lo(1):cyz_hi(1),cyz_lo(2):cyz_hi(2),cyz_lo(3):cyz_hi(3))

    ! local
    integer :: i,j,k

      ! xy edges
      do k=lo_c(3),hi_c(3)
      do j=lo_c(2),hi_c(2)+1
      do i=lo_c(1),hi_c(1)+1
         phixy_c(i,j,k) = 0.5d0*(phixy_f(2*i,2*j,2*k)+phixy_f(2*i,2*j,2*k+1))
      end do
      end do
      end do

      ! xz edges
      do k=lo_c(3),hi_c(3)+1
      do j=lo_c(2),hi_c(2)
      do i=lo_c(1),hi_c(1)+1
         phixz_c(i,j,k) =  0.5d0*(phixz_f(2*i,2*j,2*k)+phixz_f(2*i,2*j+1,2*k))
      end do
      end do
      end do

      ! yz edges
      do k=lo_c(3),hi_c(3)+1
      do j=lo_c(2),hi_c(2)+1
      do i=lo_c(1),hi_c(1)
         phiyz_c(i,j,k) =  0.5d0*(phiyz_f(2*i,2*j,2*k)+phiyz_f(2*i+1,2*j,2*k))
      end do
      end do
      end do

  end subroutine edge_restriction

  subroutine nodal_restriction(lo_c,hi_c, &
                            phi_c, c_lo, c_hi, &
                            phi_f, f_lo, f_hi) &
                            bind (C,name="nodal_restriction")

    integer         , intent(in   ) :: lo_c(3), hi_c(3)
    integer         , intent(in   ) :: c_lo(3), c_hi(3)
    integer         , intent(in   ) :: f_lo(3), f_hi(3)
    double precision, intent(inout) :: phi_c(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3))
    double precision, intent(in   ) :: phi_f(c_lo(1):c_hi(1),c_lo(2):c_hi(2),c_lo(3):c_hi(3))
    ! local
    integer :: i,j,k

#if (AMREX_SPACEDIM == 2)
    
    k = 0
    do j=lo_c(2),hi_c(2)+1
    do i=lo_c(1),hi_c(1)+1
       phi_c(i,j,k) = phi_f(2*i,2*j,k)
    end do
    end do

#elif (AMREX_SPACEDIM == 3)

      do k=lo_c(3),hi_c(3)+1
      do j=lo_c(2),hi_c(2)+1
      do i=lo_c(1),hi_c(1)+1
         phi_c(i,j,k) = phi_f(2*i,2*j,2*k)
      end do
      end do
      end do

#endif

  end subroutine nodal_restriction

end module convert_stag_module
