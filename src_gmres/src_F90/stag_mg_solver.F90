module stag_solver_module

  use amrex_error_module
  use common_namelist_module
  use gmres_namelist_module

  implicit none

  private

contains

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
    double precision, intent(inout) :: phix_c(cx_lo(1):cx_hi(1),cx_lo(2):cx_hi(2),cx_lo(3):cx_hi(3))
    double precision, intent(in   ) :: phix_f(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2),fx_lo(3):fx_hi(3))
    integer         , intent(in   ) :: cy_lo(3), cy_hi(3)
    integer         , intent(in   ) :: fy_lo(3), fy_hi(3)
    double precision, intent(inout) :: phiy_c(cy_lo(1):cy_hi(1),cy_lo(2):cy_hi(2),cy_lo(3):cy_hi(3))
    double precision, intent(in   ) :: phiy_f(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2),fy_lo(3):fy_hi(3))
#if (AMREX_SPACEDIM == 3)
    integer         , intent(in   ) :: cz_lo(3), cz_hi(3)
    integer         , intent(in   ) :: fz_lo(3), fz_hi(3)
    double precision, intent(inout) :: phiz_c(cz_lo(1):cz_hi(1),cz_lo(2):cz_hi(2),cz_lo(3):cz_hi(3))
    double precision, intent(in   ) :: phiz_f(fz_lo(1):fz_hi(1),fz_lo(2):fz_hi(2),fz_lo(3):fz_hi(3))
#endif
    integer         , intent(in   ) :: simple_stencil

    ! local
    integer :: i,j,k

#if (AMREX_SPACEDIM == 2)

    k = 0

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
          phix_c(i,j,k) = 0.25d0*( phix_f(2*i  ,2*j,k) + phix_f(2*i  ,2*j+1,k)) &
                       + 0.125d0*( phix_f(2*i+1,2*j,k) + phix_f(2*i+1,2*j+1,k) &
                                  +phix_f(2*i-1,2*j,k) + phix_f(2*i-1,2*j+1,k))
       end do
       end do

       do j=lo_c(2),hi_c(2)+1
       do i=lo_c(1),hi_c(1)
          phiy_c(i,j,k) = 0.25d0*( phiy_f(2*i,2*j  ,k) + phiy_f(2*i+1,2*j  ,k)) &
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
          phix_c(i,j,k) = 0.125d0* ( phix_f(2*i  ,2*j,2*k  ) + phix_f(2*i  ,2*j+1,2*k  ) &
                                    +phix_f(2*i  ,2*j,2*k+1) + phix_f(2*i  ,2*j+1,2*k+1) ) &
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
          phiy_c(i,j,k) = 0.125d0* ( phiy_f(2*i,2*j  ,2*k  ) + phiy_f(2*i+1,2*j  ,2*k  ) &
                                    +phiy_f(2*i,2*j  ,2*k+1) + phiy_f(2*i+1,2*j  ,2*k+1) ) &
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
          phiz_c(i,j,k) = 0.125d0* ( phiz_f(2*i,2*j  ,2*k  ) + phiz_f(2*i+1,2*j  ,2*k  ) &
                                    +phiz_f(2*i,2*j+1,2*k  ) + phiz_f(2*i+1,2*j+1,2*k  ) ) &
                       + 0.0625d0* ( phiz_f(2*i,2*j  ,2*k+1) + phiz_f(2*i+1,2*j  ,2*k+1) &
                                    +phiz_f(2*i,2*j+1,2*k+1) + phiz_f(2*i+1,2*j+1,2*k+1) ) &
                       + 0.0625d0* ( phiz_f(2*i,2*j  ,2*k-1) + phiz_f(2*i+1,2*j  ,2*k-1) &
                                    +phiz_f(2*i,2*j+1,2*k-1) + phiz_f(2*i+1,2*j+1,2*k-1) )
       end do
       end do
       end do

    end if

#endif

  end subroutine stag_restriction

  subroutine nodal_restriction(lo_c,hi_c, &
                               phi_c, c_lo, c_hi, &
                               phi_f, f_lo, f_hi) &
                               bind (C,name="nodal_restriction")

    integer         , intent(in   ) :: lo_c(3), hi_c(3)
    integer         , intent(in   ) :: c_lo(3), c_hi(3)
    integer         , intent(in   ) :: f_lo(3), f_hi(3)
    double precision, intent(inout) :: phi_c(c_lo(1):c_hi(1),c_lo(2):c_hi(2),c_lo(3):c_hi(3))
    double precision, intent(in   ) :: phi_f(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3))

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
    double precision, intent(inout) :: phixy_c(cxy_lo(1):cxy_hi(1),cxy_lo(2):cxy_hi(2),cxy_lo(3):cxy_hi(3))
    double precision, intent(in   ) :: phixy_f(fxy_lo(1):fxy_hi(1),fxy_lo(2):fxy_hi(2),fxy_lo(3):fxy_hi(3))
    integer         , intent(in   ) :: cxz_lo(3), cxz_hi(3)
    integer         , intent(in   ) :: fxz_lo(3), fxz_hi(3)
    double precision, intent(inout) :: phixz_c(cxz_lo(1):cxz_hi(1),cxz_lo(2):cxz_hi(2),cxz_lo(3):cxz_hi(3))
    double precision, intent(in   ) :: phixz_f(fxz_lo(1):fxz_hi(1),fxz_lo(2):fxz_hi(2),fxz_lo(3):fxz_hi(3))
    integer         , intent(in   ) :: cyz_lo(3), cyz_hi(3)
    integer         , intent(in   ) :: fyz_lo(3), fyz_hi(3)
    double precision, intent(inout) :: phiyz_c(cyz_lo(1):cyz_hi(1),cyz_lo(2):cyz_hi(2),cyz_lo(3):cyz_hi(3))
    double precision, intent(in   ) :: phiyz_f(fyz_lo(1):fyz_hi(1),fyz_lo(2):fyz_hi(2),fyz_lo(3):fyz_hi(3))

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

  subroutine stag_prolongation(lo_f,hi_f, &
                               phix_c, cx_lo, cx_hi, &
                               phix_f, fx_lo, fx_hi, &
                               phiy_c, cy_lo, cy_hi, &
                               phiy_f, fy_lo, fy_hi  &
#if (AMREX_SPACEDIM == 3)
                             , phiz_c, cz_lo, cz_hi, &
                               phiz_f, fz_lo, fz_hi  &
#endif
                               ) &
                              bind (C,name="stag_prolongation")

    integer         , intent(in   ) :: lo_f(3), hi_f(3)
    integer         , intent(in   ) :: cx_lo(3), cx_hi(3)
    integer         , intent(in   ) :: fx_lo(3), fx_hi(3)
    double precision, intent(in   ) :: phix_c(cx_lo(1):cx_hi(1),cx_lo(2):cx_hi(2),cx_lo(3):cx_hi(3))
    double precision, intent(inout) :: phix_f(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2),fx_lo(3):fx_hi(3))
    integer         , intent(in   ) :: cy_lo(3), cy_hi(3)
    integer         , intent(in   ) :: fy_lo(3), fy_hi(3)
    double precision, intent(in   ) :: phiy_c(cy_lo(1):cy_hi(1),cy_lo(2):cy_hi(2),cy_lo(3):cy_hi(3))
    double precision, intent(inout) :: phiy_f(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2),fy_lo(3):fy_hi(3))
#if (AMREX_SPACEDIM == 3)
    integer         , intent(in   ) :: cz_lo(3), cz_hi(3)
    integer         , intent(in   ) :: fz_lo(3), fz_hi(3)
    double precision, intent(in   ) :: phiz_c(cz_lo(1):cz_hi(1),cz_lo(2):cz_hi(2),cz_lo(3):cz_hi(3))
    double precision, intent(inout) :: phiz_f(fz_lo(1):fz_hi(1),fz_lo(2):fz_hi(2),fz_lo(3):fz_hi(3))
#endif

    ! local
    integer :: i,j,k
    integer :: ioff,joff

#if (AMREX_SPACEDIM == 2)

    k = 0

    do j=lo_f(2),hi_f(2)
    do i=lo_f(1),hi_f(1)+1

       if (mod(j,2) .eq. 0) then
          joff = -1
       else
          joff = 1
       end if

       if (mod(i,2) .eq. 0) then

          ! linear interpolation
          phix_f(i,j,k) = phix_f(i,j,k) &
               + 0.75d0*phix_c(i/2,j/2     ,k) &
               + 0.25d0*phix_c(i/2,j/2+joff,k)

       else

          ! bilinear interpolation
          phix_f(i,j,k) = phix_f(i,j,k) &
               + 0.375d0*phix_c(i/2  ,j/2     ,k) &
               + 0.125d0*phix_c(i/2  ,j/2+joff,k) &
               + 0.375d0*phix_c(i/2+1,j/2     ,k) &
               + 0.125d0*phix_c(i/2+1,j/2+joff,k)

       end if

    end do
    end do

    do j=lo_f(2),hi_f(2)+1
    do i=lo_f(1),hi_f(1)

       if (mod(i,2) .eq. 0) then
          ioff = -1
       else
          ioff = 1
       end if

       if (mod(j,2) .eq. 0) then

          ! linear interpolation
          phiy_f(i,j,k) = phiy_f(i,j,k) &
               + 0.75d0*phiy_c(i/2     ,j/2,k) &
               + 0.25d0*phiy_c(i/2+ioff,j/2,k)

       else

          ! bilinear interpolation
          phiy_f(i,j,k) = phiy_f(i,j,k) &
               + 0.375d0*phiy_c(i/2     ,j/2  ,k) &
               + 0.125d0*phiy_c(i/2+ioff,j/2  ,k) &
               + 0.375d0*phiy_c(i/2     ,j/2+1,k) &
               + 0.125d0*phiy_c(i/2+ioff,j/2+1,k)

       end if

    end do
    end do

#elif (AMREX_SPACEDIM == 3)

    integer :: koff
    double precision, parameter :: nine16 = 9.d0/16.d0
    double precision, parameter :: three16 = 3.d0/16.d0
    double precision, parameter :: one16 = 1.d0/16.d0
    double precision, parameter :: nine32 = 9.d0/32.d0
    double precision, parameter :: three32 = 3.d0/32.d0
    double precision, parameter :: one32 = 1.d0/32.d0

    do k=lo_f(3),hi_f(3)
    do j=lo_f(2),hi_f(2)
    do i=lo_f(1),hi_f(1)+1

       if (mod(j,2) .eq. 0) then
          joff = -1
       else
          joff = 1
       end if

       if (mod(k,2) .eq. 0) then
          koff = -1
       else
          koff = 1
       end if

       if (mod(i,2) .eq. 0) then
          ! bilinear in the yz plane
          phix_f(i,j,k) = phix_f(i,j,k) &
               + nine16 *phix_c(i/2,j/2     ,k/2     ) &
               + three16*phix_c(i/2,j/2+joff,k/2     ) &
               + three16*phix_c(i/2,j/2     ,k/2+koff) &
               + one16  *phix_c(i/2,j/2+joff,k/2+koff)
       else
          ! bilinear in the yz plane, linear in x
          phix_f(i,j,k) = phix_f(i,j,k) &
               + nine32 *phix_c(i/2  ,j/2     ,k/2     ) &
               + three32*phix_c(i/2  ,j/2+joff,k/2     ) &
               + three32*phix_c(i/2  ,j/2     ,k/2+koff) &
               + one32  *phix_c(i/2  ,j/2+joff,k/2+koff) &
               + nine32 *phix_c(i/2+1,j/2     ,k/2     ) &
               + three32*phix_c(i/2+1,j/2+joff,k/2     ) &
               + three32*phix_c(i/2+1,j/2     ,k/2+koff) &
               + one32  *phix_c(i/2+1,j/2+joff,k/2+koff)
       end if

    end do
    end do
    end do

    do k=lo_f(3),hi_f(3)
    do j=lo_f(2),hi_f(2)+1
    do i=lo_f(1),hi_f(1)

       if (mod(i,2) .eq. 0) then
          ioff = -1
       else
          ioff = 1
       end if

       if (mod(k,2) .eq. 0) then
          koff = -1
       else
          koff = 1
       end if

       if (mod(j,2) .eq. 0) then
          ! bilinear in the xz plane
          phiy_f(i,j,k) = phiy_f(i,j,k) &
               + nine16* phiy_c(i/2     ,j/2,k/2     ) &
               + three16*phiy_c(i/2+ioff,j/2,k/2     ) &
               + three16*phiy_c(i/2     ,j/2,k/2+koff) &
               + one16*  phiy_c(i/2+ioff,j/2,k/2+koff)
       else
          ! bilinear in the yz plane, linear in y
          phiy_f(i,j,k) = phiy_f(i,j,k) &
               + nine32* phiy_c(i/2     ,j/2  ,k/2     ) &
               + three32*phiy_c(i/2+ioff,j/2  ,k/2     ) &
               + three32*phiy_c(i/2     ,j/2  ,k/2+koff) &
               + one32*  phiy_c(i/2+ioff,j/2  ,k/2+koff) &
               + nine32* phiy_c(i/2     ,j/2+1,k/2     ) &
               + three32*phiy_c(i/2+ioff,j/2+1,k/2     ) &
               + three32*phiy_c(i/2     ,j/2+1,k/2+koff) &
               + one32*  phiy_c(i/2+ioff,j/2+1,k/2+koff)
       end if

    end do
    end do
    end do

    do k=lo_f(3),hi_f(3)+1
    do j=lo_f(2),hi_f(2)
    do i=lo_f(1),hi_f(1)

       if (mod(i,2) .eq. 0) then
          ioff = -1
       else
          ioff = 1
       end if

       if (mod(j,2) .eq. 0) then
          joff = -1
       else
          joff = 1
       end if

       if (mod(k,2) .eq. 0) then
          ! bilinear in the xy plane
          phiz_f(i,j,k) = phiz_f(i,j,k) &
               + nine16* phiz_c(i/2     ,j/2     ,k/2) &
               + three16*phiz_c(i/2+ioff,j/2     ,k/2) &
               + three16*phiz_c(i/2     ,j/2+joff,k/2) &
               + one16*  phiz_c(i/2+ioff,j/2+joff,k/2)
       else
          ! bilinear in the xy plane, linear in z
          phiz_f(i,j,k) = phiz_f(i,j,k) &
               + nine32* phiz_c(i/2     ,j/2     ,k/2  ) &
               + three32*phiz_c(i/2+ioff,j/2     ,k/2  ) &
               + three32*phiz_c(i/2     ,j/2+joff,k/2  ) &
               + one32*  phiz_c(i/2+ioff,j/2+joff,k/2  ) &
               + nine32* phiz_c(i/2     ,j/2     ,k/2+1) &
               + three32*phiz_c(i/2+ioff,j/2     ,k/2+1) &
               + three32*phiz_c(i/2     ,j/2+joff,k/2+1) &
               + one32*  phiz_c(i/2+ioff,j/2+joff,k/2+1)
       end if

    end do
    end do
    end do

#endif

  end subroutine stag_prolongation

  subroutine stag_mg_update(lo,hi, &
                            phix, px_lo, px_hi, &
                            phiy, py_lo, py_hi, &
#if (AMREX_SPACEDIM == 3)
                            phiz, pz_lo, pz_hi, &
#endif
                            rhsx, rx_lo, rx_hi, &
                            rhsy, ry_lo, ry_hi, &
#if (AMREX_SPACEDIM == 3)
                            rhsz, rz_lo, rz_hi, &
#endif
                            Lpx, Lx_lo, Lx_hi, &
                            Lpy, Ly_lo, Ly_hi, &
#if (AMREX_SPACEDIM == 3)
                            Lpz, Lz_lo, Lz_hi, &
#endif
                            alphax, ax_lo, ax_hi, &
                            alphay, ay_lo, ay_hi, &
#if (AMREX_SPACEDIM == 3)
                            alphaz, az_lo, az_hi, &
#endif
                            beta, b_lo, b_hi, &
#if (AMREX_SPACEDIM == 2)
                            beta_ed, w_lo, w_hi, &
#elif (AMREX_SPACEDIM == 3)
                            beta_xy, w_lo, w_hi, &
                            beta_xz, x_lo, x_hi, &
                            beta_yz, y_lo, y_hi, &
#endif
                            gamma, g_lo, g_hi, &
                            dx, color) &
                            bind (C,name="stag_mg_update")

    integer         , intent(in   ) :: lo(3), hi(3)
    integer         , intent(in   ) :: px_lo(3), px_hi(3)
    double precision, intent(inout) :: phix(px_lo(1):px_hi(1),px_lo(2):px_hi(2),px_lo(3):px_hi(3))
    integer         , intent(in   ) :: py_lo(3), py_hi(3)
    double precision, intent(inout) :: phiy(py_lo(1):py_hi(1),py_lo(2):py_hi(2),py_lo(3):py_hi(3))
#if (AMREX_SPACEDIM == 3)
    integer         , intent(in   ) :: pz_lo(3), pz_hi(3)
    double precision, intent(inout) :: phiz(pz_lo(1):pz_hi(1),pz_lo(2):pz_hi(2),pz_lo(3):pz_hi(3))
#endif
    integer         , intent(in   ) :: rx_lo(3), rx_hi(3)
    double precision, intent(in   ) :: rhsx(rx_lo(1):rx_hi(1),rx_lo(2):rx_hi(2),rx_lo(3):rx_hi(3))
    integer         , intent(in   ) :: ry_lo(3), ry_hi(3)
    double precision, intent(in   ) :: rhsy(ry_lo(1):ry_hi(1),ry_lo(2):ry_hi(2),ry_lo(3):ry_hi(3))
#if (AMREX_SPACEDIM == 3)
    integer         , intent(in   ) :: rz_lo(3), rz_hi(3)
    double precision, intent(in   ) :: rhsz(rz_lo(1):rz_hi(1),rz_lo(2):rz_hi(2),rz_lo(3):rz_hi(3))
#endif
    integer         , intent(in   ) :: Lx_lo(3), Lx_hi(3)
    double precision, intent(in   ) :: Lpx(Lx_lo(1):Lx_hi(1),Lx_lo(2):Lx_hi(2),Lx_lo(3):Lx_hi(3))
    integer         , intent(in   ) :: Ly_lo(3), Ly_hi(3)
    double precision, intent(in   ) :: Lpy(Ly_lo(1):Ly_hi(1),Ly_lo(2):Ly_hi(2),Ly_lo(3):Ly_hi(3))
#if (AMREX_SPACEDIM == 3)
    integer         , intent(in   ) :: Lz_lo(3), Lz_hi(3)
    double precision, intent(in   ) :: Lpz(Lz_lo(1):Lz_hi(1),Lz_lo(2):Lz_hi(2),Lz_lo(3):Lz_hi(3))
#endif
    integer         , intent(in   ) :: ax_lo(3), ax_hi(3)
    double precision, intent(in   ) :: alphax(ax_lo(1):ax_hi(1),ax_lo(2):ax_hi(2),ax_lo(3):ax_hi(3))
    integer         , intent(in   ) :: ay_lo(3), ay_hi(3)
    double precision, intent(in   ) :: alphay(ay_lo(1):ay_hi(1),ay_lo(2):ay_hi(2),ay_lo(3):ay_hi(3))
#if (AMREX_SPACEDIM == 3)
    integer         , intent(in   ) :: az_lo(3), az_hi(3)
    double precision, intent(in   ) :: alphaz(az_lo(1):az_hi(1),az_lo(2):az_hi(2),az_lo(3):az_hi(3))
#endif
    integer         , intent(in   ) :: b_lo(3), b_hi(3)
    double precision, intent(in   ) :: beta(b_lo(1):b_hi(1),b_lo(2):b_hi(2),b_lo(3):b_hi(3))
#if (AMREX_SPACEDIM == 2)
    integer         , intent(in   ) :: w_lo(3), w_hi(3)
    double precision, intent(in   ) :: beta_ed(w_lo(1):w_hi(1),w_lo(2):w_hi(2),w_lo(3):w_hi(3))
#elif (AMREX_SPACEDIM == 3)
    integer         , intent(in   ) :: w_lo(3), w_hi(3)
    double precision, intent(in   ) :: beta_xy(w_lo(1):w_hi(1),w_lo(2):w_hi(2),w_lo(3):w_hi(3))
    integer         , intent(in   ) :: x_lo(3), x_hi(3)
    double precision, intent(in   ) :: beta_xz(x_lo(1):x_hi(1),x_lo(2):x_hi(2),x_lo(3):x_hi(3))
    integer         , intent(in   ) :: y_lo(3), y_hi(3)
    double precision, intent(in   ) :: beta_yz(y_lo(1):y_hi(1),y_lo(2):y_hi(2),y_lo(3):y_hi(3))
#endif
    integer         , intent(in   ) :: g_lo(3), g_hi(3)
    double precision, intent(in   ) :: gamma(g_lo(1):g_hi(1),g_lo(2):g_hi(2),g_lo(3):g_hi(3))
    double precision, intent(in   ) :: dx(3)
    integer         , intent(in   ) :: color

#if (AMREX_SPACEDIM == 2)

    ! local
    integer :: i,j,k

    double precision :: fac, dxsq, dxsqinv, fourthirds, fourteenthirds
    double precision :: b,c

    ! coloring parameters
    logical :: do_x, do_y
    integer :: offset, ioff

    do_x = .true.
    do_y = .true.
    offset = 1

    if (color .eq. 1 .or. color .eq. 2) then
       do_y = .false.
       offset = 2
    else if (color .eq. 3 .or. color .eq. 4) then
       do_x = .false.
       offset = 2
    end if

    dxsq = dx(1)**2
    dxsqinv = 1.d0/dxsq
    fourthirds = 4.d0/3.d0
    fourteenthirds = 14.d0/3.d0

    k = 0

    if (visc_type .eq. -1) then

       if (do_x) then

          do j=lo(2),hi(2)
             ioff = 0
             if ( offset .eq. 2 .and. mod(lo(1)+j,2) .ne. mod(color+1,2) ) ioff = 1
             do i=lo(1)+ioff,hi(1)+1,offset

                fac = alphax(i,j,k) + &
                     (beta(i,j,k)+beta(i-1,j,k)+beta_ed(i,j,k)+beta_ed(i,j+1,k)) * dxsqinv

                phix(i,j,k) = phix(i,j,k) + stag_mg_omega*(rhsx(i,j,k)-Lpx(i,j,k)) / fac

             end do
          end do

       end if

       if (do_y) then

          do j=lo(2),hi(2)+1
             ioff = 0
             if ( offset .eq. 2 .and. mod(lo(1)+j,2) .ne. mod(color+1,2) ) ioff = 1
             do i=lo(1)+ioff,hi(1),offset

                fac = alphay(i,j,k) + &
                     (beta(i,j,k)+beta(i,j-1,k)+beta_ed(i,j,k)+beta_ed(i+1,j,k)) * dxsqinv

                phiy(i,j,k) = phiy(i,j,k) + stag_mg_omega*(rhsy(i,j,k)-Lpy(i,j,k)) / fac

             end do
          end do

       end if

    else if (visc_type .eq. 1) then

       b = beta(lo(1),lo(2),k)

       if (do_x) then

          do j=lo(2),hi(2)
             ioff = 0
             if ( offset .eq. 2 .and. mod(lo(1)+j,2) .ne. mod(color+1,2) ) ioff = 1
             do i=lo(1)+ioff,hi(1)+1,offset

                fac = alphax(i,j,k) + 4.d0*b * dxsqinv
                phix(i,j,k) = phix(i,j,k) + stag_mg_omega*(rhsx(i,j,k)-Lpx(i,j,k)) / fac

             end do
          end do

       end if

       if (do_y) then

          do j=lo(2),hi(2)+1
             ioff = 0
             if ( offset .eq. 2 .and. mod(lo(1)+j,2) .ne. mod(color+1,2) ) ioff = 1
             do i=lo(1)+ioff,hi(1),offset

                fac = alphay(i,j,k) + 4.d0*b * dxsqinv
                phiy(i,j,k) = phiy(i,j,k) + stag_mg_omega*(rhsy(i,j,k)-Lpy(i,j,k)) / fac

             end do
          end do

       end if

    else if (visc_type .eq. -2) then

       if (do_x) then

          do j=lo(2),hi(2)
             ioff = 0
             if ( offset .eq. 2 .and. mod(lo(1)+j,2) .ne. mod(color+1,2) ) ioff = 1
             do i=lo(1)+ioff,hi(1)+1,offset

                fac = alphax(i,j,k) + &
                     (2.d0*beta(i,j,k)+2.d0*beta(i-1,j,k)+beta_ed(i,j,k)+beta_ed(i,j+1,k)) * dxsqinv

                phix(i,j,k) = phix(i,j,k) + stag_mg_omega*(rhsx(i,j,k)-Lpx(i,j,k)) / fac

             end do
          end do

       end if

       if (do_y) then

          do j=lo(2),hi(2)+1
             ioff = 0
             if ( offset .eq. 2 .and. mod(lo(1)+j,2) .ne. mod(color+1,2) ) ioff = 1
             do i=lo(1)+ioff,hi(1),offset

                fac = alphay(i,j,k) + &
                     (2.d0*beta(i,j,k)+2.d0*beta(i,j-1,k)+beta_ed(i,j,k)+beta_ed(i+1,j,k)) * dxsqinv

                phiy(i,j,k) = phiy(i,j,k) + stag_mg_omega*(rhsy(i,j,k)-Lpy(i,j,k)) / fac

             end do
          end do

       end if

    else if (visc_type .eq. 2) then

       b = beta(lo(1),lo(2),k)

       if (do_x) then

          do j=lo(2),hi(2)
             ioff = 0
             if ( offset .eq. 2 .and. mod(lo(1)+j,2) .ne. mod(color+1,2) ) ioff = 1
             do i=lo(1)+ioff,hi(1)+1,offset

                fac = alphax(i,j,k)+6.d0*b * dxsqinv
                phix(i,j,k) = phix(i,j,k) + stag_mg_omega*(rhsx(i,j,k)-Lpx(i,j,k)) / fac

             end do
          end do

       end if

       if (do_y) then

          do j=lo(2),hi(2)+1
             ioff = 0
             if ( offset .eq. 2 .and. mod(lo(1)+j,2) .ne. mod(color+1,2) ) ioff = 1
             do i=lo(1)+ioff,hi(1),offset

                fac = alphay(i,j,k)+6.d0*b * dxsqinv
                phiy(i,j,k) = phiy(i,j,k) + stag_mg_omega*(rhsy(i,j,k)-Lpy(i,j,k)) / fac

             end do
          end do

       end if

    else if (visc_type .eq. -3) then

       if (do_x) then

          do j=lo(2),hi(2)
             ioff = 0
             if ( offset .eq. 2 .and. mod(lo(1)+j,2) .ne. mod(color+1,2) ) ioff = 1
             do i=lo(1)+ioff,hi(1)+1,offset

                fac = alphax(i,j,k) + &
                     ( fourthirds*beta(i  ,j,k)+gamma(i,j,k) &
                      +fourthirds*beta(i-1,j,k)+gamma(i-1,j,k) &
                      +beta_ed(i,j,k)+beta_ed(i,j+1,k)) * dxsqinv

                phix(i,j,k) = phix(i,j,k) + stag_mg_omega*(rhsx(i,j,k)-Lpx(i,j,k)) / fac

             end do
          end do

       end if

       if (do_y) then

          do j=lo(2),hi(2)+1
             ioff = 0
             if ( offset .eq. 2 .and. mod(lo(1)+j,2) .ne. mod(color+1,2) ) ioff = 1
             do i=lo(1)+ioff,hi(1),offset

                fac = alphay(i,j,k) + &
                     ( fourthirds*beta(i,j  ,k)+gamma(i,j,k) &
                      +fourthirds*beta(i,j-1,k)+gamma(i,j-1,k) &
                      +beta_ed(i,j,k)+beta_ed(i+1,j,k)) * dxsqinv

                phiy(i,j,k) = phiy(i,j,k) + stag_mg_omega*(rhsy(i,j,k)-Lpy(i,j,k)) / fac

             end do
          end do

       end if

    else if (visc_type .eq. 3) then

       b = beta(lo(1),lo(2),k)
       c = gamma(lo(1),lo(2),k)

       if (do_x) then

          do j=lo(2),hi(2)
             ioff = 0
             if ( offset .eq. 2 .and. mod(lo(1)+j,2) .ne. mod(color+1,2) ) ioff = 1
             do i=lo(1)+ioff,hi(1)+1,offset

                fac = alphax(i,j,k) + (fourteenthirds*b+2.d0*c) * dxsqinv
                phix(i,j,k) = phix(i,j,k) + stag_mg_omega*(rhsx(i,j,k)-Lpx(i,j,k)) / fac

             end do
          end do

       end if

       if (do_y) then

          do j=lo(2),hi(2)+1
             ioff = 0
             if ( offset .eq. 2 .and. mod(lo(1)+j,2) .ne. mod(color+1,2) ) ioff = 1
             do i=lo(1)+ioff,hi(1),offset

                fac = alphay(i,j,k) + (fourteenthirds*b+2.d0*c) * dxsqinv
                phiy(i,j,k) = phiy(i,j,k) + stag_mg_omega*(rhsy(i,j,k)-Lpy(i,j,k)) / fac

             end do
          end do

       end if

    end if

#elif (AMREX_SPACEDIM == 3)

    ! local
    integer :: i,j,k

    double precision :: fac, dxsq, dxsqinv, fourthirds, twentythirds
    double precision :: b,c

    ! coloring parameters
    logical :: do_x, do_y, do_z
    integer :: offset, ioff

    do_x = .true.
    do_y = .true.
    do_z = .true.
    offset = 1

    if (color .eq. 1 .or. color .eq. 2) then
       do_y = .false.
       do_z = .false.
       offset = 2
    else if (color .eq. 3 .or. color .eq. 4) then
       do_x = .false.
       do_z = .false.
       offset = 2
    else if (color .eq. 5 .or. color .eq. 6) then
       do_x = .false.
       do_y = .false.
       offset = 2
    end if

    dxsq = dx(1)**2
    dxsqinv = 1.d0/dxsq
    fourthirds = 4.d0/3.d0
    twentythirds = 20.d0/3.d0

    if (visc_type .eq. -1) then

       if (do_x) then

          do k=lo(3),hi(3)
             do j=lo(2),hi(2)
                ioff = 0
                if ( offset .eq. 2 .and. mod(lo(1)+j+k,2) .ne. mod(color+1,2) ) ioff = 1
                do i=lo(1)+ioff,hi(1)+1,offset

                   fac = alphax(i,j,k) + &
                        ( beta(i,j,k)+beta(i-1,j,k) &
                        +beta_xy(i,j,k)+beta_xy(i,j+1,k) &
                        +beta_xz(i,j,k)+beta_xz(i,j,k+1) ) * dxsqinv

                   phix(i,j,k) = phix(i,j,k) + stag_mg_omega*(rhsx(i,j,k)-Lpx(i,j,k)) / fac

                end do
             end do
          end do

       end if

       if (do_y) then

          do k=lo(3),hi(3)
             do j=lo(2),hi(2)+1
                ioff = 0
                if ( offset .eq. 2 .and. mod(lo(1)+j+k,2) .ne. mod(color+1,2) ) ioff = 1
                do i=lo(1)+ioff,hi(1),offset

                   fac = alphay(i,j,k) + &
                        ( beta(i,j,k)+beta(i,j-1,k) &
                        +beta_xy(i,j,k)+beta_xy(i+1,j,k) &
                        +beta_yz(i,j,k)+beta_yz(i,j,k+1) ) * dxsqinv

                   phiy(i,j,k) = phiy(i,j,k) + stag_mg_omega*(rhsy(i,j,k)-Lpy(i,j,k)) / fac

                end do
             end do
          end do

       end if

       if (do_z) then

          do k=lo(3),hi(3)+1
             do j=lo(2),hi(2)
                ioff = 0
                if ( offset .eq. 2 .and. mod(lo(1)+j+k,2) .ne. mod(color+1,2) ) ioff = 1
                do i=lo(1)+ioff,hi(1),offset

                   fac = alphaz(i,j,k) + &
                        ( beta(i,j,k)+beta(i,j,k-1) &
                        +beta_xz(i,j,k)+beta_xz(i+1,j,k) &
                        +beta_yz(i,j,k)+beta_yz(i,j+1,k) ) * dxsqinv

                   phiz(i,j,k) = phiz(i,j,k) + stag_mg_omega*(rhsz(i,j,k)-Lpz(i,j,k)) / fac

                end do
             end do
          end do

       end if

    else if (visc_type .eq. 1) then

       b = beta(lo(1),lo(2),lo(3))

       if (do_x) then

          do k=lo(3),hi(3)
             do j=lo(2),hi(2)
                ioff = 0
                if ( offset .eq. 2 .and. mod(lo(1)+j+k,2) .ne. mod(color+1,2) ) ioff = 1
                do i=lo(1)+ioff,hi(1)+1,offset

                   fac = alphax(i,j,k) + 6.d0*b * dxsqinv
                   phix(i,j,k) = phix(i,j,k) + stag_mg_omega*(rhsx(i,j,k)-Lpx(i,j,k)) / fac

                end do
             end do
          end do

       end if

       if (do_y) then

          do k=lo(3),hi(3)
             do j=lo(2),hi(2)+1
                ioff = 0
                if ( offset .eq. 2 .and. mod(lo(1)+j+k,2) .ne. mod(color+1,2) ) ioff = 1
                do i=lo(1)+ioff,hi(1),offset

                   fac = alphay(i,j,k) + 6.d0*b * dxsqinv
                   phiy(i,j,k) = phiy(i,j,k) + stag_mg_omega*(rhsy(i,j,k)-Lpy(i,j,k)) / fac

                end do
             end do
          end do

       end if

       if (do_z) then

          do k=lo(3),hi(3)+1
             do j=lo(2),hi(2)
                ioff = 0
                if ( offset .eq. 2 .and. mod(lo(1)+j+k,2) .ne. mod(color+1,2) ) ioff = 1
                do i=lo(1)+ioff,hi(1),offset

                   fac = alphaz(i,j,k) + 6.d0*b * dxsqinv
                   phiz(i,j,k) = phiz(i,j,k) + stag_mg_omega*(rhsz(i,j,k)-Lpz(i,j,k)) / fac

                end do
             end do
          end do

       end if

    else if (visc_type .eq. -2) then

       if (do_x) then

          do k=lo(3),hi(3)
             do j=lo(2),hi(2)
                ioff = 0
                if ( offset .eq. 2 .and. mod(lo(1)+j+k,2) .ne. mod(color+1,2) ) ioff = 1
                do i=lo(1)+ioff,hi(1)+1,offset

                   fac = alphax(i,j,k) + &
                        ( 2.d0*beta(i,j,k)+2.d0*beta(i-1,j,k) &
                        +beta_xy(i,j,k)+beta_xy(i,j+1,k) &
                        +beta_xz(i,j,k)+beta_xz(i,j,k+1) ) * dxsqinv

                   phix(i,j,k) = phix(i,j,k) + stag_mg_omega*(rhsx(i,j,k)-Lpx(i,j,k)) / fac

                end do
             end do
          end do

       end if

       if (do_y) then

          do k=lo(3),hi(3)
             do j=lo(2),hi(2)+1
                ioff = 0
                if ( offset .eq. 2 .and. mod(lo(1)+j+k,2) .ne. mod(color+1,2) ) ioff = 1
                do i=lo(1)+ioff,hi(1),offset

                   fac = alphay(i,j,k) + &
                        ( 2.d0*beta(i,j,k)+2.d0*beta(i,j-1,k) &
                        +beta_xy(i,j,k)+beta_xy(i+1,j,k) &
                        +beta_yz(i,j,k)+beta_yz(i,j,k+1) ) * dxsqinv

                   phiy(i,j,k) = phiy(i,j,k) + stag_mg_omega*(rhsy(i,j,k)-Lpy(i,j,k)) / fac

                end do
             end do
          end do

       end if

       if (do_z) then

          do k=lo(3),hi(3)+1
             do j=lo(2),hi(2)
                ioff = 0
                if ( offset .eq. 2 .and. mod(lo(1)+j+k,2) .ne. mod(color+1,2) ) ioff = 1
                do i=lo(1)+ioff,hi(1),offset

                   fac = alphaz(i,j,k) + &
                        ( 2.d0*beta(i,j,k)+2.d0*beta(i,j,k-1) &
                        +beta_xz(i,j,k)+beta_xz(i+1,j,k) &
                        +beta_yz(i,j,k)+beta_yz(i,j+1,k) ) * dxsqinv

                   phiz(i,j,k) = phiz(i,j,k) + stag_mg_omega*(rhsz(i,j,k)-Lpz(i,j,k)) / fac

                end do
             end do
          end do

       end if

    else if (visc_type .eq. 2) then

       b = beta(lo(1),lo(2),lo(3))

       if (do_x) then

          do k=lo(3),hi(3)
             do j=lo(2),hi(2)
                ioff = 0
                if ( offset .eq. 2 .and. mod(lo(1)+j+k,2) .ne. mod(color+1,2) ) ioff = 1
                do i=lo(1)+ioff,hi(1)+1,offset

                   fac = alphax(i,j,k) + 8.d0*b * dxsqinv
                   phix(i,j,k) = phix(i,j,k) + stag_mg_omega*(rhsx(i,j,k)-Lpx(i,j,k)) / fac

                end do
             end do
          end do

       end if

       if (do_y) then

          do k=lo(3),hi(3)
             do j=lo(2),hi(2)+1
                ioff = 0
                if ( offset .eq. 2 .and. mod(lo(1)+j+k,2) .ne. mod(color+1,2) ) ioff = 1
                do i=lo(1)+ioff,hi(1),offset

                   fac = alphay(i,j,k) + 8.d0*b * dxsqinv
                   phiy(i,j,k) = phiy(i,j,k) + stag_mg_omega*(rhsy(i,j,k)-Lpy(i,j,k)) / fac

                end do
             end do
          end do

       end if

       if (do_z) then

          do k=lo(3),hi(3)+1
             do j=lo(2),hi(2)
                ioff = 0
                if ( offset .eq. 2 .and. mod(lo(1)+j+k,2) .ne. mod(color+1,2) ) ioff = 1
                do i=lo(1)+ioff,hi(1),offset

                   fac = alphaz(i,j,k) + 8.d0*b * dxsqinv
                   phiz(i,j,k) = phiz(i,j,k) + stag_mg_omega*(rhsz(i,j,k)-Lpz(i,j,k)) / fac

                end do
             end do
          end do

       end if

    else if (visc_type .eq. -3) then

       if (do_x) then

          do k=lo(3),hi(3)
             do j=lo(2),hi(2)
                ioff = 0
                if ( offset .eq. 2 .and. mod(lo(1)+j+k,2) .ne. mod(color+1,2) ) ioff = 1
                do i=lo(1)+ioff,hi(1)+1,offset

                   fac = alphax(i,j,k) + &
                        ( fourthirds*beta(i,j,k)+gamma(i,j,k) &
                        +fourthirds*beta(i-1,j,k)+gamma(i-1,j,k) &
                        +beta_xy(i,j,k)+beta_xy(i,j+1,k) &
                        +beta_xz(i,j,k)+beta_xz(i,j,k+1) ) * dxsqinv

                   phix(i,j,k) = phix(i,j,k) + stag_mg_omega*(rhsx(i,j,k)-Lpx(i,j,k)) / fac

                end do
             end do
          end do

       end if

       if (do_y) then

          do k=lo(3),hi(3)
             do j=lo(2),hi(2)+1
                ioff = 0
                if ( offset .eq. 2 .and. mod(lo(1)+j+k,2) .ne. mod(color+1,2) ) ioff = 1
                do i=lo(1)+ioff,hi(1),offset

                   fac = alphay(i,j,k) + &
                        ( fourthirds*beta(i,j,k)+gamma(i,j,k) &
                        +fourthirds*beta(i,j-1,k)+gamma(i,j-1,k) &
                        +beta_xy(i,j,k)+beta_xy(i+1,j,k) &
                        +beta_yz(i,j,k)+beta_yz(i,j,k+1) ) * dxsqinv

                   phiy(i,j,k) = phiy(i,j,k) + stag_mg_omega*(rhsy(i,j,k)-Lpy(i,j,k)) / fac

                end do
             end do
          end do

       end if

       if (do_z) then

          do k=lo(3),hi(3)+1
             do j=lo(2),hi(2)
                ioff = 0
                if ( offset .eq. 2 .and. mod(lo(1)+j+k,2) .ne. mod(color+1,2) ) ioff = 1
                do i=lo(1)+ioff,hi(1),offset

                   fac = alphaz(i,j,k) + &
                        ( fourthirds*beta(i,j,k)+gamma(i,j,k) &
                        +fourthirds*beta(i,j,k-1)+gamma(i,j,k-1) &
                        +beta_xz(i,j,k)+beta_xz(i+1,j,k) &
                        +beta_yz(i,j,k)+beta_yz(i,j+1,k) ) * dxsqinv

                   phiz(i,j,k) = phiz(i,j,k) + stag_mg_omega*(rhsz(i,j,k)-Lpz(i,j,k)) / fac

                end do
             end do
          end do

       end if

    else if (visc_type .eq. 3) then

       b = beta(lo(1),lo(2),lo(3))
       c = gamma(lo(1),lo(2),lo(3))

       if (do_x) then

          do k=lo(3),hi(3)
             do j=lo(2),hi(2)
                ioff = 0
                if ( offset .eq. 2 .and. mod(lo(1)+j+k,2) .ne. mod(color+1,2) ) ioff = 1
                do i=lo(1)+ioff,hi(1)+1,offset

                   fac = alphax(i,j,k)+(twentythirds*b+2.d0*c) * dxsqinv
                   phix(i,j,k) = phix(i,j,k) + stag_mg_omega*(rhsx(i,j,k)-Lpx(i,j,k)) / fac

                end do
             end do
          end do

       end if

       if (do_y) then

          do k=lo(3),hi(3)
             do j=lo(2),hi(2)+1
                ioff = 0
                if ( offset .eq. 2 .and. mod(lo(1)+j+k,2) .ne. mod(color+1,2) ) ioff = 1
                do i=lo(1)+ioff,hi(1),offset

                   fac = alphay(i,j,k)+(twentythirds*b+2.d0*c) * dxsqinv
                   phiy(i,j,k) = phiy(i,j,k) + stag_mg_omega*(rhsy(i,j,k)-Lpy(i,j,k)) / fac

                end do
             end do
          end do

       end if

       if (do_z) then

          do k=lo(3),hi(3)+1
             do j=lo(2),hi(2)
                ioff = 0
                if ( offset .eq. 2 .and. mod(lo(1)+j+k,2) .ne. mod(color+1,2) ) ioff = 1
                do i=lo(1)+ioff,hi(1),offset

                   fac = alphaz(i,j,k)+(twentythirds*b+2.d0*c) * dxsqinv
                   phiz(i,j,k) = phiz(i,j,k) + stag_mg_omega*(rhsz(i,j,k)-Lpz(i,j,k)) / fac

                end do
             end do
          end do

       end if

    end if

#endif

  end subroutine stag_mg_update

end module stag_solver_module
