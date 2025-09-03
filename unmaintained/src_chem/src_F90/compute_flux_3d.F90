module compute_flux_module_3d

  implicit none

  private

  public :: compute_flux_3d

contains

  subroutine compute_flux_3d(lo, hi, dt, dx, &
                             con,ph_lo,ph_hi, &
                             umac,  u_lo,  u_hi, &
                             vmac,  v_lo,  v_hi, &
                             wmac,  w_lo,  w_hi, &
                             flxx, fx_lo, fx_hi, &
                             flxy, fy_lo, fy_hi, &
                             flxz, fz_lo, fz_hi, &
                             conx, conx_y, conx_z, &
                             cony, cony_x, cony_z, &
                             conz, conz_x, conz_y, &
                             slope, glo, ghi,nu)

    use slope_module_3d, only: slopex_3d, slopey_3d, slopez_3d

    integer, intent(in) :: lo(3), hi(3), glo(3), ghi(3)
    double precision, intent(in) :: dt, dx(3), nu
    integer, intent(in) :: ph_lo(3), ph_hi(3)
    integer, intent(in) ::  u_lo(3),  u_hi(3)
    integer, intent(in) ::  v_lo(3),  v_hi(3)
    integer, intent(in) ::  w_lo(3),  w_hi(3)
    integer, intent(in) :: fx_lo(3), fx_hi(3)
    integer, intent(in) :: fy_lo(3), fy_hi(3)
    integer, intent(in) :: fz_lo(3), fz_hi(3)
    double precision, intent(in   ) :: con (ph_lo(1):ph_hi(1),ph_lo(2):ph_hi(2),ph_lo(3):ph_hi(3))
    double precision, intent(in   ) :: umac( u_lo(1): u_hi(1), u_lo(2): u_hi(2), u_lo(3): u_hi(3))
    double precision, intent(in   ) :: vmac( v_lo(1): v_hi(1), v_lo(2): v_hi(2), v_lo(3): v_hi(3))
    double precision, intent(in   ) :: wmac( w_lo(1): w_hi(1), w_lo(2): w_hi(2), w_lo(3): w_hi(3))
    double precision, intent(  out) :: flxx(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2),fx_lo(3):fx_hi(3))
    double precision, intent(  out) :: flxy(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2),fy_lo(3):fy_hi(3))
    double precision, intent(  out) :: flxz(fz_lo(1):fz_hi(1),fz_lo(2):fz_hi(2),fz_lo(3):fz_hi(3))
    double precision, dimension(glo(1):ghi(1),glo(2):ghi(2),glo(3):ghi(3)) :: &
         conx, conx_y, conx_z, cony, cony_x, cony_z, conz, conz_x, conz_y, slope

    integer :: i, j, k
    double precision :: hdtdx(3), tdtdx(3)

    hdtdx = 0.5*(dt/dx)
    tdtdx = (1.d0/3.d0)*(dt/dx)

    call slopex_3d(glo, ghi, &
                con, ph_lo, ph_hi, &
                slope, glo, ghi)

    ! compute con on x faces using umac to upwind; ignore transverse terms
    do       k = lo(3)-1, hi(3)+1
       do    j = lo(2)-1, hi(2)+1
          do i = lo(1)  , hi(1)+1

             if (umac(i,j,k) .lt. 0.d0) then
                conx(i,j,k) = con(i  ,j,k) - (0.5d0 + hdtdx(1)*umac(i,j,k))*slope(i  ,j,k)
             else
                conx(i,j,k) = con(i-1,j,k) + (0.5d0 - hdtdx(1)*umac(i,j,k))*slope(i-1,j,k)
             end if

          end do
       end do
    end do

    call slopey_3d(glo, ghi, &
                con, ph_lo, ph_hi, &
                slope, glo, ghi)

    ! compute con on y faces using vmac to upwind; ignore transverse terms
    do       k = lo(3)-1, hi(3)+1
       do    j = lo(2)  , hi(2)+1
          do i = lo(1)-1, hi(1)+1

             if (vmac(i,j,k) .lt. 0.d0) then
                cony(i,j,k) = con(i,j  ,k) - (0.5d0 + hdtdx(2)*vmac(i,j,k))*slope(i,j  ,k)
             else
                cony(i,j,k) = con(i,j-1,k) + (0.5d0 - hdtdx(2)*vmac(i,j,k))*slope(i,j-1,k)
             end if

          end do
       end do
    end do

    call slopez_3d(glo, ghi, &
                con, ph_lo, ph_hi, &
                slope, glo, ghi)

    ! compute con on z faces using wmac to upwind; ignore transverse terms
    do       k = lo(3)  , hi(3)+1
       do    j = lo(2)-1, hi(2)+1
          do i = lo(1)-1, hi(1)+1

             if (wmac(i,j,k) .lt. 0.d0) then
                conz(i,j,k) = con(i,j,k  ) - (0.5d0 + hdtdx(3)*wmac(i,j,k))*slope(i,j,k  )
             else
                conz(i,j,k) = con(i,j,k-1) + (0.5d0 - hdtdx(3)*wmac(i,j,k))*slope(i,j,k-1)
             end if

          end do
       end do
    end do

    !!!!!!!!!!!!!!!!!!!!
    ! transverse terms
    !!!!!!!!!!!!!!!!!!!!

    ! update con on x faces by adding in y-transverse terms
    do       k=lo(3)-1, hi(3)+1
       do    j=lo(2)  , hi(2)
          do i=lo(1)  , hi(1)+1

             if (umac(i,j,k) .lt. 0.d0) then
                conx_y(i,j,k) = conx(i,j,k) &
                     - tdtdx(2) * (0.5d0*(vmac(i  ,j+1,k)+vmac(i  ,j,k)) * (cony(i  ,j+1,k)-cony(i  ,j,k)) )
             else
                conx_y(i,j,k) = conx(i,j,k) &
                     - tdtdx(2) * (0.5d0*(vmac(i-1,j+1,k)+vmac(i-1,j,k)) * (cony(i-1,j+1,k)-cony(i-1,j,k)) )
             end if

          end do
       end do
    end do

    ! update con on x faces by adding in z-transverse terms
    do       k=lo(3)  , hi(3)
       do    j=lo(2)-1, hi(2)+1
          do i=lo(1)  , hi(1)+1

             if (umac(i,j,k) .lt. 0.d0) then
                conx_z(i,j,k) = conx(i,j,k) &
                     - tdtdx(3) * (0.5d0*(wmac(i  ,j,k+1)+wmac(i  ,j,k)) * (conz(i  ,j,k+1)-conz(i  ,j,k)) )
             else
                conx_z(i,j,k) = conx(i,j,k) &
                     - tdtdx(3) * (0.5d0*(wmac(i-1,j,k+1)+wmac(i-1,j,k)) * (conz(i-1,j,k+1)-conz(i-1,j,k)) )
             end if

          end do
       end do
    end do

    ! update con on y faces by adding in x-transverse terms
    do       k = lo(3)-1, hi(3)+1
       do    j = lo(2)  , hi(2)+1
          do i = lo(1)  , hi(1)

             if (vmac(i,j,k) .lt. 0.d0) then
                cony_x(i,j,k) = cony(i,j,k) &
                     - tdtdx(1) * (0.5d0*(umac(i+1,j  ,k)+umac(i,j  ,k)) * (conx(i+1,j  ,k)-conx(i,j  ,k)) )
             else
                cony_x(i,j,k) = cony(i,j,k) &
                     - tdtdx(1) * (0.5d0*(umac(i+1,j-1,k)+umac(i,j-1,k)) * (conx(i+1,j-1,k)-conx(i,j-1,k)) )
             end if

          end do
       end do
    end do

    ! update con on y faces by adding in z-transverse terms
    do       k = lo(3)  , hi(3)
       do    j = lo(2)  , hi(2)+1
          do i = lo(1)-1, hi(1)+1

             if (vmac(i,j,k) .lt. 0.d0) then
                cony_z(i,j,k) = cony(i,j,k) &
                     - tdtdx(3) * (0.5d0*(wmac(i,j  ,k+1)+wmac(i,j  ,k)) * (conz(i,j  ,k+1)-conz(i,j  ,k)) )
             else
                cony_z(i,j,k) = cony(i,j,k) &
                     - tdtdx(3) * (0.5d0*(wmac(i,j-1,k+1)+wmac(i,j-1,k)) * (conz(i,j-1,k+1)-conz(i,j-1,k)) )
             end if

          end do
       end do
    end do

    ! update con on z faces by adding in x-transverse terms
    do       k = lo(3)  , hi(3)+1
       do    j = lo(2)-1, hi(2)+1
          do i = lo(1)  , hi(1)

             if (wmac(i,j,k) .lt. 0.d0) then
                conz_x(i,j,k) = conz(i,j,k) &
                     - tdtdx(1) * (0.5d0*(umac(i+1,j,k  )+umac(i,j,k  )) * (conx(i+1,j,k  )-conx(i,j,k  )) )
             else
                conz_x(i,j,k) = conz(i,j,k) &
                     - tdtdx(1) * (0.5d0*(umac(i+1,j,k-1)+umac(i,j,k-1)) * (conx(i+1,j,k-1)-conx(i,j,k-1)) )
             end if

          end do
       end do
    end do

    ! update con on z faces by adding in y-transverse terms
    do       k = lo(3)  , hi(3)+1
       do    j = lo(2)  , hi(2)
          do i = lo(1)-1, hi(1)+1

             if (wmac(i,j,k) .lt. 0.d0) then
                conz_y(i,j,k) = conz(i,j,k) &
                     - tdtdx(2) * (0.5d0*(vmac(i,j+1,k  )+vmac(i,j,k  )) * (cony(i,j+1,k  )-cony(i,j,k  )) )
             else
                conz_y(i,j,k) = conz(i,j,k) &
                     - tdtdx(2) * (0.5d0*(vmac(i,j+1,k-1)+vmac(i,j,k-1)) * (cony(i,j+1,k-1)-cony(i,j,k-1)) )
             end if

          end do
       end do
    end do

    !!!!!!!!!!!!!!!!!!!!
    ! final edge states
    !!!!!!!!!!!!!!!!!!!!

    ! update con on x faces by adding in yz and zy transverse terms
    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)+1

             if (umac(i,j,k) .lt. 0.d0) then
                conx(i,j,k) = conx(i,j,k) &
                     - hdtdx(2)*( 0.5d0*(vmac(i  ,j+1,k  )+vmac(i  ,j,k)) * (cony_z(i  ,j+1,k  )-cony_z(i  ,j,k)) ) &
                     - hdtdx(3)*( 0.5d0*(wmac(i  ,j  ,k+1)+wmac(i  ,j,k)) * (conz_y(i  ,j  ,k+1)-conz_y(i  ,j,k)) )
             else
                conx(i,j,k) = conx(i,j,k) &
                     - hdtdx(2)*( 0.5d0*(vmac(i-1,j+1,k  )+vmac(i-1,j,k)) * (cony_z(i-1,j+1,k  )-cony_z(i-1,j,k)) ) &
                     - hdtdx(3)*( 0.5d0*(wmac(i-1,j  ,k+1)+wmac(i-1,j,k)) * (conz_y(i-1,j  ,k+1)-conz_y(i-1,j,k)) )
             end if

             ! compute final x-fluxes
             flxx(i,j,k) = umac(i,j,k)*conx(i,j,k)-nu*(con(i,j,k)-con(i-1,j,k))/dx(1)

          end do
       end do
    end do

    ! update con on y faces by adding in xz and zx transverse terms
    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)+1
          do i = lo(1), hi(1)

             if (vmac(i,j,k) .lt. 0.d0) then
                cony(i,j,k) = cony(i,j,k) &
                     - hdtdx(1)*( 0.5d0*(umac(i+1,j  ,k  )+umac(i,j  ,k)) * (conx_z(i+1,j  ,k  )-conx_z(i,j  ,k)) ) &
                     - hdtdx(3)*( 0.5d0*(wmac(i  ,j  ,k+1)+wmac(i,j  ,k)) * (conz_x(i  ,j  ,k+1)-conz_x(i,j  ,k)) )
             else
                cony(i,j,k) = cony(i,j,k) &
                     - hdtdx(1)*( 0.5d0*(umac(i+1,j-1,k  )+umac(i,j-1,k)) * (conx_z(i+1,j-1,k  )-conx_z(i,j-1,k)) ) &
                     - hdtdx(3)*( 0.5d0*(wmac(i  ,j-1,k+1)+wmac(i,j-1,k)) * (conz_x(i  ,j-1,k+1)-conz_x(i,j-1,k)) )
             end if

             ! compute final y-fluxes
             flxy(i,j,k) = vmac(i,j,k)*cony(i,j,k)-nu*(con(i,j,k)-con(i,j-1,k))/dx(2)


          end do
       end do
    end do

    ! update con on z faces by adding in xy and yx transverse terms
    do       k = lo(3), hi(3)+1
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)

             if (wmac(i,j,k) .lt. 0.d0) then
                conz(i,j,k) = conz(i,j,k) &
                     - hdtdx(1)*( 0.5d0*(umac(i+1,j  ,k  )+umac(i  ,j,k)) * (conx_y(i+1,j  ,k  )-conx_y(i,j,k  )) ) &
                     - hdtdx(2)*( 0.5d0*(vmac(i  ,j+1,k  )+vmac(i  ,j,k)) * (cony_x(i  ,j+1,k  )-cony_x(i,j,k  )) )
             else
                conz(i,j,k) = conz(i,j,k) &
                     - hdtdx(1)*( 0.5d0*(umac(i+1,j  ,k-1)+umac(i,j,k-1)) * (conx_y(i+1,j  ,k-1)-conx_y(i,j,k-1)) ) &
                     - hdtdx(2)*( 0.5d0*(vmac(i  ,j+1,k-1)+vmac(i,j,k-1)) * (cony_x(i  ,j+1,k-1)-cony_x(i,j,k-1)) )
             end if

             ! compute final z-fluxes
             flxz(i,j,k) = wmac(i,j,k)*conz(i,j,k)-nu/dx(3)*(con(i,j,k)-con(i,j,k-1))


          end do
       end do
    end do


  end subroutine compute_flux_3d

end module compute_flux_module_3d
