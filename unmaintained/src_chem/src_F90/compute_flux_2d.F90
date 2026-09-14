module compute_flux_module_2d

  implicit none

  private

  public :: compute_flux_2d

contains

  subroutine compute_flux_2d(lo, hi, dt, dx, &
                             con,ph_lo,ph_hi, &
                             umac,  u_lo,  u_hi, &
                             vmac,  v_lo,  v_hi, &
                             flxx, fx_lo, fx_hi, &
                             flxy, fy_lo, fy_hi, &
                             conx_1d, cony_1d, conx, cony, slope, glo, ghi)

    use slope_module_2d, only: slopex_2d, slopey_2d

    integer, intent(in) :: lo(2), hi(2), glo(2), ghi(2)
    double precision, intent(in) :: dt, dx(2)
    integer, intent(in) :: ph_lo(2), ph_hi(2)
    integer, intent(in) ::  u_lo(2),  u_hi(2)
    integer, intent(in) ::  v_lo(2),  v_hi(2)
    integer, intent(in) :: fx_lo(2), fx_hi(2)
    integer, intent(in) :: fy_lo(2), fy_hi(2)
    double precision, intent(in   ) :: con (ph_lo(1):ph_hi(1),ph_lo(2):ph_hi(2))
    double precision, intent(in   ) :: umac( u_lo(1): u_hi(1), u_lo(2): u_hi(2))
    double precision, intent(in   ) :: vmac( v_lo(1): v_hi(1), v_lo(2): v_hi(2))
    double precision, intent(  out) :: flxx(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2))
    double precision, intent(  out) :: flxy(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2))
    double precision, dimension(glo(1):ghi(1),glo(2):ghi(2)) :: &
         conx_1d, cony_1d, conx, cony, slope

    integer :: i, j, k
    double precision :: hdtdx(2)

    hdtdx = 0.5*(dt/dx)

    call slopex_2d(glo, ghi, &
                con, ph_lo, ph_hi, &
                slope, glo, ghi)

    ! compute con on x faces using umac to upwind; ignore transverse terms
    do    j = lo(2)-1, hi(2)+1
       do i = lo(1)  , hi(1)+1

          if (umac(i,j) .lt. 0.d0) then
             conx_1d(i,j) = con(i  ,j) - (0.5d0 + hdtdx(1)*umac(i,j))*slope(i  ,j)
          else
             conx_1d(i,j) = con(i-1,j) + (0.5d0 - hdtdx(1)*umac(i,j))*slope(i-1,j)
          end if

       end do
    end do

    call slopey_2d(glo, ghi, &
                con, ph_lo, ph_hi, &
                slope, glo, ghi)

    ! compute con on y faces using umac to upwind; ignore transverse terms
    do    j = lo(2)  , hi(2)+1
       do i = lo(1)-1, hi(1)+1

          if (vmac(i,j) .lt. 0.d0) then
             cony_1d(i,j) = con(i,j  ) - (0.5d0 + hdtdx(2)*vmac(i,j))*slope(i,j  )
          else
             cony_1d(i,j) = con(i,j-1) + (0.5d0 - hdtdx(2)*vmac(i,j))*slope(i,j-1)
          end if

       end do
    end do

    ! update con on x faces by adding in y-transverse terms
    do    j = lo(2), hi(2)
       do i = lo(1), hi(1)+1

          if (umac(i,j) .lt. 0.d0) then
             conx(i,j) = conx_1d(i,j) &
                  - hdtdx(2)*( 0.5d0*(vmac(i  ,j+1)+vmac(i  ,j)) * (cony_1d(i  ,j+1)-cony_1d(i  ,j)) )
          else
             conx(i,j) = conx_1d(i,j) &
                  - hdtdx(2)*( 0.5d0*(vmac(i-1,j+1)+vmac(i-1,j)) * (cony_1d(i-1,j+1)-cony_1d(i-1,j)) )
          end if

          ! compute final x-fluxes
          flxx(i,j) = conx(i,j)*umac(i,j)

       end do
    end do

    ! update con on y faces by adding in x-transverse terms
    do    j = lo(2), hi(2)+1
       do i = lo(1), hi(1)

          if (vmac(i,j) .lt. 0.d0) then
             cony(i,j) = cony_1d(i,j) &
                  - hdtdx(1)*( 0.5d0*(umac(i+1,j  )+umac(i,j  )) * (conx_1d(i+1,j  )-conx_1d(i,j  )) )
          else
             cony(i,j) = cony_1d(i,j) &
                  - hdtdx(1)*( 0.5d0*(umac(i+1,j-1)+umac(i,j-1)) * (conx_1d(i+1,j-1)-conx_1d(i,j-1)) )
          end if

          ! compute final y-fluxes
          flxy(i,j) = cony(i,j)*vmac(i,j)

       end do
    end do

  end subroutine compute_flux_2d

end module compute_flux_module_2d
