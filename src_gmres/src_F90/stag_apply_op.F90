module stag_apply_op_module

  use amrex_error_module
  use common_namelist_module, only: visc_type

  implicit none

  private

contains

#if (AMREX_SPACEDIM == 2)

  subroutine stag_apply_op(lo, hi, &
                           betacc, betacclo, betacchi, &
                           gammacc, gammacclo, gammacchi, &
                           betanodal, betanodallo, betanodalhi, &  
                           velxin, velxinlo, velxinhi, &
                           velyin, velyinlo, velyinhi, &
                           velxout, velxoutlo, velxouthi, &
                           velyout, velyoutlo, velyouthi, &
                           alphax, alphaxlo, alphaxhi, &
                           alphay, alphaylo, alphayhi, &
                           dx,color) &
                           bind (C,name="stag_apply_op")

    integer         , intent(in   ) :: lo(2), hi(2), betacclo(2), betacchi(2), gammacclo(2), gammacchi(2), betanodallo(2), betanodalhi(2)
    integer         , intent(in   ) :: alphaxlo(2), alphaxhi(2), alphaylo(2), alphayhi(2)
    integer         , intent(in   ) :: velxinlo(2), velxinhi(2), velyinlo(2), velyinhi(2), velxoutlo(2), velxouthi(2), velyoutlo(2), velyouthi(2)
    double precision, intent(in   ) :: dx(2)
    double precision, intent(in   ) :: betacc(betacclo(1):betacchi(1),betacclo(2):betacchi(2))
    double precision, intent(in   ) :: gammacc(gammacclo(1):gammacchi(1),gammacclo(2):gammacchi(2))
    double precision, intent(in   ) :: betanodal(betanodallo(1):betanodalhi(1),betanodallo(2):betanodalhi(2))
    double precision, intent(in   ) :: alphax(alphaxlo(1):alphaxhi(1),alphaxlo(2):alphaxhi(2))
    double precision, intent(in   ) :: alphay(alphaylo(1):alphayhi(1),alphaylo(2):alphayhi(2))
    double precision, intent(in   ) :: velxin(velxinlo(1):velxinhi(1),velxinlo(2):velxinhi(2))
    double precision, intent(in   ) :: velyin(velyinlo(1):velyinhi(1),velyinlo(2):velyinhi(2))
    double precision, intent(inout) :: velxout(velxoutlo(1):velxouthi(1),velxoutlo(2):velxouthi(2))
    double precision, intent(inout) :: velyout(velyoutlo(1):velyouthi(1),velyoutlo(2):velyouthi(2))
    integer         , intent(in   ) :: color

    ! local
    integer :: i,j
    double precision dxsqinv, dysqinv, dxdyinv, term1, term2, term3
    double precision bt, gm

    bt = betacc(betacclo(1),betacclo(2))
    gm = gammacc(gammacclo(1),gammacclo(2))

    !note that operators are implemented as in FHDfortran, i.e. the negative of the operator

    dxsqinv = 1.d0/(dx(1)**2)
    dysqinv = 1.d0/(dx(2)**2)
    dxdyinv = 1.d0/(dx(1)*dx(2))

    !Type 1

    if (visc_type .eq. -1) then
      do j = lo(2), hi(2)
        do i = lo(1), hi(1)+1

          velxout(i,j) = +velxin(i,j)*(alphax(i,j)+(betacc(i-1,j)+betacc(i,j))*dxsqinv+(betanodal(i,j)+betanodal(i,j+1))*dysqinv) &
                         +(-velxin(i+1,j)*betacc(i,j) &
                         -velxin(i-1,j)*betacc(i-1,j))*dxsqinv &
                         +(-velxin(i,j+1)*betanodal(i,j+1) &
                         -velxin(i,j-1)*betanodal(i,j))*dysqinv
        enddo
      enddo

      do j = lo(2), hi(2)+1
        do i = lo(1), hi(1)

          velyout(i,j) = +velyin(i,j)*(alphax(i,j)+(betacc(i,j)+betacc(i,j-1))*dysqinv+(betanodal(i+1,j)+betanodal(i,j))*dxsqinv) &
                         +(-velyin(i,j+1)*betacc(i,j) &
                         -velyin(i,j-1)*betacc(i,j-1))*dysqinv &
                         +(-velyin(i+1,j)*betanodal(i+1,j) &
                         -velyin(i-1,j)*betanodal(i,j))*dxsqinv
        enddo
      enddo



    endif

    if (visc_type .eq. 1) then

      term1 = 2.d0*bt*(dxsqinv+dysqinv)
      term2 = bt*dxsqinv
      term3 = bt*dysqinv

      do j = lo(2), hi(2)
        do i = lo(1), hi(1)+1

          velxout(i,j) = +velxin(i,j)*(alphax(i,j) + term1) &
                         +(-velxin(i+1,j) &
                           -velxin(i-1,j))*term2 &
                         +(-velxin(i,j+1) &
                           -velxin(i,j-1))*term3

        enddo
      enddo

      do j = lo(2), hi(2)+1
        do i = lo(1), hi(1)

          velyout(i,j) = +velyin(i,j)*(alphay(i,j) + term1) &
                         +(-velyin(i,j+1) &
                           -velyin(i,j-1))*term3 &
                         +(-velyin(i+1,j)&
                           -velyin(i-1,j))*term2
        enddo
      enddo

    endif

    !Type 2

    if (visc_type .eq. -2) then

      do j = lo(2), hi(2)
        do i = lo(1), hi(1)+1

          velxout(i,j) = +velxin(i,j)*(alphax(i,j) + (betacc(i,j)+betacc(i-1,j))*2.d0*dxsqinv+(betanodal(i,j+1)+betanodal(i,j))*dysqinv) &
                     
                         -( (2.d0*velxin(i+1,j)*betacc(i,j) &
                           +2.d0*velxin(i-1,j)*betacc(i-1,j))*dxsqinv &
                           +(+velxin(i,j+1)*betanodal(i,j+1) &
                           +velxin(i,j-1)*betanodal(i,j))*dysqinv &
                     
                           +(+velyin(i,j+1)*betanodal(i,j+1) &
                           -velyin(i,j)*betanodal(i,j) &
                           -velyin(i-1,j+1)*betanodal(i,j+1) &
                           +velyin(i-1,j)*betanodal(i,j))*dxdyinv)
        enddo
      enddo

      do j = lo(2), hi(2)+1
        do i = lo(1), hi(1)

          velyout(i,j) = +velyin(i,j)*(alphay(i,j) + (betacc(i,j)+betacc(i,j-1))*2.d0*dysqinv+(betanodal(i+1,j)+betanodal(i,j))*dxsqinv) &
                     
                         -( (2.d0*velyin(i,j+1)*betacc(i,j) &
                           +2.d0*velyin(i,j-1)*betacc(i,j-1))*dysqinv &
                           +(+velyin(i+1,j)*betanodal(i+1,j) &
                           +velyin(i-1,j)*betanodal(i,j))*dxsqinv &
                     
                           +(+velxin(i+1,j)*betanodal(i+1,j) &
                           -velxin(i,j)*betanodal(i,j) &
                           -velxin(i+1,j-1)*betanodal(i+1,j) &
                           +velxin(i,j-1)*betanodal(i,j))*dxdyinv)
        enddo
      enddo
    endif

    if (visc_type .eq. 2) then

      term1 = 2.d0*bt*(2.d0*dxsqinv+dysqinv)

      do j = lo(2), hi(2)
        do i = lo(1), hi(1)+1

          velxout(i,j) = +velxin(i,j)*(alphax(i,j) + term1) &
                     
                         -( bt*dxsqinv*(2.d0*velxin(i+1,j) &
                                       +2.d0*velxin(i-1,j) &
                                       +velxin(i,j+1) &
                                       +velxin(i,j-1)) &
                     
                           +bt*dxdyinv*(+velyin(i,j+1) &
                                       -velyin(i,j) &
                                       -velyin(i-1,j+1) &
                                       +velyin(i-1,j)))
        enddo
      enddo

      term1 = 2.d0*bt*(2.d0*dysqinv+dxsqinv)

      do j = lo(2), hi(2)+1
        do i = lo(1), hi(1)

          velyout(i,j) = +velyin(i,j)*(alphay(i,j) + term1) &
                     
                         -( bt*dysqinv*(2.d0*velyin(i+1,j) &
                                       +2.d0*velyin(i-1,j) &
                                       +velyin(i,j+1) &
                                       +velyin(i,j-1)) &
                     
                           +bt*dxdyinv*(+velxin(i,j+1) &
                                       -velxin(i,j) &
                                       -velxin(i-1,j+1) &
                                       +velxin(i-1,j)))
        enddo
      enddo
    endif
  
  end subroutine stag_apply_op


#endif

#if (AMREX_SPACEDIM == 3)

  subroutine stag_apply_op(lo, hi, &
                           betacc, betacclo, betacchi, &
                           gammacc, gammacclo, gammacchi, &
                           betaxy, betaxylo, betaxyhi, &  
                           betaxz, betaxzlo, betaxzhi, &  
                           betayz, betayzlo, betayzhi, &  
                           velxin, velxinlo, velxinhi, &
                           velyin, velyinlo, velyinhi, &
                           velzin, velzinlo, velzinhi, &
                           velxout, velxoutlo, velxouthi, &
                           velyout, velyoutlo, velyouthi, &
                           velzout, velzoutlo, velzouthi, &
                           alphax, alphaxlo, alphaxhi, &
                           alphay, alphaylo, alphayhi, &
                           alphaz, alphazlo, alphazhi, &
                           dx,color) &
                           bind (C,name="stag_apply_op")

    integer         , intent(in   ) :: lo(3), hi(3), betacclo(3), betacchi(3), gammacclo(3), gammacchi(3)
    integer         , intent(in   ) :: betaxylo(3), betaxyhi(3), betayzlo(3), betayzhi(3), betaxzlo(3), betaxzhi(3)
    integer         , intent(in   ) :: alphaxlo(3), alphaxhi(3), alphaylo(3), alphayhi(3), alphazlo(3), alphazhi(3)
    integer         , intent(in   ) :: velxinlo(3), velxinhi(3), velyinlo(3), velyinhi(3), velzinlo(3), velzinhi(3)
    integer         , intent(in   ) :: velxoutlo(3), velxouthi(3), velyoutlo(3), velyouthi(3), velzoutlo(3), velzouthi(3)
    double precision, intent(in   ) :: dx(3)
    double precision, intent(in   ) :: betacc(betacclo(1):betacchi(1),betacclo(2):betacchi(2),betacclo(3):betacchi(3))
    double precision, intent(in   ) :: gammacc(gammacclo(1):gammacchi(1),gammacclo(2):gammacchi(2),gammacclo(3):gammacchi(3))
    double precision, intent(in   ) :: betaxy(betaxylo(1):betaxyhi(1),betaxylo(2):betaxyhi(2),betaxylo(3):betaxyhi(3))
    double precision, intent(in   ) :: betaxz(betaxylo(1):betaxzhi(1),betaxzlo(2):betaxzhi(2),betaxzlo(3):betaxzhi(3))
    double precision, intent(in   ) :: betayz(betayzlo(1):betayzhi(1),betayzlo(2):betayzhi(2),betayzlo(3):betayzhi(3))
    double precision, intent(in   ) :: alphax(alphaxlo(1):alphaxhi(1),alphaxlo(2):alphaxhi(2),alphaxlo(3):alphaxhi(3))
    double precision, intent(in   ) :: alphay(alphaylo(1):alphayhi(1),alphaylo(2):alphayhi(2),alphaylo(3):alphayhi(3))
    double precision, intent(in   ) :: alphaz(alphazlo(1):alphazhi(1),alphazlo(2):alphazhi(2),alphazlo(3):alphazhi(3))

    double precision, intent(in   ) :: velxin(velxinlo(1):velxinhi(1),velxinlo(2):velxinhi(2),velxinlo(3):velxinhi(3))
    double precision, intent(in   ) :: velyin(velyinlo(1):velyinhi(1),velyinlo(2):velyinhi(2),velyinlo(3):velyinhi(3))
    double precision, intent(in   ) :: velzin(velzinlo(1):velzinhi(1),velzinlo(2):velzinhi(2),velzinlo(3):velzinhi(3))
    double precision, intent(inout) :: velxout(velxoutlo(1):velxouthi(1),velxoutlo(2):velxouthi(2),velxoutlo(3):velxouthi(3))
    double precision, intent(inout) :: velyout(velyoutlo(1):velyouthi(1),velyoutlo(2):velyouthi(2),velyoutlo(3):velyouthi(3))
    double precision, intent(inout) :: velzout(velzoutlo(1):velzouthi(1),velzoutlo(2):velzouthi(2),velzoutlo(3):velzouthi(3))
    integer         , intent(in   ) :: color

    ! local
    integer :: i,j,k

    double precision dxsqinv, dysqinv, dzsqinv, dxdyinv, dxdzinv, dydzinv
    double precision bt, gm

    dxsqinv = 1.d0/(dx(1)**2)
    dysqinv = 1.d0/(dx(2)**2)
    dzsqinv = 1.d0/(dx(3)**2)

    dxdyinv = 1.d0/(dx(1)*dx(2))
    dxdzinv = 1.d0/(dx(1)*dx(3))
    dydzinv = 1.d0/(dx(2)*dx(3))

    !Type 1

    if (visc_type .eq. -1) then

      do k = lo(3), hi(3)      
        do j = lo(2), hi(2)
          do i = lo(1), hi(1)+1

            velxout(i,j,k) = velxin(i,j,k)*(alphax(i,j,k) + &
                            (betacc(i,j,k)+betacc(i-1,j,k))*dxsqinv+(betaxy(i,j,k)+betaxy(i,j+1,k))*dysqinv+(betaxz(i,j,k)+betaxz(i,j,k+1))*dzsqinv) &
                            -velxin(i+1,j,k)*betacc(i,j,k)*dxsqinv &
                            -velxin(i-1,j,k)*betacc(i-1,j,k)*dxsqinv &
                            -velxin(i,j+1,k)*betaxy(i,j+1,k)*dysqinv &
                            -velxin(i,j-1,k)*betaxy(i,j,k)*dysqinv &
                            -velxin(i,j,k+1)*betaxz(i,j,k+1)*dzsqinv &
                            -velxin(i,j,k-1)*betaxz(i,j,k)*dysqinv
          enddo
        enddo
      enddo

      do k = lo(3), hi(3)      
        do j = lo(2), hi(2)+1
          do i = lo(1), hi(1)

            velyout(i,j,k) = velyin(i,j,k)*(alphay(i,j,k) + &
                            (betacc(i,j,k)+betacc(i,j-1,k))*dysqinv+(betaxy(i,j,k)+betaxy(i+1,j,k))*dxsqinv+(betayz(i,j,k)+betayz(i,j,k+1))*dzsqinv ) &
                            -velyin(i,j+1,k)*betacc(i,j,k)*dysqinv &
                            -velyin(i,j-1,k)*betacc(i,j-1,k)*dysqinv &
                            -velyin(i+1,j,k)*betaxy(i+1,j,k)*dxsqinv &
                            -velyin(i-1,j,k)*betaxy(i,j,k)*dxsqinv &
                            -velyin(i,j,k+1)*betayz(i,j,k+1)*dzsqinv &
                            -velyin(i,j,k-1)*betayz(i,j,k)*dzsqinv 

          enddo
        enddo
      enddo

      do k = lo(3), hi(3)+1      
        do j = lo(2), hi(2)
          do i = lo(1), hi(1)

            velzout(i,j,k) = velzin(i,j,k)*(alphaz(i,j,k) + &
                            (betacc(i,j,k)+betacc(i,j,k-1))*dzsqinv+(betaxz(i,j,k)+betaxz(i+1,j,k))*dxsqinv+(betayz(i,j,k)+betayz(i,j+1,k))*dysqinv) &
                            -velzin(i,j,k+1)*betacc(i,j,k)*dzsqinv &
                            -velzin(i,j,k-1)*betacc(i,j,k-1)*dzsqinv &
                            -velzin(i+1,j,k)*betaxz(i+1,j,k)*dxsqinv &
                            -velzin(i-1,j,k)*betaxz(i,j,k)*dxsqinv &
                            -velzin(i,j+1,k)*betayz(i,j+1,k)*dysqinv &
                            -velzin(i,j-1,k)*betayz(i,j,k)*dysqinv 
          enddo
        enddo
      enddo

    endif

    !Type 2

    if (visc_type .eq. -2) then

      do k = lo(3), hi(3)      
        do j = lo(2), hi(2)
          do i = lo(1), hi(1)+1


           velxout(i,j,k) = velxin(i,j,k)*( alphax(i,j,k) + & 
                             2.d0*(betacc(i,j,k)+betacc(i-1,j,k))*dxsqinv + (betaxy(i,j,k)+betaxy(i,j+1,k))*dxsqinv+(betaxz(i,j,k)+betaxz(i,j,k+1))*dzsqinv  ) &
                        
                            -2.d0*velxin(i+1,j,k)*betacc(i,j,k)*dxsqinv &
                            -2.d0*velxin(i-1,j,k)*betacc(i-1,j,k)*dxsqinv &
                            -velxin(i,j+1,k)*betaxy(i,j+1,k)*dysqinv &
                            -velxin(i,j-1,k)*betaxy(i,j,k)*dysqinv &
                            -velxin(i,j,k+1)*betaxz(i,j,k+1)*dzsqinv &
                            -velxin(i,j,k-1)*betaxz(i,j,k)*dzsqinv &
                        
                            -velyin(i,j+1,k)*betaxy(i,j+1,k)*dxdyinv &
                            +velyin(i,j,k)*betaxy(i,j,k)*dxdyinv &
                            +velyin(i-1,j+1,k)*betaxy(i,j+1,k)*dxdyinv &
                            -velyin(i-1,j,k)*betaxy(i,j,k)*dxdyinv &
                        
                            -velzin(i,j,k+1)*betaxz(i,j,k+1)*dxdzinv &
                            +velzin(i,j,k)*betaxz(i,j,k)*dxdzinv &
                            +velzin(i-1,j,k+1)*betaxz(i,j,k+1)*dxdzinv &
                            -velzin(i-1,j,k)*betaxz(i,j,k)*dxdzinv
          enddo
        enddo
      enddo

      do k = lo(3), hi(3)      
        do j = lo(2), hi(2)+1
          do i = lo(1), hi(1)

            velyout(i,j,k) = velyin(i,j,k)*( alphay(i,j,k) + &
                             2.d0*(betacc(i,j,k)+betacc(i,j-1,k))*dysqinv + (betaxy(i,j,k)+betaxy(i+1,j,k))*dxsqinv + (betayz(i,j,k)+betayz(i,j,k+1))*dzsqinv ) &
                        
                            -2.d0*velyin(i,j+1,k)*betacc(i,j,k)*dysqinv &
                            -2.d0*velyin(i,j-1,k)*betacc(i,j-1,k)*dysqinv &
                            -velyin(i+1,j,k)*betaxy(i+1,j,k)*dxsqinv &
                            -velyin(i-1,j,k)*betaxy(i,j,k)*dxsqinv &
                            -velyin(i,j,k+1)*betayz(i,j,k+1)*dzsqinv &
                            -velyin(i,j,k-1)*betayz(i,j,k)*dzsqinv &
                            
                            -velxin(i+1,j,k)*betaxy(i+1,j,k)*dxdyinv &
                            +velxin(i,j,k)*betaxy(i,j,k)*dxdyinv &
                            +velxin(i+1,j-1,k)*betaxy(i+1,j,k)*dxdyinv &
                            -velxin(i,j-1,k)*betaxy(i,j,k)*dxdyinv &
                            
                            -velzin(i,j,k+1)*betayz(i,j,k+1)*dydzinv &
                            +velzin(i,j,k)*betayz(i,j,k)*dydzinv &
                            +velzin(i,j-1,k+1)*betayz(i,j,k+1)*dydzinv &
                            -velzin(i,j-1,k)*betayz(i,j,k)*dydzinv

          enddo
        enddo
      enddo

      do k = lo(3), hi(3)+1      
        do j = lo(2), hi(2)
          do i = lo(1), hi(1)

            velzout(i,j,k) = velzin(i,j,k)*( alphaz(i,j,k) + &
                             2.d0*(betacc(i,j,k)+betacc(i,j,k-1))*dzsqinv + (betaxz(i,j,k)+betaxz(i+1,j,k))*dxsqinv + (betayz(i,j,k)+betayz(i,j+1,k))*dysqinv ) &
                        
                            -2.d0*velzin(i,j,k+1)*betacc(i,j,k)*dzsqinv &
                            -2.d0*velzin(i,j,k-1)*betacc(i,j,k-1)*dzsqinv &
                            -velzin(i+1,j,k)*betaxz(i+1,j,k)*dxsqinv &
                            -velzin(i-1,j,k)*betaxz(i,j,k)*dxsqinv &
                            -velzin(i,j+1,k)*betayz(i,j+1,k)*dysqinv &
                            -velzin(i,j-1,k)*betayz(i,j,k)*dysqinv &
                            
                            -velxin(i+1,j,k)*betaxz(i+1,j,k)*dxdzinv &
                            +velxin(i,j,k)*betaxz(i,j,k)*dxdzinv &
                            +velxin(i+1,j,k-1)*betaxz(i+1,j,k)*dxdzinv &
                            -velxin(i,j,k-1)*betaxz(i,j,k)*dxdzinv &
                            
                            -velyin(i,j+1,k)*betayz(i,j+1,k)*dydzinv &
                            +velyin(i,j,k)*betayz(i,j,k)*dydzinv &
                            +velyin(i,j+1,k-1)*betayz(i,j+1,k)*dydzinv &
                            -velyin(i,j,k-1)*betayz(i,j,k)*dydzinv
          enddo
        enddo
      enddo

    endif

  
  end subroutine stag_apply_op


#endif

end module stag_apply_op_module
