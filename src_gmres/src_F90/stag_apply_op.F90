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
                           betaxy, betaxylo, betaxyhi, &
                           phix, phixlo, phixhi, &
                           phiy, phiylo, phiyhi, &
                           Lphix, Lphixlo, Lphixhi, &
                           Lphiy, Lphiylo, Lphiyhi, &
                           alphax, alphaxlo, alphaxhi, &
                           alphay, alphaylo, alphayhi, &
                           dx, color) &
                           bind (C,name="stag_apply_op")

    integer         , intent(in   ) :: lo(2), hi(2)
    integer         , intent(in   ) :: betacclo(2), betacchi(2), betaxylo(2), betaxyhi(2)
    integer         , intent(in   ) :: gammacclo(2), gammacchi(2)
    integer         , intent(in   ) :: alphaxlo(2), alphaxhi(2), alphaylo(2), alphayhi(2)
    integer         , intent(in   ) :: phixlo(2), phixhi(2), phiylo(2), phiyhi(2)
    integer         , intent(in   ) :: Lphixlo(2), Lphixhi(2), Lphiylo(2), Lphiyhi(2)
    double precision, intent(in   ) :: betacc(betacclo(1):betacchi(1),betacclo(2):betacchi(2))
    double precision, intent(in   ) :: gammacc(gammacclo(1):gammacchi(1),gammacclo(2):gammacchi(2))
    double precision, intent(in   ) :: betaxy(betaxylo(1):betaxyhi(1),betaxylo(2):betaxyhi(2))
    double precision, intent(in   ) :: alphax(alphaxlo(1):alphaxhi(1),alphaxlo(2):alphaxhi(2))
    double precision, intent(in   ) :: alphay(alphaylo(1):alphayhi(1),alphaylo(2):alphayhi(2))
    double precision, intent(in   ) :: phix(phixlo(1):phixhi(1),phixlo(2):phixhi(2))
    double precision, intent(in   ) :: phiy(phiylo(1):phiyhi(1),phiylo(2):phiyhi(2))
    double precision, intent(inout) :: Lphix(Lphixlo(1):Lphixhi(1),Lphixlo(2):Lphixhi(2))
    double precision, intent(inout) :: Lphiy(Lphiylo(1):Lphiyhi(1),Lphiylo(2):Lphiyhi(2))
    double precision, intent(in   ) :: dx(2)
    integer         , intent(in   ) :: color

    ! local
    integer :: i,j
    double precision dxsqinv, dysqinv, dxdyinv

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

    dxsqinv = 1.d0/(dx(1)**2)
    dysqinv = 1.d0/(dx(2)**2)
    dxdyinv = 1.d0/(dx(1)*dx(2))

    if (visc_type .eq. -2) then

       if (do_x) then

          do j = lo(2), hi(2)
             ioff = 0
             if ( offset .eq. 2 .and. mod(lo(1)+j,2) .ne. mod(color+1,2) ) ioff = 1
             do i=lo(1)+ioff,hi(1)+1,offset

                Lphix(i,j) = phix(i,j)*(alphax(i,j) + &
                     2.d0*(betacc(i,j)+betacc(i-1,j))*dxsqinv&
                     + (betaxy(i,j+1)+betaxy(i,j))*dysqinv) &
                     
                     -2.d0*phix(i+1,j)*betacc(i,j)*dxsqinv &
                     -2.d0*phix(i-1,j)*betacc(i-1,j)*dxsqinv &
                     -phix(i,j+1)*betaxy(i,j+1)*dysqinv &
                     -phix(i,j-1)*betaxy(i,j)*dysqinv &
                     
                     -phiy(i,j+1)*betaxy(i,j+1)*dxdyinv &
                     +phiy(i,j)*betaxy(i,j)*dxdyinv &
                     +phiy(i-1,j+1)*betaxy(i,j+1)*dxdyinv &
                     -phiy(i-1,j)*betaxy(i,j)*dxdyinv 
             enddo
          enddo

       end if

       if (do_y) then

          do j = lo(2), hi(2)+1
             ioff = 0
             if ( offset .eq. 2 .and. mod(lo(1)+j,2) .ne. mod(color+1,2) ) ioff = 1
             do i=lo(1)+ioff,hi(1),offset

                Lphiy(i,j) = phiy(i,j)*(alphay(i,j) + &
                     2.d0*(betacc(i,j)+betacc(i,j-1))*dysqinv &
                     + (betaxy(i+1,j)+betaxy(i,j))*dxsqinv) &
                        
                     -2.d0*phiy(i,j+1)*betacc(i,j)*dysqinv &
                     -2.d0*phiy(i,j-1)*betacc(i,j-1)*dysqinv &
                     -phiy(i+1,j)*betaxy(i+1,j)*dxsqinv &
                     -phiy(i-1,j)*betaxy(i,j)*dxsqinv &
                        
                     -phix(i+1,j)*betaxy(i+1,j)*dxdyinv &
                     +phix(i,j)*betaxy(i,j)*dxdyinv &
                     +phix(i+1,j-1)*betaxy(i+1,j)*dxdyinv &
                     -phix(i,j-1)*betaxy(i,j)*dxdyinv
             enddo
          enddo

       end if

    end if

  end subroutine stag_apply_op


#endif

#if (AMREX_SPACEDIM == 3)

  subroutine stag_apply_op(lo, hi, &
                           betacc, betacclo, betacchi, &
                           gammacc, gammacclo, gammacchi, &
                           betaxy, betaxylo, betaxyhi, &
                           betaxz, betaxzlo, betaxzhi, &
                           betayz, betayzlo, betayzhi, &
                           phix, phixlo, phixhi, &
                           phiy, phiylo, phiyhi, &
                           phiz, phizlo, phizhi, &
                           Lphix, Lphixlo, Lphixhi, &
                           Lphiy, Lphiylo, Lphiyhi, &
                           Lphiz, Lphizlo, Lphizhi, &
                           alphax, alphaxlo, alphaxhi, &
                           alphay, alphaylo, alphayhi, &
                           alphaz, alphazlo, alphazhi, &
                           dx,color) &
                           bind (C,name="stag_apply_op")

    integer         , intent(in   ) :: lo(3), hi(3), betacclo(3), betacchi(3), gammacclo(3), gammacchi(3)
    integer         , intent(in   ) :: betaxylo(3), betaxyhi(3), betayzlo(3), betayzhi(3), betaxzlo(3), betaxzhi(3)
    integer         , intent(in   ) :: alphaxlo(3), alphaxhi(3), alphaylo(3), alphayhi(3), alphazlo(3), alphazhi(3)
    integer         , intent(in   ) :: phixlo(3), phixhi(3), phiylo(3), phiyhi(3), phizlo(3), phizhi(3)
    integer         , intent(in   ) :: Lphixlo(3), Lphixhi(3), Lphiylo(3), Lphiyhi(3), Lphizlo(3), Lphizhi(3)
    double precision, intent(in   ) :: dx(3)
    double precision, intent(in   ) :: betacc(betacclo(1):betacchi(1),betacclo(2):betacchi(2),betacclo(3):betacchi(3))
    double precision, intent(in   ) :: gammacc(gammacclo(1):gammacchi(1),gammacclo(2):gammacchi(2),gammacclo(3):gammacchi(3))
    double precision, intent(in   ) :: betaxy(betaxylo(1):betaxyhi(1),betaxylo(2):betaxyhi(2),betaxylo(3):betaxyhi(3))
    double precision, intent(in   ) :: betaxz(betaxylo(1):betaxzhi(1),betaxzlo(2):betaxzhi(2),betaxzlo(3):betaxzhi(3))
    double precision, intent(in   ) :: betayz(betayzlo(1):betayzhi(1),betayzlo(2):betayzhi(2),betayzlo(3):betayzhi(3))
    double precision, intent(in   ) :: alphax(alphaxlo(1):alphaxhi(1),alphaxlo(2):alphaxhi(2),alphaxlo(3):alphaxhi(3))
    double precision, intent(in   ) :: alphay(alphaylo(1):alphayhi(1),alphaylo(2):alphayhi(2),alphaylo(3):alphayhi(3))
    double precision, intent(in   ) :: alphaz(alphazlo(1):alphazhi(1),alphazlo(2):alphazhi(2),alphazlo(3):alphazhi(3))

    double precision, intent(in   ) :: phix(phixlo(1):phixhi(1),phixlo(2):phixhi(2),phixlo(3):phixhi(3))
    double precision, intent(in   ) :: phiy(phiylo(1):phiyhi(1),phiylo(2):phiyhi(2),phiylo(3):phiyhi(3))
    double precision, intent(in   ) :: phiz(phizlo(1):phizhi(1),phizlo(2):phizhi(2),phizlo(3):phizhi(3))
    double precision, intent(inout) :: Lphix(Lphixlo(1):Lphixhi(1),Lphixlo(2):Lphixhi(2),Lphixlo(3):Lphixhi(3))
    double precision, intent(inout) :: Lphiy(Lphiylo(1):Lphiyhi(1),Lphiylo(2):Lphiyhi(2),Lphiylo(3):Lphiyhi(3))
    double precision, intent(inout) :: Lphiz(Lphizlo(1):Lphizhi(1),Lphizlo(2):Lphizhi(2),Lphizlo(3):Lphizhi(3))
    integer         , intent(in   ) :: color

    ! local
    integer :: i,j,k

    double precision dxsqinv, dysqinv, dzsqinv, dxdyinv, dxdzinv, dydzinv
    
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
    
    dxsqinv = 1.d0/(dx(1)**2)
    dysqinv = 1.d0/(dx(2)**2)
    dzsqinv = 1.d0/(dx(3)**2)

    dxdyinv = 1.d0/(dx(1)*dx(2))
    dxdzinv = 1.d0/(dx(1)*dx(3))
    dydzinv = 1.d0/(dx(2)*dx(3))

    if (visc_type .eq. -2) then

       if (do_x) then

          do k = lo(3), hi(3)
             do j = lo(2), hi(2)
                ioff = 0
                if ( offset .eq. 2 .and. mod(lo(1)+j+k,2) .ne. mod(color+1,2) ) ioff = 1
                do i=lo(1)+ioff,hi(1)+1,offset
                   
                   Lphix(i,j,k) = phix(i,j,k)*( alphax(i,j,k) + &
                        2.d0*(betacc(i,j,k)+betacc(i-1,j,k))*dxsqinv &
                        + (betaxy(i,j,k)+betaxy(i,j+1,k))*dysqinv&
                        + (betaxz(i,j,k)+betaxz(i,j,k+1))*dzsqinv  ) &
                        
                        -2.d0*phix(i+1,j,k)*betacc(i,j,k)*dxsqinv &
                        -2.d0*phix(i-1,j,k)*betacc(i-1,j,k)*dxsqinv &
                        -phix(i,j+1,k)*betaxy(i,j+1,k)*dysqinv &
                        -phix(i,j-1,k)*betaxy(i,j,k)*dysqinv &
                        -phix(i,j,k+1)*betaxz(i,j,k+1)*dzsqinv &
                        -phix(i,j,k-1)*betaxz(i,j,k)*dzsqinv &
                        
                        -phiy(i,j+1,k)*betaxy(i,j+1,k)*dxdyinv &
                        +phiy(i,j,k)*betaxy(i,j,k)*dxdyinv &
                        +phiy(i-1,j+1,k)*betaxy(i,j+1,k)*dxdyinv &
                        -phiy(i-1,j,k)*betaxy(i,j,k)*dxdyinv &
                        
                        -phiz(i,j,k+1)*betaxz(i,j,k+1)*dxdzinv &
                        +phiz(i,j,k)*betaxz(i,j,k)*dxdzinv &
                        +phiz(i-1,j,k+1)*betaxz(i,j,k+1)*dxdzinv &
                        -phiz(i-1,j,k)*betaxz(i,j,k)*dxdzinv
                enddo
             enddo
          enddo

       end if

       if (do_y) then

          do k = lo(3), hi(3)
             do j = lo(2), hi(2)+1
                ioff = 0
                if ( offset .eq. 2 .and. mod(lo(1)+j+k,2) .ne. mod(color+1,2) ) ioff = 1
                do i=lo(1)+ioff,hi(1),offset
                   
                   Lphiy(i,j,k) = phiy(i,j,k)*( alphay(i,j,k) + &
                        2.d0*(betacc(i,j,k)+betacc(i,j-1,k))*dysqinv &
                        + (betaxy(i,j,k)+betaxy(i+1,j,k))*dxsqinv &
                        + (betayz(i,j,k)+betayz(i,j,k+1))*dzsqinv ) &
                        
                        -2.d0*phiy(i,j+1,k)*betacc(i,j,k)*dysqinv &
                        -2.d0*phiy(i,j-1,k)*betacc(i,j-1,k)*dysqinv &
                        -phiy(i+1,j,k)*betaxy(i+1,j,k)*dxsqinv &
                        -phiy(i-1,j,k)*betaxy(i,j,k)*dxsqinv &
                        -phiy(i,j,k+1)*betayz(i,j,k+1)*dzsqinv &
                        -phiy(i,j,k-1)*betayz(i,j,k)*dzsqinv &
                        
                        -phix(i+1,j,k)*betaxy(i+1,j,k)*dxdyinv &
                        +phix(i,j,k)*betaxy(i,j,k)*dxdyinv &
                        +phix(i+1,j-1,k)*betaxy(i+1,j,k)*dxdyinv &
                        -phix(i,j-1,k)*betaxy(i,j,k)*dxdyinv &
                        
                        -phiz(i,j,k+1)*betayz(i,j,k+1)*dydzinv &
                        +phiz(i,j,k)*betayz(i,j,k)*dydzinv &
                        +phiz(i,j-1,k+1)*betayz(i,j,k+1)*dydzinv &
                        -phiz(i,j-1,k)*betayz(i,j,k)*dydzinv
                enddo
             enddo
          enddo

       end if

       if (do_z) then

          do k = lo(3), hi(3)+1
             do j = lo(2), hi(2)
                ioff = 0
                if ( offset .eq. 2 .and. mod(lo(1)+j+k,2) .ne. mod(color+1,2) ) ioff = 1
                do i=lo(1)+ioff,hi(1),offset

                   Lphiz(i,j,k) = phiz(i,j,k)*( alphaz(i,j,k) + &
                        2.d0*(betacc(i,j,k)+betacc(i,j,k-1))*dzsqinv &
                        + (betaxz(i,j,k)+betaxz(i+1,j,k))*dxsqinv &
                        + (betayz(i,j,k)+betayz(i,j+1,k))*dysqinv ) &
                        
                        -2.d0*phiz(i,j,k+1)*betacc(i,j,k)*dzsqinv &
                        -2.d0*phiz(i,j,k-1)*betacc(i,j,k-1)*dzsqinv &
                        -phiz(i+1,j,k)*betaxz(i+1,j,k)*dxsqinv &
                        -phiz(i-1,j,k)*betaxz(i,j,k)*dxsqinv &
                        -phiz(i,j+1,k)*betayz(i,j+1,k)*dysqinv &
                        -phiz(i,j-1,k)*betayz(i,j,k)*dysqinv &
                        
                        -phix(i+1,j,k)*betaxz(i+1,j,k)*dxdzinv &
                        +phix(i,j,k)*betaxz(i,j,k)*dxdzinv &
                        +phix(i+1,j,k-1)*betaxz(i+1,j,k)*dxdzinv &
                        -phix(i,j,k-1)*betaxz(i,j,k)*dxdzinv &
                        
                        -phiy(i,j+1,k)*betayz(i,j+1,k)*dydzinv &
                        +phiy(i,j,k)*betayz(i,j,k)*dydzinv &
                        +phiy(i,j+1,k-1)*betayz(i,j+1,k)*dydzinv &
                        -phiy(i,j,k-1)*betayz(i,j,k)*dydzinv
                enddo
             enddo
          enddo

       end if

    end if

  end subroutine stag_apply_op

#endif

end module stag_apply_op_module
