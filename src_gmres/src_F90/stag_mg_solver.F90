module stag_solver_module

  use amrex_error_module
  use common_namelist_module
  use gmres_namelist_module

  implicit none

  private

contains

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
                            beta_xy, w_lo, w_hi, &
#if (AMREX_SPACEDIM == 3)
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
    integer         , intent(in   ) :: w_lo(3), w_hi(3)
    double precision, intent(in   ) :: beta_xy(w_lo(1):w_hi(1),w_lo(2):w_hi(2),w_lo(3):w_hi(3))
#if (AMREX_SPACEDIM == 3)
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

    if (visc_type .eq. 1) then

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
                     (2.d0*beta(i,j,k)+2.d0*beta(i-1,j,k)+beta_xy(i,j,k)+beta_xy(i,j+1,k)) * dxsqinv

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
                     (2.d0*beta(i,j,k)+2.d0*beta(i,j-1,k)+beta_xy(i,j,k)+beta_xy(i+1,j,k)) * dxsqinv

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
                      +beta_xy(i,j,k)+beta_xy(i,j+1,k)) * dxsqinv

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
                      +beta_xy(i,j,k)+beta_xy(i+1,j,k)) * dxsqinv

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

    if (visc_type .eq. 1) then

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
