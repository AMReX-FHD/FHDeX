subroutine get_congrad_3d( lo, hi, &
     &            con , p_lo, p_hi, &
     &            cons_x, px_lo, px_hi, &
     &            cons_y, py_lo, py_hi,&
     &            cons_z, pz_lo, pz_hi,&
     &            magDc, magDc_lo, magDc_hi,&
     &            iface, if_lo, if_hi,&
     &            ls, ls_lo, ls_hi,&
     &            ib_cen_x, ib_cen_y, ib_cen_z,&
     &            dx, prob_lo) bind(C, name="get_congrad_3d")
  

  implicit none

  integer, intent(in) :: lo(3), hi(3)
  double precision, intent(in) :: dx(3), prob_lo(3), ib_cen_x, ib_cen_y, ib_cen_z
  integer, intent(in) :: p_lo(3), p_hi(3)
  integer, intent(in) :: px_lo(3), px_hi(3)
  integer, intent(in) :: py_lo(3), py_hi(3)
  integer, intent(in) :: pz_lo(3), pz_hi(3)
  integer, intent(in) :: magDc_lo(3), magDc_hi(3)
  integer, intent(in) :: if_lo(3), if_hi(3)
  integer, intent(in) :: ls_lo(3), ls_hi(3)

  double precision, intent(in   ) :: con (p_lo(1):p_hi(1),p_lo(2):p_hi(2),p_lo(3):p_hi(3))
  double precision, intent(out) :: cons_x(px_lo(1):px_hi(1),px_lo(2):px_hi(2),px_lo(3):px_hi(3))
  double precision, intent(out) :: cons_y(py_lo(1):py_hi(1),py_lo(2):py_hi(2),py_lo(3):py_hi(3))
  double precision, intent(out) :: cons_z(pz_lo(1):pz_hi(1),pz_lo(2):pz_hi(2),pz_lo(3):pz_hi(3))
  double precision, intent(out) :: magDc(magDc_lo(1):magDc_hi(1),magDc_lo(2):magDc_hi(2),magDc_lo(3):magDc_hi(3))

  integer, intent(in) :: iface(if_lo(1):if_hi(1),if_lo(2):if_hi(2),if_lo(3):if_hi(3))
  double precision, intent(in) :: ls(ls_lo(1):ls_hi(1),ls_lo(2):ls_hi(2),ls_lo(3):ls_hi(3))

  integer :: i, j, k
  double precision :: con_xp, con_xm, con_xc, con_x,  con_yp, con_ym, con_yc, con_y, con_zp, con_zm, con_zc, con_z, nDcon, ls_x, ls_y, ls_z, x, y, z, psi, theta
  double precision ::  n1, n2, n3, t1, t2, t3, t_t1, t_t2, t_t3, t_p1, t_p2, t_p3

           print *, "ConGrad_3D max ls ", maxval(abs(ls)), " min ls ", minval(abs(ls))

  ! Do a conservative update
  do       k = lo(3), hi(3)
     do    j = lo(2), hi(2)
        do i = lo(1), hi(1)


           z = prob_lo(3) + (dble(k)+0.5d0) * dx(3)-ib_cen_z
           y = prob_lo(2) + (dble(j)+0.5d0) * dx(2)-ib_cen_y
           x = prob_lo(1) + (dble(i)+0.5d0) * dx(1)-ib_cen_x



!           if (iface(i,j,k) .lt. 2) then
           if (iface(i,j,k) .eq. 1) then

           ! one sided and centered difference of con
           con_xp=(con(i+1,j,k)-con(i,j,k))/dx(1)
           con_xm=(con(i,j,k)-con(i-1,j,k))/dx(1)
           con_xc=(con_xp+con_xm)/2

           con_yp=(con(i,j+1,k)-con(i,j,k))/dx(2)
           con_ym=(con(i,j,k)-con(i,j-1,k))/dx(2)
           con_yc=(con_yp+con_ym)/2

           con_zp=(con(i,j,k+1)-con(i,j,k))/dx(3)
           con_zm=(con(i,j,k)-con(i,j,k-1))/dx(3)
           con_zc=(con_zp+con_zm)/2

           ls_x=(ls(i+1,j,k)-ls(i-1,j,k))/(2*dx(1))
           ls_y=(ls(i,j+1,k)-ls(i,j-1,k))/(2*dx(2))
           ls_z=(ls(i,j,k+1)-ls(i,j,k-1))/(2*dx(3))

           
           ! if one side is interior to IB then use other one sided derivative  
!              if (iface(i,j,k) .eq. 1) then
                 if (iface(i+1,j,k) .eq. 2)then
                 con_x=con_xm 
                 else if (iface(i-1,j,k) .eq. 2) then
                 con_x=con_xp
                 else
                 con_x=con_xc
                 end if
                 
                 if (iface(i,j+1,k) .eq. 2)then
                 con_y=con_ym 
                 else if (iface(i,j-1,k) .eq. 2) then
                 con_y=con_yp
                 else
                 con_y=con_yc
                 end if

                 if (iface(i,j,k+1) .eq. 2)then
                 con_z=con_zm
                 else if (iface(i,j,k-1) .eq. 2) then
                 con_z=con_zp
                 else
                 con_z=con_zc
                 end if
!              else
!              con_x(i,j,k)=con_xc !con_x*n1+con_y*n2+con_z*n3
!              con_y(i,j,k)=con_yc !con_x*t_t1+con_y*t_t2+con_z*t_t3
!              con_z(i,j,k)=con_zc
!               end if
!              else
!              n1=x/(sqrt(x**2+y**2+z**2))
!              n2=y/(sqrt(x**2+y**2+z**2))
!              n3=z/(sqrt(x**2+y**2+z**2))
!
!              t1=con_x- n1*(con_x*n1+con_y*n2+con_z*n3)
!              t2=con_y- n2*(con_x*n1+con_y*n2+con_z*n3)
!              t3=con_z- n3*(con_x*n1+con_y*n2+con_z*n3)
!             
!              t_t1=-y/(sqrt(x**2+y**2))
!              t_t2=x/(sqrt(x**2+y**2))
!
!              t_p1=(x/sqrt(x**2+y**2))*(z/sqrt(x**2+y**2+z**2))
!              t_p2=(z/sqrt(x**2+y**2+z**2))*(y/sqrt(x**2+y**2))
!              t_p3=-sqrt(x**2+y**2)/sqrt(x**2+y*2+z**2)
!
!              con_n(i,j,k)=con_x*n1+con_y*n2+con_z*n3
!              con_t(i,j,k)=t1*t_t1+t2*t_t2 
!              con_p(i,j,k)=t1*t_p1+t2*t_p2+t3*t_p3
!
!
!              theta=atan2(y,x)
!              psi =acos(z/sqrt(x**2+y**2+z**2))
!
!              n1=sin(psi)*cos(theta)
!              n2=sin(psi)*sin(theta)
!              n3= cos(psi)
!
!              t_p1=cos(psi)*cos(theta)
!              t_p2=cos(psi)*sin(theta)
!              t_p3=-sin(psi)
!
!              t_t1=-sin(theta)
!              t_t2=cos(theta)
!              t_t3=0.d0
!
!             ! t1=con_x- n1*(con_x*n1+con_y*n2+con_z*n3)
!             ! t2=con_y- n2*(con_x*n1+con_y*n2+con_z*n3)
!             ! t3=con_z- n3*(con_x*n1+con_y*n2+con_z*n3)
!
!              con_n(i,j,k)=con_xc !con_x*n1+con_y*n2+con_z*n3
!              con_t(i,j,k)=con_yc !con_x*t_t1+con_y*t_t2+con_z*t_t3
!              con_p(i,j,k)=con_zc !con_x*t_p1+con_y*t_p2+con_z*t_p3 
!              con_t(i,j,k)=t1*t_t1+t2*t_t2+t3*t_t3 
!              con_p(i,j,k)=t1*t_p1+t2*t_p2+t3*t_p3


            !  end if
            nDcon=ls_x*con_x+ls_y*con_y+ls_z*con_z
            cons_x(i,j,k)=con_x-ls_x*nDcon
            cons_y(i,j,k)=con_y-ls_y*nDcon
            cons_z(i,j,k)=con_z-ls_z*nDcon
            magDc(i,j,k)=sqrt( cons_z(i,j,k)*cons_z(i,j,k)+ cons_y(i,j,k)*cons_y(i,j,k)+cons_x(i,j,k)*cons_x(i,j,k))
          ! print *, "dlsdx = ", ls_x, " dCdx = ", con_x(i,j,k), " dlsdy = ", ls_y,  "dCdy = ", con_y(i,j,k), " dlsdz = ", ls_z, " dCdz = ", con_z(i,j,k)
           !print *, "ls_c = ", ls(i,j,k), " ls_r = ", ls(i-1,j,k), " ls_l = ", ls(i+1,j,k), "dx", dx

           else
           cons_x(i,j,k)=0.d0
           cons_y(i,j,k)=0.d0
           cons_z(i,j,k)=0.d0
           magDc(i,j,k)=0.d0
          endif
        enddo
     enddo
  enddo

end subroutine get_congrad_3d
