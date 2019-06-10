subroutine get_phigrad_3d( lo, hi, &
     &            phi , p_lo, p_hi, &
     &            phi_x, px_lo, px_hi, &
     &            phi_y, py_lo, py_hi,&
     &            phi_z, pz_lo, pz_hi,&
     &            iface, if_lo, if_hi,&
     &            ib_cen_x, ib_cen_y, ib_cen_z,&
     &            dx, prob_lo) bind(C, name="get_phigrad_3d")
  

  implicit none

  integer, intent(in) :: lo(3), hi(3)
  double precision, intent(in) :: dx(3), prob_lo(3), ib_cen_x, ib_cen_y, ib_cen_z
  integer, intent(in) :: p_lo(3), p_hi(3)
  integer, intent(in) :: px_lo(3), px_hi(3)
  integer, intent(in) :: py_lo(3), py_hi(3)
  integer, intent(in) :: pz_lo(3), pz_hi(3)
  integer, intent(in) :: if_lo(3), if_hi(3)
  double precision, intent(in   ) :: phi (p_lo(1):p_hi(1),p_lo(2):p_hi(2),p_lo(3):p_hi(3))
  double precision, intent(out) :: phi_x(px_lo(1):px_hi(1),px_lo(2):px_hi(2),px_lo(3):px_hi(3))
  double precision, intent(out) :: phi_y(py_lo(1):py_hi(1),py_lo(2):py_hi(2),py_lo(3):py_hi(3))
  double precision, intent(out) :: phi_z(pz_lo(1):pz_hi(1),pz_lo(2):pz_hi(2),pz_lo(3):pz_hi(3))
  integer, intent(in) :: iface(if_lo(1):if_hi(1),if_lo(2):if_hi(2),if_lo(3):if_hi(3))

  integer :: i, j, k
  double precision :: phi_xp, phi_xm, phi_xc, phi_yp, phi_ym, phi_yc, phi_zp, phi_zm, phi_zc, x, y, z, psi, theta
  double precision ::  n1, n2, n3, t1, t2, t3, t_t1, t_t2, t_t3, t_p1, t_p2, t_p3

  ! Do a conservative update
  do       k = lo(3), hi(3)
     do    j = lo(2), hi(2)
        do i = lo(1), hi(1)


           z = prob_lo(3) + (dble(k)+0.5d0) * dx(3)-ib_cen_z
           y = prob_lo(2) + (dble(j)+0.5d0) * dx(2)-ib_cen_y
           x = prob_lo(1) + (dble(i)+0.5d0) * dx(1)-ib_cen_x



           if (iface(i,j,k) .lt. 2) then

           ! one sided and centered difference of phi
           phi_xp=(phi(i+1,j,k)-phi(i,j,k))/dx(1)
           phi_xm=(phi(i,j,k)-phi(i-1,j,k))/dx(1)
           phi_xc=(phi_xp+phi_xm)/2

           phi_yp=(phi(i,j+1,k)-phi(i,j,k))/dx(2)
           phi_ym=(phi(i,j,k)-phi(i,j-1,k))/dx(2)
           phi_yc=(phi_yp+phi_ym)/2

           phi_zp=(phi(i,j,k+1)-phi(i,j,k))/dx(3)
           phi_zm=(phi(i,j,k)-phi(i,j,k-1))/dx(3)
           phi_zc=(phi_zp+phi_zm)/2

           ! if one side is interior to IB then use other one sided derivative  
              if (iface(i,j,k) .eq. 1) then
                 if (iface(i+1,j,k) .eq. 2)then
                 phi_x(i,j,k)=phi_xm 
                 else if (iface(i-1,j,k) .eq. 2) then
                 phi_x(i,j,k)=phi_xp
                 else
                 phi_x(i,j,k)=phi_xc
                 end if
                 
                 if (iface(i,j+1,k) .eq. 2)then
                 phi_y(i,j,k)=phi_ym 
                 else if (iface(i,j-1,k) .eq. 2) then
                 phi_y(i,j,k)=phi_yp
                 else
                 phi_y(i,j,k)=phi_yc
                 end if

                 if (iface(i,j,k+1) .eq. 2)then
                 phi_z(i,j,k)=phi_zm
                 else if (iface(i,j,k-1) .eq. 2) then
                 phi_z(i,j,k)=phi_zp
                 else
                 phi_z(i,j,k)=phi_zc
                 end if
              else
              phi_x(i,j,k)=phi_xc !phi_x*n1+phi_y*n2+phi_z*n3
              phi_y(i,j,k)=phi_yc !phi_x*t_t1+phi_y*t_t2+phi_z*t_t3
              phi_z(i,j,k)=phi_zc
               end if
!              else
!              n1=x/(sqrt(x**2+y**2+z**2))
!              n2=y/(sqrt(x**2+y**2+z**2))
!              n3=z/(sqrt(x**2+y**2+z**2))
!
!              t1=phi_x- n1*(phi_x*n1+phi_y*n2+phi_z*n3)
!              t2=phi_y- n2*(phi_x*n1+phi_y*n2+phi_z*n3)
!              t3=phi_z- n3*(phi_x*n1+phi_y*n2+phi_z*n3)
!             
!              t_t1=-y/(sqrt(x**2+y**2))
!              t_t2=x/(sqrt(x**2+y**2))
!
!              t_p1=(x/sqrt(x**2+y**2))*(z/sqrt(x**2+y**2+z**2))
!              t_p2=(z/sqrt(x**2+y**2+z**2))*(y/sqrt(x**2+y**2))
!              t_p3=-sqrt(x**2+y**2)/sqrt(x**2+y*2+z**2)
!
!              phi_n(i,j,k)=phi_x*n1+phi_y*n2+phi_z*n3
!              phi_t(i,j,k)=t1*t_t1+t2*t_t2 
!              phi_p(i,j,k)=t1*t_p1+t2*t_p2+t3*t_p3
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
!             ! t1=phi_x- n1*(phi_x*n1+phi_y*n2+phi_z*n3)
!             ! t2=phi_y- n2*(phi_x*n1+phi_y*n2+phi_z*n3)
!             ! t3=phi_z- n3*(phi_x*n1+phi_y*n2+phi_z*n3)
!
!              phi_n(i,j,k)=phi_xc !phi_x*n1+phi_y*n2+phi_z*n3
!              phi_t(i,j,k)=phi_yc !phi_x*t_t1+phi_y*t_t2+phi_z*t_t3
!              phi_p(i,j,k)=phi_zc !phi_x*t_p1+phi_y*t_p2+phi_z*t_p3 
!              phi_t(i,j,k)=t1*t_t1+t2*t_t2+t3*t_t3 
!              phi_p(i,j,k)=t1*t_p1+t2*t_p2+t3*t_p3


            !  end if
           else
           phi_x(i,j,k)=0.d0
           phi_y(i,j,k)=0.d0
           phi_z(i,j,k)=0.d0
          endif
        enddo
     enddo
  enddo

end subroutine get_phigrad_3d
