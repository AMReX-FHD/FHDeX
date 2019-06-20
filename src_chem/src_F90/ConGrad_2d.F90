subroutine get_congrad_2d( lo, hi, &
     &            con , p_lo, p_hi, &
     &            con_x, px_lo, px_hi, &
     &            con_y, py_lo, py_hi,&
     &            iface, if_lo, if_hi,&
     &            ib_cen_x, ib_cen_y, &
     &            dx, prob_lo) bind(C, name="get_congrad_2d")
  

  implicit none

  integer, intent(in) :: lo(2), hi(2)
  double precision, intent(in) :: dx(2), prob_lo(2), ib_cen_x, ib_cen_y
  integer, intent(in) :: p_lo(2), p_hi(2)
  integer, intent(in) :: px_lo(2), px_hi(2)
  integer, intent(in) :: py_lo(2), py_hi(2)
  integer, intent(in) :: if_lo(2), if_hi(2)
  double precision, intent(in   ) :: con (p_lo(1):p_hi(1),p_lo(2):p_hi(2))
  double precision, intent(out) :: con_x(px_lo(1):px_hi(1),px_lo(2):px_hi(2))
  double precision, intent(out) :: con_y(py_lo(1):py_hi(1),py_lo(2):py_hi(2))
  integer, intent(in) :: iface(if_lo(1):if_hi(1),if_lo(2):if_hi(2))

  integer :: i, j
  double precision :: con_xp, con_xm, con_xc, con_yp, con_ym, con_yc, x, y, theta
  double precision ::  n1, n2, t_t1, t_t2

     do    j = lo(2), hi(2)
        do i = lo(1), hi(1)


           y = prob_lo(2) + (dble(j)+0.5d0) * dx(2)-ib_cen_y
           x = prob_lo(1) + (dble(i)+0.5d0) * dx(1)-ib_cen_x



           if (iface(i,j) .lt. 2) then

           ! one sided and centered difference of con
           con_xp=(con(i+1,j)-con(i,j))/dx(1)
           con_xm=(con(i,j)-con(i-1,j))/dx(1)
           con_xc=(con_xp+con_xm)/2

           con_yp=(con(i,j+1)-con(i,j))/dx(2)
           con_ym=(con(i,j)-con(i,j-1))/dx(2)
           con_yc=(con_yp+con_ym)/2


           ! if one side is interior to IB then use other one sided derivative  
              if (iface(i,j) .eq. 1) then
                 if (iface(i+1,j) .eq. 2)then
                 con_x(i,j)=con_xm 
                 else if (iface(i-1,j) .eq. 2) then
                 con_x(i,j)=con_xp
                 else
                 con_x(i,j)=con_xc
                 end if
                 
                 if (iface(i,j+1) .eq. 2)then
                 con_y(i,j)=con_ym 
                 else if (iface(i,j-1) .eq. 2) then
                 con_y(i,j)=con_yp
                 else 
                 con_y(i,j)=con_yc
                 end if

             else 
             con_x(i,j)=con_xc
             con_y(i,j)=con_yc
             end if
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


              theta=atan2(y,x)

              n1=cos(theta)
              n2=sin(theta)


              t_t1=-sin(theta)
              t_t2=cos(theta)


             ! con_n(i,j,k)= con_x !con_x*n1+con_y*n2+con_z*n3
             ! con_t(i,j,k)= con_y !con_x*t_t1+con_y*t_t2+con_z*t_t3
             ! con_p(i,j,k)= con_z !con_x*t_p1+con_y*t_p2+con_z*t_p3 

            !  end if
           else
           con_x(i,j)=0.d0
           con_y(i,j)=0.d0
          endif
        enddo
     enddo

end subroutine get_congrad_2d
