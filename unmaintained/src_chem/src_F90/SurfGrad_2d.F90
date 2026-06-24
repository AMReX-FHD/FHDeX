subroutine get_surfgrad_2d( lo, hi, &
     &            con , c_lo, c_hi, &
     &            cons_x, csx_lo, csx_hi, &
     &            cons_y, csy_lo, csy_hi,&
     &            MagDconc, m_lo, m_hi, &
     &            ls,    ls_lo, ls_hi, &
     &            iface, if_lo, if_hi,&
     &            dx, prob_lo) bind(C, name="get_surfgrad_2d")


  use amrex_eb_levelset_module, only : amrex_eb_normal_levelset

  implicit none
  ! problem specifications
  double precision, intent(in) :: dx(2), prob_lo(2)
  ! work region
  integer, intent(in) :: lo(2), hi(2)
  integer, intent(in) :: c_lo(2), c_hi(2)
  integer, intent(in) :: m_lo(2),  m_hi(2)

  integer, intent(in) :: if_lo(2), if_hi(2)
  integer, intent(in) :: csx_lo(2), csx_hi(2)
  integer, intent(in) :: csy_lo(2), csy_hi(2)

  integer, intent(in) :: ls_lo(2), ls_hi(2)
  ! ** IN: con   - concentration
  !        ls    - levelset
  !        iface - location of interface on grid

double precision, intent(in   ) :: con(c_lo(1):c_hi(1),c_lo(2):c_hi(2))
  double precision, intent(in) :: ls(ls_lo(1):ls_hi(1),ls_lo(2):ls_hi(2))
 integer, intent(in) :: iface(if_lo(1):if_hi(1),if_lo(2):if_hi(2))

  ! ** OUT: MagDconc - magnitude of the surface gradient of con
  !         cons_x   - x component of the surface gradient of con
  !         cons_y   - y component of the surface gradient of con
  !         cons_z   - z component of the surface gradient of con

 double precision, intent(out) :: MagDconc(m_lo(1):m_hi(1),m_lo(2):m_hi(2))
  double precision, intent(out) :: cons_x(csx_lo(1):csx_hi(1),csx_lo(2):csx_hi(2))
  double precision, intent(out) :: cons_y(csy_lo(1):csy_hi(1),csy_lo(2):csy_hi(2))


  ! norm - normal of level set defined on cell centers
  ! cx, cy, cz  - x, y, and z componets of con gradient

  double precision  :: norm(3), cx, cy
 double precision :: con_xp, con_xm, con_xc, con_yp, con_ym, con_yc, con_x, con_y
  integer :: i, j
  double precision ::  pos(3), x, y, z, s
     do    j = lo(2), hi(2)
        do i = lo(1), hi(1)
         ! cell centered position in grid
           y = prob_lo(2) + (dble(j)+0.5d0) * dx(2)!-ib_cen_y
           x = prob_lo(1) + (dble(i)+0.5d0) * dx(1)!-ib_cen_xi
           z = 0.
          pos=(/x, y, z /)

!           end if
           if (iface(i,j) .lt. 2) then

           ! one sided and centered difference of con
           con_xp=(con(i+1,j)-con(i,j))/(dx(1))
           con_xm=(con(i,j)-con(i-1,j))/(dx(1))
           con_xc=(con_xp+con_xm)/2

           con_yp=(con(i,j+1)-con(i,j))/(dx(2))
           con_ym=(con(i,j)-con(i,j-1))/(dx(2))
           con_yc=(con_yp+con_ym)/2


           ! if one side is interior to IB then use other one sided derivative
           ! here is is implicity decided that as a boundary condition on the interface we are assuming that con_xc =(con(i+1,j,k)-con(i-1,j,k))/(2dx), but if con(i-1,j,k) is interior to the immersed boundary then con(i-1,j,k)=c(i,j,k) )

              if (iface(i,j) .eq. 1) then
                 if (iface(i+1,j) .eq. 2)then
                 con_x=con_xm/2
                 else if (iface(i-1,j) .eq. 2) then
                 con_x=con_xp/2
                 else
                 con_x=con_xc
                 end if

                 if (iface(i,j+1) .eq. 2)then
                 con_y=con_ym/2
                 else if (iface(i,j-1) .eq. 2) then
                 con_y=con_yp/2
                 else
                 con_y=con_yc
                 end if

              else
              con_x=con_xc
              con_y=con_yc
               end if
          else
          ! interior to the IB gradient is zero
              con_x=0.
              con_y=0.
          end if

           ! Currently there is no 2D normal level set function so when the code gets to this point, print out an error and abort
           print *, " We need a 2D normal level set function, aborting "
           call abort()
           ! level set normal defined on cell centers
           call amrex_eb_normal_levelset(pos, prob_lo,  1, &
                                             ls, ls_lo,  ls_hi,     &
                                             dx,  norm )
           ! (grad C dot norm)/(norm dot norm)
           s=(norm(1)*con_x+norm(2)*con_y)/(norm(1)**2+norm(2)**2)
           ! Tangential gradients
           cons_x(i,j)=con_x-norm(1)*s
           cons_y(i,j)=con_y-norm(2)*s
           ! Magnitude of tangential gradients
           MagDconc(i,j)=sqrt(cons_x(i,j)**2+cons_y(i,j)**2)
        enddo
     enddo

end subroutine get_surfgrad_2d
