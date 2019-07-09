subroutine get_surfgrad_3d( lo, hi, &
     &            con , c_lo, c_hi, &
     &            cons_x, csx_lo, csx_hi, &
     &            cons_y, csy_lo, csy_hi,&
     &            cons_z, csz_lo, csz_hi,&
     &            MagDconc, m_lo, m_hi, & 
     &            ls,    ls_lo, ls_hi, &
     &            iface, if_lo, if_hi,&
     &            dx, prob_lo) bind(C, name="get_surfgrad_3d")
  

  use amrex_eb_levelset_module, only : amrex_eb_normal_levelset

  implicit none
 

  integer, intent(in) :: lo(3), hi(3)
  double precision, intent(in) :: dx(3), prob_lo(3)
  integer, intent(in) :: c_lo(3), c_hi(3)
  integer, intent(in) :: m_lo(3),  m_hi(3)

  integer, intent(in) :: if_lo(3), if_hi(3) 
  integer, intent(in) :: csx_lo(3), csx_hi(3)
  integer, intent(in) :: csy_lo(3), csy_hi(3)
  integer, intent(in) :: csz_lo(3), csz_hi(3)

  
  integer, intent(in) :: ls_lo(3), ls_hi(3)
  
double precision, intent(in   ) :: con(c_lo(1):c_hi(1),c_lo(2):c_hi(2),c_lo(3):c_hi(3))
  double precision, intent(out) :: cons_x(csx_lo(1):csx_hi(1),csx_lo(2):csx_hi(2),csx_lo(3):csx_hi(3))
  double precision, intent(out) :: cons_y(csy_lo(1):csy_hi(1),csy_lo(2):csy_hi(2),csy_lo(3):csy_hi(3))
  double precision, intent(out) :: cons_z(csz_lo(1):csz_hi(1),csz_lo(2):csz_hi(2),csz_lo(3):csz_hi(3))
  
  double precision, intent(in) :: ls(ls_lo(1):ls_hi(1),ls_lo(2):ls_hi(2),ls_lo(3):ls_hi(3))

 double precision, intent(out) :: MagDconc(m_lo(1):m_hi(1),m_lo(2):m_hi(2),m_lo(3):m_hi(3))

 integer, intent(in) :: iface(if_lo(1):if_hi(1),if_lo(2):if_hi(2),if_lo(3):if_hi(3))

  double precision  :: norm(3), cx, cy, cz
  integer :: i, j, k
  double precision ::  pos(3), x, y, z, s
  ! Do a conservative update
do       k = lo(3), hi(3)
     do    j = lo(2), hi(2)
        do i = lo(1), hi(1)
           
           z = prob_lo(3) + (dble(k)+0.5d0) * dx(3)!-ib_cen_z
           y = prob_lo(2) + (dble(j)+0.5d0) * dx(2)!-ib_cen_y
           x = prob_lo(1) + (dble(i)+0.5d0) * dx(1)!-ib_cen_x
          pos=(/x, y, z/)
  
           if ((iface(i,j,k) .eq. 2) .or. (iface(i,j,k).eq. 1))then
           cx =0.0
           else
           cx=(con(i+1,j,k)-con(i-1,j,k))/(2*dx(1))
           end if

            if ((iface(i,j,k) .eq. 2) .or. (iface(i,j,k).eq. 1))then
           cy =0.0
           else
           cy=(con(i,j+1,k)-con(i,j-1,k))/(2*dx(2))
           end if

           if ((iface(i,j,k) .eq. 2) .or. (iface(i,j,k).eq. 1))then
           cz =0.0
           else
           cz=(con(i,j+1,k)-con(i,j-1,k))/(2*dx(3))
           end if

           call amrex_eb_normal_levelset(pos, prob_lo,  1, &
                                             ls, ls_lo,  ls_hi,     &
                                             dx,  norm )
           s=(norm(1)*cx+norm(2)*cy+norm(3)*cz)/(norm(1)**2+norm(2)**2+norm(3)**2)

           cons_x(i,j,k)=cx-norm(1)*s
           cons_y(i,j,k)=cy-norm(2)*s
           cons_z(i,j,k)=cz-norm(3)*s
           MagDconc(i,j,k)=sqrt(cons_x(i,j,k)**2+cons_y(i,j,k)**2+cons_z(i,j,k)**2)

        enddo
     enddo
  enddo

end subroutine get_surfgrad_3d
