subroutine get_surfgrad_3d( lo, hi, &
     &            conx , cx_lo, cx_hi, &
     &            cony , cy_lo, cy_hi, &
     &            conz , cz_lo, cz_hi, &
     &            sx   , sx_lo, sx_hi, &
     &            sy   , sy_lo, sy_hi, &
     &            sz   , sz_lo, sz_hi, &
     &            cons_x, csx_lo, csx_hi, &
     &            cons_y, csy_lo, csy_hi,&
     &            cons_z, csz_lo, csz_hi,&
     &            ls,    ls_lo, ls_hi, &
     &            xfc, xfc_lo, xfc_hi, &
     &            yfc, yfc_lo, yfc_hi, &
     &            zfc, zfc_lo, zfc_hi, &
     &            dx, prob_lo) bind(C, name="get_surfgrad_3d")
  

  use amrex_eb_levelset_module, only : amrex_eb_normal_levelset

  implicit none
 

  integer, intent(in) :: lo(3), hi(3)
  double precision, intent(in) :: dx(3), prob_lo(3)
  integer, intent(in) :: cx_lo(3), cx_hi(3)
  integer, intent(in) :: cy_lo(3), cy_hi(3)
  integer, intent(in) :: cz_lo(3), cz_hi(3)
  
  integer, intent(in) :: csx_lo(3), csx_hi(3)
  integer, intent(in) :: csy_lo(3), csy_hi(3)
  integer, intent(in) :: csz_lo(3), csz_hi(3)

  integer, intent(in) :: sx_lo(3), sx_hi(3)
  integer, intent(in) :: sy_lo(3), sy_hi(3)
  integer, intent(in) :: sz_lo(3), sz_hi(3)
  
  integer, intent(in) :: ls_lo(3), ls_hi(3)
  integer, intent(in) :: xfc_lo(3), xfc_hi(3) 
  integer, intent(in) :: yfc_lo(3), yfc_hi(3) 
  integer, intent(in) :: zfc_lo(3), zfc_hi(3)
  
double precision, intent(in   ) :: conx(cx_lo(1):cx_hi(1),cx_lo(2):cx_hi(2),cx_lo(3):cx_hi(3))
  double precision, intent(in   ) :: cony(cy_lo(1):cy_hi(1),cy_lo(2):cy_hi(2),cy_lo(3):cy_hi(3))
  double precision, intent(in   ) :: conz(cz_lo(1):cz_hi(1),cz_lo(2):cz_hi(2),cz_lo(3):cz_hi(3))

  double precision, intent(in   ) :: sx(sx_lo(1):sx_hi(1),sx_lo(2):sx_hi(2),sx_lo(3):sx_hi(3))
  double precision, intent(in   ) :: sy(sy_lo(1):sy_hi(1),sy_lo(2):sy_hi(2),sy_lo(3):sy_hi(3))
  double precision, intent(in   ) :: sz(sz_lo(1):sz_hi(1),sz_lo(2):sz_hi(2),sz_lo(3):sz_hi(3))

  double precision, intent(out) :: cons_x(csx_lo(1):csx_hi(1),csx_lo(2):csx_hi(2),csx_lo(3):csx_hi(3))
  double precision, intent(out) :: cons_y(csy_lo(1):csy_hi(1),csy_lo(2):csy_hi(2),csy_lo(3):csy_hi(3))
  double precision, intent(out) :: cons_z(csz_lo(1):csz_hi(1),csz_lo(2):csz_hi(2),csz_lo(3):csz_hi(3))
  
  double precision, intent(in) :: ls(ls_lo(1):ls_hi(1),ls_lo(2):ls_hi(2),ls_lo(3):ls_hi(3))
  double precision, intent(in) :: xfc(xfc_lo(1):xfc_hi(1),xfc_lo(2):xfc_hi(2),xfc_lo(3):xfc_hi(3),3)
  double precision, intent(in) :: yfc(yfc_lo(1):yfc_hi(1),yfc_lo(2):yfc_hi(2),yfc_lo(3):yfc_hi(3),3)
  double precision, intent(in) :: zfc(zfc_lo(1):zfc_hi(1),zfc_lo(2):zfc_hi(2),zfc_lo(3):zfc_hi(3),3)

  double precision  :: norm(3)
  integer :: i, j, k
  double precision ::  pos(3), x, y, z
  ! Do a conservative update
do       k = lo(3), hi(3)
     do    j = lo(2), hi(2)
        do i = lo(1), hi(1)
           call amrex_eb_normal_levelset(xfc(i,j,k,:), prob_lo,  1, &
                                             ls, ls_lo,  ls_hi,     &
                                             dx,  norm )
          cons_x(i,j,k)=conx(i,j,k)-norm(1)*sx(i,j,k)

           call amrex_eb_normal_levelset(yfc(i,j,k,:), prob_lo,  1, &
                                             ls, ls_lo,  ls_hi,     &
                                             dx,  norm )
           cons_y(i,j,k)=cony(i,j,k)-norm(2)*sy(i,j,k)
 
           call amrex_eb_normal_levelset(zfc(i,j,k,:), prob_lo,  1, &
                                             ls, ls_lo,  ls_hi,     &
                                             dx,  norm )
           cons_z(i,j,k)=conz(i,j,k)-norm(3)*sz(i,j,k)
        enddo
     enddo
  enddo

end subroutine get_surfgrad_3d
