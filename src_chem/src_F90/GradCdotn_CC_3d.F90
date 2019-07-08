subroutine get_gradCdotncc_3d( lo, hi, &
     &            concx , cx_lo, cx_hi, &
     &            concy , cy_lo, cy_hi, &
     &            concz , cz_lo, cz_hi, &
     &            s, s_lo, s_hi, &
     &            MagDconc, m_lo, m_hi, &
     &            ls, ls_lo, ls_hi,&
     &            dx, prob_lo) bind(C, name="get_gradCdotncc_3d")
  

  use amrex_eb_levelset_module, only : amrex_eb_normal_levelset

  implicit none
 

  integer, intent(in) :: lo(3), hi(3)
  double precision, intent(in) :: dx(3), prob_lo(3)
  integer, intent(in) :: cx_lo(3), cx_hi(3)
  integer, intent(in) :: cy_lo(3), cy_hi(3)
  integer, intent(in) :: cz_lo(3), cz_hi(3)

  integer, intent(in) :: s_lo(3), s_hi(3)
  integer, intent(in) :: m_lo(3),  m_hi(3)
  integer, intent(in) :: ls_lo(3), ls_hi(3)

  double precision, intent(in   ) :: concx(cx_lo(1):cx_hi(1),cx_lo(2):cx_hi(2),cx_lo(3):cx_hi(3))
  double precision, intent(in   ) :: concy(cy_lo(1):cy_hi(1),cy_lo(2):cy_hi(2),cy_lo(3):cy_hi(3))
  double precision, intent(in   ) :: concz(cz_lo(1):cz_hi(1),cz_lo(2):cz_hi(2),cz_lo(3):cz_hi(3))


  double precision, intent(out) ::  s(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3))
  
  double precision, intent(out) :: MagDconc(m_lo(1):m_hi(1),m_lo(2):m_hi(2),m_lo(3):m_hi(3))

  double precision, intent(in) :: ls(ls_lo(1):ls_hi(1),ls_lo(2):ls_hi(2),ls_lo(3):ls_hi(3))

  double precision  :: normc(3), normx(3), normy(3), normz(3)
  integer :: i, j, k
  double precision ::  pos(3), x, y, z
  do       k = lo(3), hi(3)
     do    j = lo(2), hi(2)
        do i = lo(1), hi(1)

           z = prob_lo(3) + (dble(k)+0.5d0) * dx(3)!-ib_cen_z
           y = prob_lo(2) + (dble(j)+0.5d0) * dx(2)!-ib_cen_y
           x = prob_lo(1) + (dble(i)+0.5d0) * dx(1)!-ib_cen_x
           pos=(/x, y, z/)
           call amrex_eb_normal_levelset(pos, prob_lo,  1, &
                                             ls, ls_lo,  ls_hi,     &
                                             dx,  normc ) 
           s(i,j,k)=normc(1)*concx(i,j,k)+normc(2)*concy(i,j,k)+normc(3)*concz(i,j,k)
           MagDconc(i,j,k)=sqrt(concx(i,j,k)**2+concy(i,j,k)**2+concz(i,j,k)**2)
   
        enddo
     enddo
  enddo

end subroutine get_gradCdotncc_3d  
