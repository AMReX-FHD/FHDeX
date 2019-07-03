subroutine get_surfgrad_3d( lo, hi, &
     &            conx , cx_lo, cx_hi, &
     &            cony , cy_lo, cy_hi, &
     &            conz , cz_lo, cz_hi, &
     &            cons_x, px_lo, px_hi, &
     &            cons_y, py_lo, py_hi,&
     &            cons_z, pz_lo, pz_hi,&
     &            MagDCon, m_lo, m_hi, &
     &            ls, ls_lo, ls_hi,&
     &            xfc, xfc_lo, xfc_hi,&
     &            yfc, yfc_lo, yfc_hi,&
     &            zfc, zfc_lo, zfc_hi,&
     &            ib_cen_x, ib_cen_y, ib_cen_z,&
     &            dx, prob_lo) bind(C, name="get_surfgrad_3d")
  

  use amrex_eb_levelset_module, only : amrex_eb_normal_levelset
  use convert_stag_module,      only : average_face_to_cc, average_cc_to_face

  implicit none
 

  integer, intent(in) :: lo(3), hi(3)
  double precision, intent(in) :: dx(3), prob_lo(3), ib_cen_x, ib_cen_y, ib_cen_z
  integer, intent(in) :: cx_lo(3), cx_hi(3)
  integer, intent(in) :: cy_lo(3), cy_hi(3)
  integer, intent(in) :: cz_lo(3), cz_hi(3)

  integer, intent(in) :: px_lo(3), px_hi(3)
  integer, intent(in) :: py_lo(3), py_hi(3)
  integer, intent(in) :: pz_lo(3), pz_hi(3)
  integer, intent(in) :: m_lo(3),  m_hi(3)
  integer, intent(in) :: ls_lo(3), ls_hi(3)
  integer, intent(in) :: xfc_lo(3), xfc_hi(3)
  integer, intent(in) :: yfc_lo(3), yfc_hi(3)
  integer, intent(in) :: zfc_lo(3), zfc_hi(3)

  double precision, intent(in   ) :: conx(cx_lo(1):cx_hi(1),cx_lo(2):cx_hi(2),cx_lo(3):cx_hi(3))
  double precision, intent(in   ) :: cony(cy_lo(1):cy_hi(1),cy_lo(2):cy_hi(2),cy_lo(3):cy_hi(3))
  double precision, intent(in   ) :: conz(cz_lo(1):cz_hi(1),cz_lo(2):cz_hi(2),cz_lo(3):cz_hi(3))


  double precision, intent(out) :: cons_x(px_lo(1):px_hi(1),px_lo(2):px_hi(2),px_lo(3):px_hi(3))
  double precision, intent(out) :: cons_y(py_lo(1):py_hi(1),py_lo(2):py_hi(2),py_lo(3):py_hi(3))
  double precision, intent(out) :: cons_z(pz_lo(1):pz_hi(1),pz_lo(2):pz_hi(2),pz_lo(3):pz_hi(3))
  
  double precision, intent(out) :: MagDcon(m_lo(1):m_hi(1),m_lo(2):m_hi(2),m_lo(3):m_hi(3))

  double precision, intent(in) :: ls(ls_lo(1):ls_hi(1),ls_lo(2):ls_hi(2),ls_lo(3):ls_hi(3))
  double precision, intent(in) :: xfc(xfc_lo(1):xfc_hi(1),xfc_lo(2):xfc_hi(2),xfc_lo(3):xfc_hi(3),3)
  double precision, intent(in) :: yfc(yfc_lo(1):yfc_hi(1),yfc_lo(2):yfc_hi(2),yfc_lo(3):yfc_hi(3),3)
  double precision, intent(in) :: zfc(zfc_lo(1):zfc_hi(1),zfc_lo(2):zfc_hi(2),zfc_lo(3):zfc_hi(3),3)

  double precision  :: normc(3), normx(3), normy(3), normz(3)
  double precision  :: ndotGC_c(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
  double precision  :: ndotGC_x(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
  double precision  :: ndotGC_y(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
  double precision  :: ndotGC_z(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
  double precision  :: conx_c(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
  double precision  :: cony_c(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
  double precision  :: conz_c(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
  integer :: i, j, k
  double precision ::  pos(3), x, y, z

  call average_face_to_cc( lo, hi, conx, cx_lo, cx_hi, 1, conx_c, lo, hi, 1, 1, 1 ,1, 0)
  call average_face_to_cc( lo, hi, cony, cy_lo, cy_hi, 1, cony_c, lo, hi, 1, 1, 1 ,1, 1)
  call average_face_to_cc( lo, hi, conz, cz_lo, cz_hi, 1, conz_c, lo, hi, 1, 1, 1 ,1, 2)
  ! Do a conservative update
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
           ndotGC_c(i,j,k)=normc(1)*conx_c(i,j,k)+normc(2)*cony_c(i,j,k)+normc(3)*conz_c(i,j,k)
           MagDcon(i,j,k)=sqrt(conx_c(i,j,k)**2+cony_c(i,j,k)**2+conz_c(i,j,k)**2)
   
        enddo
     enddo
  enddo
call average_cc_to_face(lo, hi, ndotGC_c, lo, hi, 1, ndotGC_x, lo, hi, 1,ndotGC_y, lo, hi, 1, ndotGC_z, lo,hi, 1 ,1, 1, 1) 
  do       k = lo(3), hi(3)
     do    j = lo(2), hi(2)
        do i = lo(1), hi(1)

           call amrex_eb_normal_levelset(xfc(i,j,k,:), prob_lo,  1, &
                                             ls, ls_lo,  ls_hi,     &
                                             dx,  normx )
           print *, " norm x ", normx 
      !     cons_x(i,j,k)=conx(i,j,k)-normx(1)*ndotGC_x

           call amrex_eb_normal_levelset(yfc(i,j,k,:), prob_lo,  1, &
                                             ls, ls_lo,  ls_hi,     &
                                             dx,  normy ) 
      !     cons_y(i,j,k)=cony(i,j,k)-normy(2)*ndotGC_y

           call amrex_eb_normal_levelset(zfc(i,j,k,:), prob_lo,  1, &
                                             ls, ls_lo,  ls_hi,     &
                                             dx,  normz ) 
       !    cons_z(i,j,k)=conz(i,j,k)-normz(3)*ndotGC_z
            cons_x(i,j,k)=0.
            cons_y(i,j,k)=0.
            cons_z(i,j,k)=0.
       
        enddo
     enddo
  enddo

end subroutine get_surfgrad_3d
