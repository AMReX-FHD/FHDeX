subroutine get_congrad_3d( lo, hi, &
     &            con , p_lo, p_hi, &
     &            con_x, px_lo, px_hi, &
     &            con_y, py_lo, py_hi,&
     &            con_z, pz_lo, pz_hi,&
     &            iface, if_lo, if_hi,&
     &            ib_cen_x, ib_cen_y, ib_cen_z,&
     &            dx, prob_lo) bind(C, name="get_congrad_3d")
  

  use amrex_eb_levelset_module, only : amrex_eb_normal_levelset
  implicit none
 

  integer, intent(in) :: lo(3), hi(3)
  double precision, intent(in) :: dx(3), prob_lo(3), ib_cen_x, ib_cen_y, ib_cen_z
  integer, intent(in) :: p_lo(3), p_hi(3)
  integer, intent(in) :: px_lo(3), px_hi(3)
  integer, intent(in) :: py_lo(3), py_hi(3)
  integer, intent(in) :: pz_lo(3), pz_hi(3)
  integer, intent(in) :: if_lo(3), if_hi(3)

  double precision, intent(in   ) :: con (p_lo(1):p_hi(1),p_lo(2):p_hi(2),p_lo(3):p_hi(3))
  double precision, intent(out) :: con_x(px_lo(1):px_hi(1),px_lo(2):px_hi(2),px_lo(3):px_hi(3))
  double precision, intent(out) :: con_y(py_lo(1):py_hi(1),py_lo(2):py_hi(2),py_lo(3):py_hi(3))
  double precision, intent(out) :: con_z(pz_lo(1):pz_hi(1),pz_lo(2):pz_hi(2),pz_lo(3):pz_hi(3))

  integer, intent(in) :: iface(if_lo(1):if_hi(1),if_lo(2):if_hi(2),if_lo(3):if_hi(3))

  double precision  :: norm(3)! (lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))

  integer :: i, j, k
  double precision :: x, y, z

  ! Do a conservative update
  do       k = lo(3), hi(3)
     do    j = lo(2), hi(2)
        do i = lo(1), hi(1)


           z = prob_lo(3) + (dble(k)+0.5d0) * dx(3)-ib_cen_z
           y = prob_lo(2) + (dble(j)+0.5d0) * dx(2)-ib_cen_y
           x = prob_lo(1) + (dble(i)+0.5d0) * dx(1)-ib_cen_x



!           if (iface(i,j,k) .lt. 2) then
           if (iface(i,j,k) .eq. 0) then
             con_x(i,j,k) =0.
             con_y(i,j,k) =0.
             con_z(i,j,k) =0.
           else
           if ((iface(i,j,k) .eq. 0) .or. (iface(i-1,j,k).eq. 0))then
           con_x(i,j,k) =0.0
           else
           con_x(i,j,k)=(con(i,j,k)-con(i-1,j,k))/dx(1)
           end if
           
            if ((iface(i,j,k) .eq. 0) .or. (iface(i,j-1,k).eq. 0))then
           con_y(i,j,k) =0.0
           else
           con_y(i,j,k)=(con(i,j,k)-con(i,j-1,k))/dx(2)
           end if

           if ((iface(i,j,k) .eq. 0) .or. (iface(i,j,k-1).eq. 0))then
           con_z(i,j,k) =0.0
           else
           con_z(i,j,k)=(con(i,j,k)-con(i,j,k-1))/dx(3)
           end if

           end if     
        enddo
     enddo
  enddo

end subroutine get_congrad_3d
