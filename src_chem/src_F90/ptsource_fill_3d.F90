subroutine get_ptsource_3d( lo, hi, &
     &            iface, if_lo, if_hi, &
     &            ptS, pts_lo, pts_hi, &
     &            strength, dx, prob_lo) bind(C, name="get_ptsource_3d")
  
  use amrex_mempool_module, only : bl_allocate, bl_deallocate

  implicit none
  double precision, intent(in) :: strength, prob_lo(3), dx(3)
  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) :: if_lo(3), if_hi(3)
  integer, intent(in) :: pts_lo(3), pts_hi(3)
  integer, intent(in) :: iface(if_lo(1):if_hi(1),if_lo(2):if_hi(2),if_lo(3):if_hi(3))
  double precision, intent(out) :: ptS(pts_lo(1):pts_hi(1),pts_lo(2):pts_hi(2),pts_lo(3):pts_hi(3))

  integer :: i, j, k
  double precision:: y, z


 do       k = lo(3), hi(3)
    do    j = lo(2), hi(2)
       do i = lo(1), hi(1)

           ! z = prob_lo(3) + (dble(k)+0.5d0) * dx(3)
           y = prob_lo(2) + (dble(j)+0.5d0) * dx(2)

           if (iface(i,j,k).eq.1) then
           pts(i,j,k)=strength*dx(1)*dx(2)*dx(3)
           endif
        enddo
     enddo
  enddo
end subroutine get_ptsource_3d
