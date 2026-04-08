subroutine get_ptsource_3d( lo, hi, &
     &            ctag, c_lo, c_hi, &
     &            ptS, pts_lo, pts_hi, &
     &            strength, dx, prob_lo, Num_loc) bind(C, name="get_ptsource_3d")

  use amrex_mempool_module, only : bl_allocate, bl_deallocate

  implicit none
  ! * work region
  double precision, intent(in) :: prob_lo(3), dx(3)
  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) :: c_lo(3), c_hi(3)
  integer, intent(in) :: pts_lo(3), pts_hi(3)
  ! * IN stength - source strength
  !      ctag    - location of catalyst (1 if catlyst present 0 otherwise)
  double precision, intent(in) :: strength
  integer, intent(in) :: Num_loc
  integer, intent(in) :: ctag(c_lo(1):c_hi(1),c_lo(2):c_hi(2),c_lo(3):c_hi(3))
  ! * OUT ptS    - created point sources
  double precision, intent(out) :: ptS(pts_lo(1):pts_hi(1),pts_lo(2):pts_hi(2),pts_lo(3):pts_hi(3))

  integer :: i, j, k
  double precision:: numloc
      if (Num_loc ==0)then
          numloc =1
       else
          numloc=Num_loc
      end if
 do       k = lo(3), hi(3)
    do    j = lo(2), hi(2)
       do i = lo(1), hi(1)

           ! if catalyst is present then create a point source, strength scaled by grid size so that solution isn't grid dependent
           if (ctag(i,j,k).eq.1) then
           pts(i,j,k)=strength!(dx(1)*dx(2)*dx(3))
          ! print *, pts(i,j,k)
           endif
        enddo
     enddo
  enddo
end subroutine get_ptsource_3d
