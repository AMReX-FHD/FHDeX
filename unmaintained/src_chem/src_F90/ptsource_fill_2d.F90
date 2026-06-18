subroutine get_ptsource_2d( lo, hi, &
     &            ctag, c_lo, c_hi, &
     &            ptS, pts_lo, pts_hi, &
     &            strength, dx, prob_lo, Num_loc ) bind(C, name="get_ptsource_2d")

  use amrex_mempool_module, only : bl_allocate, bl_deallocate

  implicit none
  ! * work region
  double precision, intent(in) :: prob_lo(2), dx(2)
  integer, intent(in) :: lo(2), hi(2)
  integer, intent(in) :: c_lo(2), c_hi(2)
  integer, intent(in) :: pts_lo(2), pts_hi(2)
  ! * IN stength - source strength
  !      ctag    - location of catalyst (1 if catlyst present 0 otherwise)
  double precision, intent(in) :: strength, Num_loc
  integer, intent(in) :: ctag(c_lo(1):c_hi(1),c_lo(2):c_hi(2))
  ! * OUT ptS    - created point sources
  double precision, intent(out) :: ptS(pts_lo(1):pts_hi(1),pts_lo(2):pts_hi(2))

  integer :: i, j
  double precision :: y


     do    j = lo(2), hi(2)
        do i = lo(1), hi(1)
           y = prob_lo(2)+(dble(j)+0.5d0)*dx(2)
           ! if catalyst is present then create a point source, strength scaled by grid size so that solution isn't grid dependent
           if (ctag(i,j).eq.1) then
           pts(i,j)=strength/Num_loc
           endif
        enddo
     enddo

end subroutine get_ptsource_2d
