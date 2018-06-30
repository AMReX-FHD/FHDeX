subroutine find_nodal_values(lo, hi, fabin, inlo, inhi, fabout, outlo, outhi, xcheck, ycheck, zcheck) bind(C, name="find_nodal_values")

  use amrex_fort_module, only : amrex_real
  implicit none

  integer, intent(in) :: lo(3), hi(3), inlo(3), inhi(3), outlo(3), outhi(3)
  integer, intent(in) :: xcheck, ycheck, zcheck

  real(amrex_real), intent(in) :: fabin(inlo(1):inhi(1),inlo(2):inhi(2),inlo(3):inhi(3))
  real(amrex_real), intent(inout) :: fabout(outlo(1):outhi(1),outlo(2):outhi(2),outlo(3):outhi(3))

  integer i, j
#if (AMREX_SPACEDIM == 3) 
  integer k
#endif

!print *, fabin(1,1,1)

#if (AMREX_SPACEDIM == 3) 
  if (xCheck .eq. 0 .and. yCheck .eq. 0 .and. zCheck .eq. 0)  then
    do k = lo(3), hi(3)+1
      do j = lo(2), hi(2)+1
        do i = lo(1), hi(1)+1

          fabout(i,j,k) = 0.125*(fabin(i,j,k) + fabin(i-1,j,k) + fabin(i,j-1,k) + fabin(i-1,j-1,k) + fabin(i,j,k-1) + fabin(i-1,j,k-1) + fabin(i,j-1,k-1) + fabin(i-1,j-1,k-1))

        enddo
      enddo
    enddo
  endif

#endif

#if (AMREX_SPACEDIM == 3) 
  if (xCheck .eq. 1 .and. yCheck .eq. 0 .and. zCheck .eq. 0)  then
    do k = lo(3), hi(3)+1
      do j = lo(2), hi(2)+1
        do i = lo(1), hi(1)+1

          fabout(i,j,k) = 0.25*(fabin(i,j,k) + fabin(i,j-1,k) + fabin(i,j,k-1) + fabin(i,j-1,k-1)) 

        enddo
      enddo
    enddo
  endif

#endif

#if (AMREX_SPACEDIM == 3) 
  if (xCheck .eq. 0 .and. yCheck .eq. 1 .and. zCheck .eq. 0)  then
    do k = lo(3), hi(3)+1
      do j = lo(2), hi(2)+1
        do i = lo(1), hi(1)+1

          fabout(i,j,k) = 0.25*(fabin(i,j,k) + fabin(i-1,j,k) + fabin(i,j,k-1) + fabin(i-1,j,k-1))

        enddo
      enddo
    enddo
  endif

#endif

#if (AMREX_SPACEDIM == 3) 
  if (xCheck .eq. 0 .and. yCheck .eq. 0 .and. zCheck .eq. 1)  then
    do k = lo(3), hi(3)+1
      do j = lo(2), hi(2)+1
        do i = lo(1), hi(1)+1

          fabout(i,j,k) = 0.25*(fabin(i,j,k) + fabin(i-1,j,k) + fabin(i,j-1,k) + fabin(i-1,j-1,k))

        enddo
      enddo
    enddo
  endif

#endif

#if (AMREX_SPACEDIM == 2) 
  if (xCheck .eq. 0 .and. yCheck .eq. 0)  then
      do j = lo(2), hi(2)+1
        do i = lo(1), hi(1)+1

          fabout(i,j,0) = 0.25*(fabin(i,j,0) + fabin(i,j-1,0) + fabin(i-1,j,0) + fabin(i-1,j-1,0))

        enddo
      enddo
  endif

#endif

#if (AMREX_SPACEDIM == 2) 
  if (xCheck .eq. 1 .and. yCheck .eq. 0)  then  
      do j = lo(2), hi(2)+1
        do i = lo(1), hi(1)+1

          fabout(i,j,0) = 0.5*(fabin(i,j,0) + fabin(i,j-1,0))

        enddo
      enddo
  endif

#endif

#if (AMREX_SPACEDIM == 2) 
  if (xCheck .eq. 0 .and. yCheck .eq. 1)  then 
      do j = lo(2), hi(2)+1
        do i = lo(1), hi(1)+1

          fabout(i,j,0) = 0.5*(fabin(i,j,0) + fabin(i-1,j,0))

        enddo
      enddo
  endif

#endif


end subroutine find_nodal_values


