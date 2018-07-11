module norm_inner_prod_module

  use amrex_error_module
  ! use common_namelist_module, only: visc_type

  implicit none

  private

contains

#if (AMREX_SPACEDIM == 2)

  subroutine sum_stag(lo, hi, &
                      m1x, m1xlo, m1xhi, &
                      m1y, m1ylo, m1yhi, &
                      sum         ) &
                      bind (C,name="sum_stag")

    integer         , intent(in   ) :: lo(2), hi(2)
    integer         , intent(in   ) :: m1xlo(2),m1xhi(2)
    integer         , intent(in   ) :: m1ylo(2),m1yhi(2)
    double precision, intent(in   ) :: m1x(m1xlo(1):m1xhi(1),m1xlo(2):m1xhi(2))
    double precision, intent(in   ) :: m1y(m1ylo(1):m1yhi(1),m1ylo(2):m1yhi(2))
    double precision, intent(inout) :: sum(2)

    ! local
    integer :: i,j

    ! Inner product of x-component
    ! Loop over interior values & cell-centered boundaries
    do j = lo(2), hi(2)
       do i = lo(1)+1, hi(1)
          sum(1) = sum(1) + m1x(i,j)
       enddo
    enddo
    ! Average nodal boundaries
    do j = lo(2), hi(2)
       sum(1) = sum(1) + 0.5d0*( & 
            m1x(lo(1)  ,j) + &
            m1x(hi(1)+1,j))
    enddo

    ! Inner product of y-component
    ! Loop over interior values & cell-centered boundaries
    do j = lo(2)+1, hi(2)
       do i = lo(1), hi(1)
          sum(2) = sum(2) + m1y(i,j)
       enddo
    enddo
    ! Average nodal boundaries
    do i = lo(1), hi(1)
       sum(2) = sum(2) + 0.5d0*( & 
            m1y(i,lo(2)  ) + &
            m1y(i,hi(2)+1))
    enddo
    
  end subroutine sum_stag


#endif

#if (AMREX_SPACEDIM == 3)

  subroutine sum_stag(lo, hi, &
                      m1x, m1xlo, m1xhi, &
                      m1y, m1ylo, m1yhi, &
                      m1z, m1zlo, m1zhi, &
                      sum         ) &
                      bind (C,name="sum_stag")

    integer         , intent(in   ) :: lo(3), hi(3)
    integer         , intent(in   ) :: m1xlo(3),m1xhi(3)
    integer         , intent(in   ) :: m1ylo(3),m1yhi(3)
    integer         , intent(in   ) :: m1zlo(3),m1zhi(3)
    double precision, intent(in   ) :: m1x(m1xlo(1):m1xhi(1),m1xlo(2):m1xhi(2),m1xlo(3):m1xhi(3))
    double precision, intent(in   ) :: m1y(m1ylo(1):m1yhi(1),m1ylo(2):m1yhi(2),m1ylo(3):m1yhi(3)) 
    double precision, intent(in   ) :: m1z(m1zlo(1):m1zhi(1),m1zlo(2):m1zhi(2),m1zlo(3):m1zhi(3)) 
    double precision, intent(inout) :: sum(3)

    ! local
    integer :: i,j,k

    ! Inner product of x-component
    ! Loop over interior values & cell-centered boundaries
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1)+1, hi(1)
             sum(1) = sum(1) + m1x(i,j,k)
          enddo
       enddo
    enddo
    ! Average nodal boundaries
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          sum(1) = sum(1) + 0.5d0*( & 
               m1x(lo(1)  ,j,k) + &
               m1x(hi(1)+1,j,k))
       enddo
    enddo

    ! Inner product of y-component
    ! Loop over interior values & cell-centered boundaries
    do k = lo(3), hi(3)
       do j = lo(2)+1, hi(2)
          do i = lo(1), hi(1)
             sum(2) = sum(2) + m1y(i,j,k)
          enddo
       enddo
    enddo
    ! Average nodal boundaries
    do k = lo(3), hi(3)
       do i = lo(1), hi(1)
          sum(2) = sum(2) + 0.5d0*( & 
               m1y(i,lo(2)  ,k) + &
               m1y(i,hi(2)+1,k))
       enddo
    enddo

    ! Inner product of z-component
    ! Loop over interior values & cell-centered boundaries
    do k = lo(3)+1, hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             sum(3) = sum(3) + m1z(i,j,k)
          enddo
       enddo
    enddo
    ! Average nodal boundaries
    do j = lo(2), hi(2)
       do i = lo(1), hi(1)
          sum(3) = sum(3) + 0.5d0*( & 
               m1z(i,j,lo(3)  ) + &
               m1z(i,j,hi(3)+1))
       enddo
    enddo
  
  end subroutine sum_stag

#endif

end module norm_inner_prod_module
