module stag_inner_prod_module

  use amrex_error_module
  ! use common_namelist_module, only: visc_type

  implicit none

  private

contains

#if (AMREX_SPACEDIM == 2)

  subroutine stag_inner_prod(lo, hi, &
                           m1x, m1xlo, m1xhi, &
                           m1y, m1ylo, m1yhi, &
                           m2x, m2xlo, m2xhi, &
                           m2y, m2ylo, m2yhi, &
                           prod_val) &
                           bind (C,name="stag_inner_prod")
    
    integer         , intent(in   ) :: lo(2), hi(2)
    integer         , intent(in   ) :: m1xlo(2),m1xhi(2)
    integer         , intent(in   ) :: m1ylo(2),m1yhi(2)
    integer         , intent(in   ) :: m2xlo(2),m2xhi(2)
    integer         , intent(in   ) :: m2ylo(2),m2yhi(2)
    double precision, intent(in   ) :: m1x(m1xlo(1):m1xhi(1),m1xlo(2):m1xhi(2))
    double precision, intent(in   ) :: m1y(m1ylo(1):m1yhi(1),m1ylo(2):m1yhi(2))
    double precision, intent(in   ) :: m2x(m2xlo(1):m2xhi(1),m2xlo(2):m2xhi(2))
    double precision, intent(in   ) :: m2y(m2ylo(1):m2yhi(1),m2ylo(2):m2yhi(2))
    double precision, intent(inout) :: prod_val(2)

    ! local
    integer :: i,j

    ! Check bounds
    if ((m1xlo(1).ne.m2xlo(1)).or.(m1xlo(2).ne.m2xlo(2)).or.&
        (m1xhi(1).ne.m2xhi(1)).or.(m1xhi(2).ne.m2xhi(2))) then
       print*, "Error: x-component bounds mismatch"
    end if
    if ((m1ylo(1).ne.m2ylo(1)).or.(m1ylo(2).ne.m2ylo(2)).or.&
        (m1yhi(1).ne.m2yhi(1)).or.(m1yhi(2).ne.m2yhi(2))) then
       print*, "Error: y-component bounds mismatch"
    end if

    ! Inner product of x-component
    ! Loop over interior values & cell-centered boundaries
    do j = lo(2), hi(2)
       do i = lo(1)+1, hi(1)
          prod_val(1) = prod_val(1) + m1x(i,j)*m2x(i,j)
       enddo
    enddo
    ! Average nodal boundaries
    do j = lo(2), hi(2)
       prod_val(1) = prod_val(1) + 0.5d0*( & 
            m1x(lo(1)  ,j)*m2x(lo(1)  ,j) + &
            m1x(hi(1)+1,j)*m2x(hi(1)+1,j))
    enddo

    ! Inner product of y-component
    ! Loop over interior values & cell-centered boundaries
    do j = lo(2)+1, hi(2)
       do i = lo(1), hi(1)
          prod_val(2) = prod_val(2) + m1y(i,j)*m2y(i,j)
       enddo
    enddo
    ! Average nodal boundaries
    do i = lo(1), hi(1)
       prod_val(2) = prod_val(2) + 0.5d0*( & 
            m1y(i,lo(2)  )*m2y(i,lo(2)  ) + &
            m1y(i,hi(2)+1)*m2y(i,hi(2)+1) )
    enddo
    
  end subroutine stag_inner_prod


#endif

#if (AMREX_SPACEDIM == 3)

  subroutine stag_inner_prod(lo, hi, &
                           m1x, m1xlo, m1xhi, &
                           m1y, m1ylo, m1yhi, &
                           m1z, m1zlo, m1zhi, &
                           m2x, m2xlo, m2xhi, &
                           m2y, m2ylo, m2yhi, &
                           m2z, m2zlo, m2zhi, &
                           prod_val) &
                           bind (C,name="stag_inner_prod")

    integer         , intent(in   ) :: lo(3), hi(3)
    integer         , intent(in   ) :: m1xlo(3),m1xhi(3)
    integer         , intent(in   ) :: m1ylo(3),m1yhi(3)
    integer         , intent(in   ) :: m1zlo(3),m1zhi(3)
    integer         , intent(in   ) :: m2xlo(3),m2xhi(3)
    integer         , intent(in   ) :: m2ylo(3),m2yhi(3)
    integer         , intent(in   ) :: m2zlo(3),m2zhi(3)
    double precision, intent(in   ) :: m1x(m1xlo(1):m1xhi(1),m1xlo(2):m1xhi(2),m1xlo(3):m1xhi(3))
    double precision, intent(in   ) :: m1y(m1ylo(1):m1yhi(1),m1ylo(2):m1yhi(2),m1ylo(3):m1yhi(3)) 
    double precision, intent(in   ) :: m1z(m1zlo(1):m1zhi(1),m1zlo(2):m1zhi(2),m1zlo(3):m1zhi(3)) 
    double precision, intent(in   ) :: m2x(m2xlo(1):m2xhi(1),m2xlo(2):m2xhi(2),m2xlo(3):m2xhi(3))
    double precision, intent(in   ) :: m2y(m2ylo(1):m2yhi(1),m2ylo(2):m2yhi(2),m2ylo(3):m2yhi(3)) 
    double precision, intent(in   ) :: m2z(m2zlo(1):m2zhi(1),m2zlo(2):m2zhi(2),m2zlo(3):m2zhi(3))
    double precision, intent(inout) :: prod_val(3)

    ! local
    integer :: i,j,k

    ! Check bounds
    if ((m1xlo(1).ne.m2xlo(1)).or.(m1xlo(2).ne.m2xlo(2)).or.(m1xlo(3).ne.m2xlo(3)).or.&
        (m1xhi(1).ne.m2xhi(1)).or.(m1xhi(2).ne.m2xhi(2)).or.(m1xhi(3).ne.m2xhi(3))) then
       print*, "Error: x-component bounds mismatch"
    end if
    if ((m1ylo(1).ne.m2ylo(1)).or.(m1ylo(2).ne.m2ylo(2)).or.(m1ylo(3).ne.m2ylo(3)).or.&
         (m1yhi(1).ne.m2yhi(1)).or.(m1yhi(2).ne.m2yhi(2)).or.(m1yhi(3).ne.m2yhi(3))) then
       print*, "Error: y-component bounds mismatch"
    end if
    if ((m1zlo(1).ne.m2zlo(1)).or.(m1zlo(2).ne.m2zlo(2)).or.(m1zlo(3).ne.m2zlo(3)).or.&
        (m1zhi(1).ne.m2zhi(1)).or.(m1zhi(2).ne.m2zhi(2)).or.(m1zhi(3).ne.m2zhi(3))) then
       print*, "Error: z-component bounds mismatch"
    end if


    ! Inner product of x-component
    ! Loop over interior values & cell-centered boundaries
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1)+1, hi(1)
             prod_val(1) = prod_val(1) + m1x(i,j,k)*m2x(i,j,k)
          enddo
       enddo
    enddo
    ! Average nodal boundaries
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          prod_val(1) = prod_val(1) + 0.5d0*( & 
               m1x(lo(1)  ,j,k)*m2x(lo(1)  ,j,k) + &
               m1x(hi(1)+1,j,k)*m2x(hi(1)+1,j,k))
       enddo
    enddo

    ! Inner product of y-component
    ! Loop over interior values & cell-centered boundaries
    do k = lo(3), hi(3)
       do j = lo(2)+1, hi(2)
          do i = lo(1), hi(1)
             prod_val(2) = prod_val(2) + m1y(i,j,k)*m2y(i,j,k)
          enddo
       enddo
    enddo
    ! Average nodal boundaries
    do k = lo(3), hi(3)
       do i = lo(1), hi(1)
          prod_val(2) = prod_val(2) + 0.5d0*( & 
               m1y(i,lo(2)  ,k)*m2y(i,lo(2)  ,k) + &
               m1y(i,hi(2)+1,k)*m2y(i,hi(2)+1,k))
       enddo
    enddo

    ! Inner product of z-component
    ! Loop over interior values & cell-centered boundaries
    do k = lo(3)+1, hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             prod_val(3) = prod_val(3) + m1z(i,j,k)*m2z(i,j,k)
          enddo
       enddo
    enddo
    ! Average nodal boundaries
    do j = lo(2), hi(2)
       do i = lo(1), hi(1)
          prod_val(3) = prod_val(3) + 0.5d0*( & 
               m1z(i,j,lo(3)  )*m2z(i,j,lo(3)  ) + &
               m1z(i,j,hi(3)+1)*m2z(i,j,hi(3)+1))
       enddo
    enddo

  
  end subroutine stag_inner_prod

#endif

end module stag_inner_prod_module
