module stag_inner_prod_module

  use amrex_error_module
  ! use common_namelist_module, only: visc_type

  implicit none

  private

contains

#if (AMREX_SPACEDIM == 2)

  subroutine stag_inner_prod(m1x, m1xlo, m1xhi, &
                           m1y, m1ylo, m1yhi, &
                           m2x, m2xlo, m2xhi, &
                           m2y, m2ylo, m2yhi, &
                           prod_val) &
                           bind (C,name="stag_inner_prod")

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
    do j = m1xlo(2), m1xhi(2)
       do i = m1xlo(1)+1, m1xhi(1)-1
          ! print*, "Hack: prod_val = ", prod_val
          prod_val(1) = prod_val(1) + m1x(i,j)*m2x(i,j)
       enddo
    enddo
    ! Average nodal boundaries
    do j = m1xlo(2), m1xhi(2)
       ! print*, "Hack: prod_val = ", prod_val
       prod_val(1) = prod_val(1) + 0.5d0*( & 
            m1x(m1xlo(1),j)*m2x(m1xlo(1),j) + &
            m1x(m1xhi(1),j)*m2x(m1xhi(1),j))
    enddo

    ! Inner product of y-component
    ! Loop over interior values & cell-centered boundaries
    do j = m1ylo(2)+1, m1yhi(2)-1
       do i = m1ylo(1), m1yhi(1)
          ! print*, "Hack: prod_val = ", prod_val
          prod_val(2) = prod_val(2) + m1y(i,j)*m2y(i,j)
       enddo
    enddo
    ! Average nodal boundaries
    do i = m1ylo(1), m1yhi(1)
       ! print*, "Hack: prod_val = ", prod_val
       prod_val(2) = prod_val(2) + 0.5d0*( & 
            m1y(i,m1ylo(2))*m2y(i,m1ylo(2)) + &
            m1y(i,m1yhi(2))*m2y(i,m1yhi(2)) )
    enddo
    
    print*, "Hack: prod_val = ", prod_val

  end subroutine stag_inner_prod


#endif

#if (AMREX_SPACEDIM == 3)

  subroutine stag_inner_prod() &
                           bind (C,name="stag_inner_prod")


    ! local
    ! integer :: i,j,k

  
  end subroutine stag_inner_prod

#endif

end module stag_inner_prod_module
