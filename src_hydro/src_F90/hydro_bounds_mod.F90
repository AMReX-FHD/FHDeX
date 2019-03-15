module hydro_bounds_module

  use amrex_fort_module,      only : amrex_real
  use common_namelist_module, only : bc_lo, bc_hi, p_lo, p_hi

  implicit none

  private

  public :: set_pressure_bc

contains

#if (AMREX_SPACEDIM == 2)

  pure subroutine set_pressure_bc(lo,         hi,            & ! dim(lo) == dim(hi) == 3
       &                          dom_lo,     dom_hi,        &
       &                          data, d_lo, d_hi, d_ncomp, & ! dim(d_lo) == dim(d_hi) == 3
       &                          ngc)                       &
       &                          bind(C, name="set_pressure_bc")

    integer,          intent(in   ) :: lo(3), hi(3), dom_lo(3), dom_hi(3), &
         &                             d_lo(3), d_hi(3), d_ncomp
    integer, value,   intent(in   ) :: ngc
    real(amrex_real), intent(inout) :: data(d_lo(1):d_hi(1), &
         &                                  d_lo(2):d_hi(2), d_ncomp)

    ! ** loop indices
    integer :: i,j


    !____________________________________________________________________________
    ! Apply BC to X faces

    if (lo(1) .eq. dom_lo(1)) then ! lower bound
       if (bc_lo(1) .eq. -2) then ! pressure inflow

          do j = lo(2)-ngc, hi(2)+ngc
             do i = 1, ngc

                data(lo(1)-i, j, :) = p_lo(1)

             end do
          end do

       end if
    end if

    if (hi(1) .eq. dom_hi(1)) then ! upper bound
       if (bc_hi(1) .eq. -2) then ! pressure inflow

          do j = lo(2)-ngc, hi(2)+ngc
             do i = 1, ngc

                data(hi(1)+i, j, :) = p_hi(1)

             end do
          end do

       end if
    end if


    !____________________________________________________________________________
    ! Apply BC to Y faces

    if (lo(2) .eq. dom_lo(2)) then ! lower bound
       if (bc_lo(2) .eq. -2) then ! pressure inflow

          do j = 1, ngc
             do i = lo(1)-ngc, hi(1)+ngc

                data(i, lo(2)-j, :) = p_lo(2)

             end do
          end do

       end if
    end if

    if (hi(2) .eq. dom_hi(2)) then ! upper bound
       if (bc_hi(2) .eq. -2) then ! pressure inflow

          do j = 1, ngc
             do i = lo(1)-ngc, hi(1)+ngc

                data(i, hi(2)+j, :) = p_hi(2)

             end do
          end do

       end if
    end if

  end subroutine set_pressure_bc

#elif (AMREX_SPACEDIM == 3)

  pure subroutine set_pressure_bc(lo,         hi,            & ! dim(lo) == dim(hi) == 3
       &                          dom_lo,     dom_hi,        &
       &                          data, d_lo, d_hi, d_ncomp, & ! dim(d_lo) == dim(d_hi) == 3
       &                          ngc)                       &
       &                          bind(C, name="set_pressure_bc")

    integer,          intent(in   ) :: lo(3), hi(3), dom_lo(3), dom_hi(3), &
         &                             d_lo(3), d_hi(3), d_ncomp
    integer, value,   intent(in   ) :: ngc
    real(amrex_real), intent(inout) :: data(d_lo(1):d_hi(1), &
         &                                  d_lo(2):d_hi(2), &
         &                                  d_lo(3):d_hi(3), d_ncomp)

    ! ** loop indices
    integer :: i,j,k


    !____________________________________________________________________________
    ! Apply BC to X faces

    if (lo(1) .eq. dom_lo(1)) then ! lower bound
       if (bc_lo(1) .eq. -2) then ! pressure inflow

          do k = lo(3)-ngc, hi(3)+ngc
             do j = lo(2)-ngc, hi(2)+ngc
                do i = 1, ngc

                   data(lo(1)-i, j, k, :) = p_lo(1)

                end do
             end do
          end do

       end if
    end if

    if (hi(1) .eq. dom_hi(1)) then ! upper bound
       if (bc_hi(1) .eq. -2) then ! pressure inflow

          do k = lo(3)-ngc, hi(3)+ngc
             do j = lo(2)-ngc, hi(2)+ngc
                do i = 1, ngc

                   data(hi(1)+i, j, k, :) = p_hi(1)

                end do
             end do
          end do

       end if
    end if


    !____________________________________________________________________________
    ! Apply BC to Y faces

    if (lo(2) .eq. dom_lo(2)) then ! lower bound
       if (bc_lo(2) .eq. -2) then ! pressure inflow

          do k = lo(3)-ngc, hi(3)+ngc
             do j = 1, ngc
                do i = lo(1)-ngc, hi(1)+ngc

                   data(i, lo(2)-j, k, :) = p_lo(2)

                end do
             end do
          end do

       end if
    end if

    if (hi(2) .eq. dom_hi(2)) then ! upper bound
       if (bc_hi(2) .eq. -2) then ! pressure inflow

          do k = lo(3)-ngc, hi(3)+ngc
             do j = 1, ngc
                do i = lo(1)-ngc, hi(1)+ngc

                   data(i, hi(2)+j, k, :) = p_hi(2)

                end do
             end do
          end do

       end if
    end if


    !____________________________________________________________________________
    ! Apply BC to Z faces

    if (lo(3) .eq. dom_lo(3)) then ! lower bound
       if (bc_lo(3) .eq. -2) then ! pressure inflow

          do k = 1, ngc
             do j = lo(2)-ngc, hi(2)+ngc
                do i = lo(1)-ngc, hi(1)+ngc

                   data(i, j, lo(3)-k, :) = p_lo(3)

                end do
             end do
          end do

       end if
    end if

    if (hi(3) .eq. dom_hi(3)) then ! upper bound
       if (bc_hi(3) .eq. -2) then ! pressure inflow

          do k = 1, ngc
             do j = lo(2)-ngc, hi(2)+ngc
                do i = lo(1)-ngc, hi(1)+ngc

                   data(i, j, hi(3)+k, :) = p_hi(3)

                end do
             end do
          end do

       end if
    end if

  end subroutine set_pressure_bc

#endif

end module hydro_bounds_module
