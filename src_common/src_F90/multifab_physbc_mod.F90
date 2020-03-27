module multifab_physbc_module

  use amrex_fort_module,      only : amrex_real
  use common_namelist_module, only : bc_vel_lo, bc_vel_hi, bc_es_lo, bc_es_hi, potential_lo, potential_hi

  implicit none

  private

  public :: fab_physbc

contains


#if (AMREX_SPACEDIM == 1)

subroutine fab_physbc

end subroutine fab_physbc

#endif



#if (AMREX_SPACEDIM == 2)

  pure subroutine fab_physbc(lo,     hi,                  & ! dim(lo) == dim(hi) == 3
       &                     dom_lo, dom_hi,              &
       &                     data,   d_lo, d_hi, d_ncomp, & ! dim(d_lo) == dim(d_hi) == 3
       &                     ngc, dim_fill_ghost)         &
       &                     bind(C, name="fab_physbc")

    integer,          intent(in   ) :: lo(3), hi(3), dom_lo(3), dom_hi(3), &
         &                             d_lo(3), d_hi(3), d_ncomp
    integer,          intent(in   ) :: dim_fill_ghost(2)
    integer, value,   intent(in   ) :: ngc
    real(amrex_real), intent(inout) :: data(d_lo(1):d_hi(1), &
         &                                  d_lo(2):d_hi(2), d_ncomp)

    ! ** loop indices
    integer :: i,j

    ! ** number of ghost cells to fill in each dimension
    integer, dimension(2) :: ngc_eff

    ngc_eff(:) = ngc*dim_fill_ghost(:)


    !____________________________________________________________________________
    ! Apply BC to X faces

    if (lo(1) .eq. dom_lo(1)) then ! lower bound
       if(bc_vel_lo(1) .eq. 2) then ! no slip thermal

          do j = lo(2)-ngc_eff(2), hi(2)+ngc_eff(2)
             do i = 1, ngc ! always fill the ghost cells at the bc face

                data(lo(1)-i, j, :) = data(lo(1)-1+i, j, :)

             end do
          end do

       end if
    end if

    if(hi(1) .eq. dom_hi(1)) then ! upper bound
       if(bc_vel_hi(1) .eq. 2) then ! no slip thermal

          do j = lo(2)-ngc_eff(2), hi(2)+ngc_eff(2)
             do i = 1, ngc ! always fill the ghost cells at the bc face

                data(hi(1)+i, j, :) = data(hi(1)+1-i, j, :)

             end do
          end do

       end if
    end if


    !____________________________________________________________________________
    ! Apply BC to Y faces

    if(lo(2) .eq. dom_lo(2)) then ! lower bound
       if(bc_vel_lo(2) .eq. 2) then ! no slip thermal

          do j = 1, ngc ! always fill the ghost cells at the bc face
             do i = lo(1)-ngc_eff(1), hi(1)+ngc_eff(1)

                data(i, lo(2)-j, :) = data(i, lo(2)-1+j, :)

             end do
          end do

       end if
    end if

    if(hi(2) .eq. dom_hi(2)) then ! upper bound
       if(bc_vel_hi(2) .eq. 2) then ! no slip thermal

          do j = 1, ngc ! always fill the ghost cells at the bc face
             do i = lo(1)-ngc_eff(1), hi(1)+ngc_eff(1)

                data(i, hi(2)+j, :) = data(i, hi(2)+1-j, :)

             end do
          end do

       end if
    end if

  end subroutine fab_physbc

#elif (AMREX_SPACEDIM == 3)

  pure subroutine fab_physbc(lo,     hi,                  & ! dim(lo) == dim(hi) == 3
       &                     dom_lo, dom_hi,              &
       &                     data,   d_lo, d_hi, d_ncomp, & ! dim(d_lo) == dim(d_hi) == 3
       &                     ngc, dim_fill_ghost)         &
       &                     bind(C, name="fab_physbc")

    integer,          intent(in   ) :: lo(3), hi(3), dom_lo(3), dom_hi(3), &
         &                             d_lo(3), d_hi(3), d_ncomp
    integer,          intent(in   ) :: dim_fill_ghost(3)
    integer, value,   intent(in   ) :: ngc
    real(amrex_real), intent(inout) :: data(d_lo(1):d_hi(1), &
         &                                  d_lo(2):d_hi(2), &
         &                                  d_lo(3):d_hi(3), d_ncomp)

    ! ** loop indices
    integer :: i,j,k

    ! ** number of ghost cells to fill in each dimension
    integer, dimension(3) :: ngc_eff

    ngc_eff(:) = ngc*dim_fill_ghost(:)


    !____________________________________________________________________________
    ! Apply BC to X faces

    if (lo(1) .eq. dom_lo(1)) then ! lower bound
       if (bc_vel_lo(1) .eq. 2) then ! no slip thermal

          do k = lo(3)-ngc_eff(3), hi(3)+ngc_eff(3)
             do j = lo(2)-ngc_eff(2), hi(2)+ngc_eff(2)
                do i = 1, ngc ! always fill the ghost cells at the bc face

                   data(lo(1)-i, j, k, :) = data(lo(1)-1+i, j, k, :)

                end do
             end do
          end do

       end if
    end if

    if (hi(1) .eq. dom_hi(1)) then ! upper bound
       if (bc_vel_hi(1) .eq. 2) then ! no slip thermal

          do k = lo(3)-ngc_eff(3), hi(3)+ngc_eff(3)
             do j = lo(2)-ngc_eff(2), hi(2)+ngc_eff(2)
                do i = 1, ngc ! always fill the ghost cells at the bc face

                   data(hi(1)+i, j, k, :) = data(hi(1)+1-i, j, k, :)

                end do
             end do
          end do

       end if
    end if


    !____________________________________________________________________________
    ! Apply BC to Y faces

    if (lo(2) .eq. dom_lo(2)) then ! lower bound
       if (bc_vel_lo(2) .eq. 2) then ! no slip thermal

          do k = lo(3)-ngc_eff(3), hi(3)+ngc_eff(3)
             do j = 1, ngc ! always fill the ghost cells at the bc face
                do i = lo(1)-ngc_eff(1), hi(1)+ngc_eff(1)

                   data(i, lo(2)-j, k, :) = data(i, lo(2)-1+j, k, :)

                end do
             end do
          end do

       end if
    end if

    if (hi(2) .eq. dom_hi(2)) then ! upper bound
       if (bc_vel_hi(2) .eq. 2) then ! no slip thermal

          do k = lo(3)-ngc_eff(3), hi(3)+ngc_eff(3)
             do j = 1, ngc ! always fill the ghost cells at the bc face
                do i = lo(1)-ngc_eff(1), hi(1)+ngc_eff(1)

                   data(i, hi(2)+j, k, :) = data(i, hi(2)+1-j, k, :)

                end do
             end do
          end do

       end if
    end if


    !____________________________________________________________________________
    ! Apply BC to Z faces

    if (lo(3) .eq. dom_lo(3)) then ! lower bound
       if (bc_vel_lo(3) .eq. 2) then ! no slip thermal

          do k = 1, ngc ! always fill the ghost cells at the bc face
             do j = lo(2)-ngc_eff(2), hi(2)+ngc_eff(2)
                do i = lo(1)-ngc_eff(1), hi(1)+ngc_eff(1)

                   data(i, j, lo(3)-k, :) = data(i, j, lo(3)-1+k, :)

                end do
             end do
          end do

       end if
    end if

    if (hi(3) .eq. dom_hi(3)) then ! upper bound
       if (bc_vel_hi(3) .eq. 2) then ! no slip thermal

          do k = 1, ngc ! always fill the ghost cells at the bc face
             do j = lo(2)-ngc_eff(2), hi(2)+ngc_eff(2)
                do i = lo(1)-ngc_eff(1), hi(1)+ngc_eff(1)

                   data(i, j, hi(3)+k, :) = data(i, j, hi(3)+1-k, :)

                end do
             end do
          end do

       end if
    end if

  end subroutine fab_physbc

#endif

end module multifab_physbc_module
