module multifab_physbc_module

  use amrex_fort_module,      only : amrex_real
  use common_namelist_module, only : bc_vel_lo, bc_vel_hi, bc_es_lo, bc_es_hi, potential_lo, potential_hi

  implicit none

  private

  public :: fab_physbc
  public :: fab_electricbc
  public :: fab_chargebc
  public :: fab_physbc_domainvel
  public :: fab_physbc_macvel
  public :: fab_physbc_domainstress
  public :: fab_physbc_macstress


contains


#if (AMREX_SPACEDIM == 1)

subroutine fab_physbc

end subroutine fab_physbc


subroutine fab_electricbc

end subroutine fab_electricbc

subroutine fab_chargebc

end subroutine fab_chargebc

subroutine fab_physbc_domainvel

end subroutine fab_physbc_domainvel

subroutine fab_physbc_macvel

end subroutine fab_physbc_macvel

subroutine fab_physbc_domainstress

end subroutine fab_physbc_domainstress

subroutine fab_physbc_macstress

end subroutine fab_physbc_macstress

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

  pure subroutine fab_electricbc(lo,     hi,                  & ! dim(lo) == dim(hi) == 3
       &                     dom_lo, dom_hi,              &
       &                     data,   d_lo, d_hi, d_ncomp, & ! dim(d_lo) == dim(d_hi) == 3
       &                     ngc, dim_fill_ghost)         &
       &                     bind(C, name="fab_electricbc")

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
       if(bc_es_lo(1) .eq. 2) then ! Neumann
          do j = lo(2)-ngc_eff(2), hi(2)+ngc_eff(2)
             do i = 1, ngc ! always fill the ghost cells at the bc face

                data(lo(1)-i, j, :) = -data(lo(1)-1+i, j, :)

             end do
          end do
       elseif(bc_es_lo(1) .eq. 1) then ! Dirichlet
          do j = lo(2)-ngc_eff(2), hi(2)+ngc_eff(2)
             do i = 1, ngc ! always fill the ghost cells at the bc face

                data(lo(1)-i, j, :) = data(lo(1)-1+i, j, :)

             end do
          end do
       end if
    end if

    if(hi(1) .eq. dom_hi(1)) then ! upper bound
       if(bc_es_hi(1) .eq. 2) then ! Neumann
          do j = lo(2)-ngc_eff(2), hi(2)+ngc_eff(2)
             do i = 1, ngc ! always fill the ghost cells at the bc face

                data(hi(1)+i, j, :) = -data(hi(1)+1-i, j, :)

             end do
          end do
       elseif(bc_es_hi(1) .eq. 1) then ! Dirichlet

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
       if(bc_es_lo(2) .eq. 2) then ! Neumann

          do j = 1, ngc ! always fill the ghost cells at the bc face
             do i = lo(1)-ngc_eff(1), hi(1)+ngc_eff(1)

                data(i, lo(2)-j, :) = -data(i, lo(2)-1+j, :)

             end do
          end do

       elseif(bc_es_lo(2) .eq. 1) then ! Dirichlet
          do j = 1, ngc ! always fill the ghost cells at the bc face
             do i = lo(1)-ngc_eff(1), hi(1)+ngc_eff(1)

                data(i, lo(2)-j, :) = data(i, lo(2)-1+j, :)

             end do
          end do

       end if
    end if

    if(hi(2) .eq. dom_hi(2)) then ! upper bound
       if(bc_es_hi(2) .eq. 2) then ! Neumann
          do j = 1, ngc ! always fill the ghost cells at the bc face
             do i = lo(1)-ngc_eff(1), hi(1)+ngc_eff(1)

                data(i, hi(2)+j, :) = -data(i, hi(2)+1-j, :)

             end do
          end do
       elseif(bc_es_hi(2) .eq. 1) then ! Direchlet
          do j = 1, ngc ! always fill the ghost cells at the bc face
             do i = lo(1)-ngc_eff(1), hi(1)+ngc_eff(1)

                data(i, hi(2)+j, :) = data(i, hi(2)+1-j, :)

             end do
          end do
       end if
    end if

  end subroutine fab_electricbc

  pure subroutine fab_potentialbc(lo,     hi,                  & ! dim(lo) == dim(hi) == 3
                                  dom_lo, dom_hi,              &
                                  data,   d_lo, d_hi, d_ncomp, & ! dim(d_lo) == dim(d_hi) == 3
                                  ngc, dim_fill_ghost)         &
                                  bind(C, name="fab_potentialbc")

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
       if(bc_es_lo(1) .eq. 2) then ! Neumann
          do j = lo(2)-ngc_eff(2), hi(2)+ngc_eff(2)
             do i = 1, 1 ! always fill the ghost cells at the bc face

                data(lo(1)-i, j, :) = data(lo(1)-1+i, j, :)

             end do
          end do
       elseif(bc_es_lo(1) .eq. 1) then ! Dirichlet
          do j = lo(2)-ngc_eff(2), hi(2)+ngc_eff(2)
             do i = 1, 1 ! always fill the ghost cells at the bc face

                data(lo(1)-i, j, :) = -data(lo(1)-1+i, j, :) + 2*potential_lo(1)

             end do
          end do
       end if
    end if

    if(hi(1) .eq. dom_hi(1)) then ! upper bound
       if(bc_es_hi(1) .eq. 2) then ! Neumann
          do j = lo(2)-ngc_eff(2), hi(2)+ngc_eff(2)
             do i = 1, 1 ! always fill the ghost cells at the bc face

                data(hi(1)+i, j, :) = data(hi(1)+1-i, j, :)

             end do
          end do
       elseif(bc_es_hi(1) .eq. 1) then ! Dirichlet

          do j = lo(2)-ngc_eff(2), hi(2)+ngc_eff(2)
             do i = 1, 1 ! always fill the ghost cells at the bc face

                data(hi(1)+i, j, :) = -data(hi(1)+1-i, j, :) + 2*potential_hi(1)

             end do
          end do
       end if

    end if


    !____________________________________________________________________________
    ! Apply BC to Y faces

    if(lo(2) .eq. dom_lo(2)) then ! lower bound
       if(bc_es_lo(2) .eq. 2) then ! Neumann

          do j = 1, 1 ! always fill the ghost cells at the bc face
             do i = lo(1)-ngc_eff(1), hi(1)+ngc_eff(1)

                data(i, lo(2)-j, :) = data(i, lo(2)-1+j, :)

             end do
          end do

       elseif(bc_es_lo(2) .eq. 1) then ! Dirichlet
          do j = 1, 1 ! always fill the ghost cells at the bc face
             do i = lo(1)-ngc_eff(1), hi(1)+ngc_eff(1)

                data(i, lo(2)-j, :) = -data(i, lo(2)-1+j, :) + 2*potential_lo(2)

             end do
          end do

       end if
    end if

    if(hi(2) .eq. dom_hi(2)) then ! upper bound
       if(bc_es_hi(2) .eq. 2) then ! Neumann
          do j = 1, 1 ! always fill the ghost cells at the bc face
             do i = lo(1)-ngc_eff(1), hi(1)+ngc_eff(1)

                data(i, hi(2)+j, :) = data(i, hi(2)+1-j, :)

             end do
          end do
       elseif(bc_es_hi(2) .eq. 1) then ! Direchlet
          do j = 1, 1 ! always fill the ghost cells at the bc face
             do i = lo(1)-ngc_eff(1), hi(1)+ngc_eff(1)

                data(i, hi(2)+j, :) = -data(i, hi(2)+1-j, :) + 2*potential_hi(2)

             end do
          end do
       end if
    end if

  end subroutine fab_potentialbc

  pure subroutine fab_potentialbc_solver(lo,     hi,                  & ! dim(lo) == dim(hi) == 3
                                         dom_lo, dom_hi,              &
                                         data,   d_lo, d_hi, d_ncomp, & ! dim(d_lo) == dim(d_hi) == 3
                                         ngc, dim_fill_ghost)         &
                                         bind(C, name="fab_potentialbc_solver")

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
       if(bc_es_lo(1) .eq. 2) then ! Neumann
          do j = lo(2)-ngc_eff(2), hi(2)+ngc_eff(2)
             do i = 1, 1 ! always fill the ghost cells at the bc face

                data(lo(1)-i, j, :) = potential_lo(1)

             end do
          end do
       elseif(bc_es_lo(1) .eq. 1) then ! Dirichlet
          do j = lo(2)-ngc_eff(2), hi(2)+ngc_eff(2)
             do i = 1, 1 ! always fill the ghost cells at the bc face

                data(lo(1)-i, j, :) = potential_lo(1)

             end do
          end do
       end if
    end if

    if(hi(1) .eq. dom_hi(1)) then ! upper bound
       if(bc_es_hi(1) .eq. 2) then ! Neumann
          do j = lo(2)-ngc_eff(2), hi(2)+ngc_eff(2)
             do i = 1, 1 ! always fill the ghost cells at the bc face

                data(hi(1)+i, j, :) = potential_hi(1)

             end do
          end do
       elseif(bc_es_hi(1) .eq. 1) then ! Dirichlet

          do j = lo(2)-ngc_eff(2), hi(2)+ngc_eff(2)
             do i = 1, 1 ! always fill the ghost cells at the bc face

                data(hi(1)+i, j, :) = potential_hi(1)

             end do
          end do
       end if

    end if


    !____________________________________________________________________________
    ! Apply BC to Y faces

    if(lo(2) .eq. dom_lo(2)) then ! lower bound
       if(bc_es_lo(2) .eq. 2) then ! Neumann

          do j = 1, 1 ! always fill the ghost cells at the bc face
             do i = lo(1)-ngc_eff(1), hi(1)+ngc_eff(1)

                data(i, lo(2)-j, :) = potential_lo(2)

             end do
          end do

       elseif(bc_es_lo(2) .eq. 1) then ! Dirichlet
          do j = 1, 1 ! always fill the ghost cells at the bc face
             do i = lo(1)-ngc_eff(1), hi(1)+ngc_eff(1)

                data(i, lo(2)-j, :) = potential_lo(2)

             end do
          end do

       end if
    end if

    if(hi(2) .eq. dom_hi(2)) then ! upper bound
       if(bc_es_hi(2) .eq. 2) then ! Neumann
          do j = 1, 1 ! always fill the ghost cells at the bc face
             do i = lo(1)-ngc_eff(1), hi(1)+ngc_eff(1)

                data(i, hi(2)+j, :) = potential_hi(2)

             end do
          end do
       elseif(bc_es_hi(2) .eq. 1) then ! Direchlet
          do j = 1, 1 ! always fill the ghost cells at the bc face
             do i = lo(1)-ngc_eff(1), hi(1)+ngc_eff(1)

                data(i, hi(2)+j, :) = potential_hi(2)

             end do
          end do
       end if
    end if

  end subroutine fab_potentialbc_solver

  pure subroutine fab_chargebc(lo,     hi,                  & ! dim(lo) == dim(hi) == 3
       &                     dom_lo, dom_hi,              &
       &                     data,   d_lo, d_hi, d_ncomp, & ! dim(d_lo) == dim(d_hi) == 3
       &                     ngc, dim_fill_ghost)         &
       &                     bind(C, name="fab_chargebc")

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
       if(bc_es_lo(1) .eq. 2) then ! Neumann
          do j = lo(2)-ngc_eff(2), hi(2)+ngc_eff(2)
             do i = 1, ngc

                data(lo(1)-1+i, j, :) = data(lo(1)-1+i, j, :) + data(lo(1)-i, j, :)

             end do
          end do
       elseif(bc_es_lo(1) .eq. 1) then ! Dirichlet
          do j = lo(2)-ngc_eff(2), hi(2)+ngc_eff(2)
             do i = 1, ngc ! always fill the ghost cells at the bc face

                data(lo(1)-1+i, j, :) = data(lo(1)-1+i, j, :) - data(lo(1)-i, j, :)

             end do
          end do
       end if


    end if

    if(hi(1) .eq. dom_hi(1)) then ! upper bound
       if(bc_es_hi(1) .eq. 2) then ! Neumann
          do j = lo(2)-ngc_eff(2), hi(2)+ngc_eff(2)
             do i = 1, ngc ! always fill the ghost cells at the bc face

                data(hi(1)+1-i, j, :) = data(hi(1)+1-i, j, :) + data(hi(1)+i, j, :)

             end do
          end do
       elseif(bc_es_hi(1) .eq. 1) then ! Dirichlet

          do j = lo(2)-ngc_eff(2), hi(2)+ngc_eff(2)
             do i = 1, ngc ! always fill the ghost cells at the bc face

                data(hi(1)+1-i, j, :) = data(hi(1)+1-i, j, :) - data(hi(1)+i, j, :)

             end do
          end do
       end if

    end if


    !____________________________________________________________________________
    ! Apply BC to Y faces

    if(lo(2) .eq. dom_lo(2)) then ! lower bound
       if(bc_es_lo(2) .eq. 2) then ! Neumann

          do j = 1, ngc ! always fill the ghost cells at the bc face
             do i = lo(1)-ngc_eff(1), hi(1)+ngc_eff(1)

                data(i, lo(2)-1+j, :) = data(i, lo(2)-1+j, :) + data(i, lo(2)-j, :)

             end do
          end do

       elseif(bc_es_lo(2) .eq. 1) then ! Dirichlet
          do j = 1, ngc ! always fill the ghost cells at the bc face
             do i = lo(1)-ngc_eff(1), hi(1)+ngc_eff(1)

                data(i, lo(2)-1+j, :) = data(i, lo(2)-1+j, :) - data(i, lo(2)-j, :)

             end do
          end do

       end if
    end if

    if(hi(2) .eq. dom_hi(2)) then ! upper bound
       if(bc_es_hi(2) .eq. 2) then ! Neumann
          do j = 1, ngc ! always fill the ghost cells at the bc face
             do i = lo(1)-ngc_eff(1), hi(1)+ngc_eff(1)

                data(i, hi(2)+1-j, :) = -data(i, hi(2)+1-j, :) + data(i, hi(2)+j, :)

             end do
          end do
       elseif(bc_es_hi(2) .eq. 1) then ! Dirichlet
          do j = 1, ngc ! always fill the ghost cells at the bc face
             do i = lo(1)-ngc_eff(1), hi(1)+ngc_eff(1)

                data(i, hi(2)+1-j, :) = data(i, hi(2)+1-j, :) - data(i, hi(2)+j, :)

             end do
          end do
       end if
    end if

  end subroutine fab_chargebc

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

  pure subroutine fab_electricbc(lo,     hi,                  & ! dim(lo) == dim(hi) == 3
       &                     dom_lo, dom_hi,              &
       &                     data,   d_lo, d_hi, d_ncomp, & ! dim(d_lo) == dim(d_hi) == 3
       &                     ngc, dim_fill_ghost)         &
       &                     bind(C, name="fab_electricbc")

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
       if (bc_es_lo(1) .eq. 2) then ! Neumann
          do k = lo(3)-ngc_eff(3), hi(3)+ngc_eff(3)
             do j = lo(2)-ngc_eff(2), hi(2)+ngc_eff(2)
                do i = 1, ngc ! always fill the ghost cells at the bc face

                   data(lo(1)-i, j, k, :) = -data(lo(1)-1+i, j, k, :)

                end do
             end do
          end do
       elseif (bc_es_lo(1) .eq. 1) then !Dirichlet
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
       if (bc_es_hi(1) .eq. 2) then !Neumann
          do k = lo(3)-ngc_eff(3), hi(3)+ngc_eff(3)
             do j = lo(2)-ngc_eff(2), hi(2)+ngc_eff(2)
                do i = 1, ngc ! always fill the ghost cells at the bc face

                   data(hi(1)+i, j, k, :) = -data(hi(1)+1-i, j, k, :)

                end do
             end do
          end do
       elseif (bc_es_hi(1) .eq. 1) then ! Dirichlet
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
       if (bc_es_lo(2) .eq. 2) then ! Neumann
          do k = lo(3)-ngc_eff(3), hi(3)+ngc_eff(3)
             do j = 1, ngc ! always fill the ghost cells at the bc face
                do i = lo(1)-ngc_eff(1), hi(1)+ngc_eff(1)

                   data(i, lo(2)-j, k, :) = -data(i, lo(2)-1+j, k, :)

                end do
             end do
          end do
       elseif (bc_es_lo(2) .eq. 1) then ! Dirichlet
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
       if (bc_es_hi(2) .eq. 2) then ! Neumann
          do k = lo(3)-ngc_eff(3), hi(3)+ngc_eff(3)
             do j = 1, ngc ! always fill the ghost cells at the bc face
                do i = lo(1)-ngc_eff(1), hi(1)+ngc_eff(1)

                   data(i, hi(2)+j, k, :) = -data(i, hi(2)+1-j, k, :)

                end do
             end do
          end do
       elseif (bc_es_hi(2) .eq. 1) then ! Dirichlet
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
       if (bc_es_lo(3) .eq. 2) then ! Neumann
          do k = 1, ngc ! always fill the ghost cells at the bc face
             do j = lo(2)-ngc_eff(2), hi(2)+ngc_eff(2)
                do i = lo(1)-ngc_eff(1), hi(1)+ngc_eff(1)

                   data(i, j, lo(3)-k, :) = -data(i, j, lo(3)-1+k, :)

                end do
             end do
          end do
       elseif (bc_es_lo(3) .eq. 1) then ! Dirichlet
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
       if (bc_es_hi(3) .eq. 2) then ! no slip thermal
          do k = 1, ngc ! always fill the ghost cells at the bc face
             do j = lo(2)-ngc_eff(2), hi(2)+ngc_eff(2)
                do i = lo(1)-ngc_eff(1), hi(1)+ngc_eff(1)

                   data(i, j, hi(3)+k, :) = -data(i, j, hi(3)+1-k, :)

                end do
             end do
          end do
       elseif (bc_es_hi(3) .eq. 1) then ! no slip thermal
          do k = 1, ngc ! always fill the ghost cells at the bc face
             do j = lo(2)-ngc_eff(2), hi(2)+ngc_eff(2)
                do i = lo(1)-ngc_eff(1), hi(1)+ngc_eff(1)

                   data(i, j, hi(3)+k, :) = data(i, j, hi(3)+1-k, :)

                end do
             end do
          end do
       end if
    end if

  end subroutine fab_electricbc

  pure subroutine fab_potentialbc(lo,     hi,                  & ! dim(lo) == dim(hi) == 3
                                  dom_lo, dom_hi,              &
                                  data,   d_lo, d_hi, d_ncomp, & ! dim(d_lo) == dim(d_hi) == 3
                                  ngc, dim_fill_ghost)         &
                                  bind(C, name="fab_potentialbc")

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
       if (bc_es_lo(1) .eq. 2) then ! Neumann
          do k = lo(3)-ngc_eff(3), hi(3)+ngc_eff(3)
             do j = lo(2)-ngc_eff(2), hi(2)+ngc_eff(2)
                do i = 1, 1 ! always fill the ghost cells at the bc face

                   data(lo(1)-i, j, k, :) = data(lo(1)-1+i, j, k, :)

                end do
             end do
          end do
       elseif (bc_es_lo(1) .eq. 1) then !Dirichlet
          do k = lo(3)-ngc_eff(3), hi(3)+ngc_eff(3)
             do j = lo(2)-ngc_eff(2), hi(2)+ngc_eff(2)
                do i = 1, 1 ! always fill the ghost cells at the bc face

                   data(lo(1)-i, j, k, :) = -data(lo(1)-1+i, j, k, :) + 2*potential_lo(1)

                end do
             end do
          end do
       end if
    end if

    if (hi(1) .eq. dom_hi(1)) then ! upper bound
       if (bc_es_hi(1) .eq. 2) then !Neumann
          do k = lo(3)-ngc_eff(3), hi(3)+ngc_eff(3)
             do j = lo(2)-ngc_eff(2), hi(2)+ngc_eff(2)
                do i = 1, 1 ! always fill the ghost cells at the bc face

                   data(hi(1)+i, j, k, :) = data(hi(1)+1-i, j, k, :)

                end do
             end do
          end do
       elseif (bc_es_hi(1) .eq. 1) then ! Dirichlet
          do k = lo(3)-ngc_eff(3), hi(3)+ngc_eff(3)
             do j = lo(2)-ngc_eff(2), hi(2)+ngc_eff(2)
                do i = 1, 1 ! always fill the ghost cells at the bc face

                   data(hi(1)+i, j, k, :) = -data(hi(1)+1-i, j, k, :) + 2*potential_hi(1)

                end do
             end do
          end do
       end if

    end if


    !____________________________________________________________________________
    ! Apply BC to Y faces

    if (lo(2) .eq. dom_lo(2)) then ! lower bound
       if (bc_es_lo(2) .eq. 2) then ! Neumann
          do k = lo(3)-ngc_eff(3), hi(3)+ngc_eff(3)
             do j = 1, 1 ! always fill the ghost cells at the bc face
                do i = lo(1)-ngc_eff(1), hi(1)+ngc_eff(1)

                   data(i, lo(2)-j, k, :) = data(i, lo(2)-1+j, k, :)

                end do
             end do
          end do
       elseif (bc_es_lo(2) .eq. 1) then ! Dirichlet
          do k = lo(3)-ngc_eff(3), hi(3)+ngc_eff(3)
             do j = 1, 1 ! always fill the ghost cells at the bc face
                do i = lo(1)-ngc_eff(1), hi(1)+ngc_eff(1)

                   data(i, lo(2)-j, k, :) = -data(i, lo(2)-1+j, k, :) + 2*potential_lo(2)

                end do
             end do
          end do
       end if
    end if

    if (hi(2) .eq. dom_hi(2)) then ! upper bound
       if (bc_es_hi(2) .eq. 2) then ! Neumann
          do k = lo(3)-ngc_eff(3), hi(3)+ngc_eff(3)
             do j = 1, 1 ! always fill the ghost cells at the bc face
                do i = lo(1)-ngc_eff(1), hi(1)+ngc_eff(1)

                   data(i, hi(2)+j, k, :) = data(i, hi(2)+1-j, k, :)

                end do
             end do
          end do
       elseif (bc_es_hi(2) .eq. 1) then ! Dirichlet
          do k = lo(3)-ngc_eff(3), hi(3)+ngc_eff(3)
             do j = 1, 1 ! always fill the ghost cells at the bc face
                do i = lo(1)-ngc_eff(1), hi(1)+ngc_eff(1)

                   data(i, hi(2)+j, k, :) = -data(i, hi(2)+1-j, k, :) + 2*potential_hi(2)

                end do
             end do
          end do
       end if
    end if


    !____________________________________________________________________________
    ! Apply BC to Z faces

    if (lo(3) .eq. dom_lo(3)) then ! lower bound
       if (bc_es_lo(3) .eq. 2) then ! Neumann
          do k = 1, 1 ! always fill the ghost cells at the bc face
             do j = lo(2)-ngc_eff(2), hi(2)+ngc_eff(2)
                do i = lo(1)-ngc_eff(1), hi(1)+ngc_eff(1)

                   data(i, j, lo(3)-k, :) = data(i, j, lo(3)-1+k, :)

                end do
             end do
          end do
       elseif (bc_es_lo(3) .eq. 1) then ! Dirichlet
          do k = 1, 1 ! always fill the ghost cells at the bc face
             do j = lo(2)-ngc_eff(2), hi(2)+ngc_eff(2)
                do i = lo(1)-ngc_eff(1), hi(1)+ngc_eff(1)

                   data(i, j, lo(3)-k, :) = -data(i, j, lo(3)-1+k, :) + 2*potential_lo(3)

                end do
             end do
          end do
       end if

    end if

    if (hi(3) .eq. dom_hi(3)) then ! upper bound
       if (bc_es_hi(3) .eq. 2) then ! no slip thermal
          do k = 1, 1 ! always fill the ghost cells at the bc face
             do j = lo(2)-ngc_eff(2), hi(2)+ngc_eff(2)
                do i = lo(1)-ngc_eff(1), hi(1)+ngc_eff(1)

                   data(i, j, hi(3)+k, :) = data(i, j, hi(3)+1-k, :)

                end do
             end do
          end do
       elseif (bc_es_hi(3) .eq. 1) then ! no slip thermal
          do k = 1, 1 ! always fill the ghost cells at the bc face
             do j = lo(2)-ngc_eff(2), hi(2)+ngc_eff(2)
                do i = lo(1)-ngc_eff(1), hi(1)+ngc_eff(1)

                   data(i, j, hi(3)+k, :) = -data(i, j, hi(3)+1-k, :) + 2*potential_hi(3)

                end do
             end do
          end do
       end if
    end if

  end subroutine fab_potentialbc

  pure subroutine fab_potentialbc_solver(lo,     hi,                  & ! dim(lo) == dim(hi) == 3
                                         dom_lo, dom_hi,              &
                                         data,   d_lo, d_hi, d_ncomp, & ! dim(d_lo) == dim(d_hi) == 3
                                         ngc, dim_fill_ghost)         &
                                         bind(C, name="fab_potentialbc_solver")

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
       if (bc_es_lo(1) .eq. 2) then ! Neumann
          do k = lo(3)-ngc_eff(3), hi(3)+ngc_eff(3)
             do j = lo(2)-ngc_eff(2), hi(2)+ngc_eff(2)
                do i = 1, 1 ! always fill the ghost cells at the bc face

                   data(lo(1)-i, j, k, :) = potential_lo(1)

                end do
             end do
          end do
       elseif (bc_es_lo(1) .eq. 1) then !Dirichlet
          do k = lo(3)-ngc_eff(3), hi(3)+ngc_eff(3)
             do j = lo(2)-ngc_eff(2), hi(2)+ngc_eff(2)
                do i = 1, 1 ! always fill the ghost cells at the bc face

                   data(lo(1)-i, j, k, :) = potential_lo(1)

                end do
             end do
          end do
       end if
    end if

    if (hi(1) .eq. dom_hi(1)) then ! upper bound
       if (bc_es_hi(1) .eq. 2) then !Neumann
          do k = lo(3)-ngc_eff(3), hi(3)+ngc_eff(3)
             do j = lo(2)-ngc_eff(2), hi(2)+ngc_eff(2)
                do i = 1, 1 ! always fill the ghost cells at the bc face

                   data(hi(1)+i, j, k, :) = potential_hi(1)

                end do
             end do
          end do
       elseif (bc_es_hi(1) .eq. 1) then ! Dirichlet
          do k = lo(3)-ngc_eff(3), hi(3)+ngc_eff(3)
             do j = lo(2)-ngc_eff(2), hi(2)+ngc_eff(2)
                do i = 1, 1 ! always fill the ghost cells at the bc face

                   data(hi(1)+i, j, k, :) = potential_hi(1)

                end do
             end do
          end do
       end if

    end if


    !____________________________________________________________________________
    ! Apply BC to Y faces

    if (lo(2) .eq. dom_lo(2)) then ! lower bound
       if (bc_es_lo(2) .eq. 2) then ! Neumann
          do k = lo(3)-ngc_eff(3), hi(3)+ngc_eff(3)
             do j = 1, 1 ! always fill the ghost cells at the bc face
                do i = lo(1)-ngc_eff(1), hi(1)+ngc_eff(1)
                   
                   data(i, lo(2)-j, k, :) = potential_lo(2)

                end do
             end do
          end do
       elseif (bc_es_lo(2) .eq. 1) then ! Dirichlet
          do k = lo(3)-ngc_eff(3), hi(3)+ngc_eff(3)
             do j = 1, 1 ! always fill the ghost cells at the bc face
                do i = lo(1)-ngc_eff(1), hi(1)+ngc_eff(1)

                   data(i, lo(2)-j, k, :) = potential_lo(2)

                end do
             end do
          end do
       end if
    end if

    if (hi(2) .eq. dom_hi(2)) then ! upper bound
       if (bc_es_hi(2) .eq. 2) then ! Neumann
          do k = lo(3)-ngc_eff(3), hi(3)+ngc_eff(3)
             do j = 1, 1 ! always fill the ghost cells at the bc face
                do i = lo(1)-ngc_eff(1), hi(1)+ngc_eff(1)

                   data(i, hi(2)+j, k, :) = potential_hi(2)

                end do
             end do
          end do
       elseif (bc_es_hi(2) .eq. 1) then ! Dirichlet
          do k = lo(3)-ngc_eff(3), hi(3)+ngc_eff(3)
             do j = 1, 1 ! always fill the ghost cells at the bc face
                do i = lo(1)-ngc_eff(1), hi(1)+ngc_eff(1)

                   data(i, hi(2)+j, k, :) = potential_hi(2)

                end do
             end do
          end do
       end if
    end if


    !____________________________________________________________________________
    ! Apply BC to Z faces

    if (lo(3) .eq. dom_lo(3)) then ! lower bound
       if (bc_es_lo(3) .eq. 2) then ! Neumann
          do k = 1, 1 ! always fill the ghost cells at the bc face
             do j = lo(2)-ngc_eff(2), hi(2)+ngc_eff(2)
                do i = lo(1)-ngc_eff(1), hi(1)+ngc_eff(1)

                   data(i, j, lo(3)-k, :) = potential_lo(3)

                end do
             end do
          end do
       elseif (bc_es_lo(3) .eq. 1) then ! Dirichlet
          do k = 1, 1 ! always fill the ghost cells at the bc face
             do j = lo(2)-ngc_eff(2), hi(2)+ngc_eff(2)
                do i = lo(1)-ngc_eff(1), hi(1)+ngc_eff(1)

                   data(i, j, lo(3)-k, :) = potential_lo(3)

                end do
             end do
          end do
       end if

    end if

    if (hi(3) .eq. dom_hi(3)) then ! upper bound
       if (bc_es_hi(3) .eq. 2) then ! no slip thermal
          do k = 1, 1 ! always fill the ghost cells at the bc face
             do j = lo(2)-ngc_eff(2), hi(2)+ngc_eff(2)
                do i = lo(1)-ngc_eff(1), hi(1)+ngc_eff(1)

                   data(i, j, hi(3)+k, :) = potential_hi(3)

                end do
             end do
          end do
       elseif (bc_es_hi(3) .eq. 1) then ! no slip thermal
          do k = 1, 1 ! always fill the ghost cells at the bc face
             do j = lo(2)-ngc_eff(2), hi(2)+ngc_eff(2)
                do i = lo(1)-ngc_eff(1), hi(1)+ngc_eff(1)

                   data(i, j, hi(3)+k, :) = potential_hi(3)

                end do
             end do
          end do
       end if
    end if

  end subroutine fab_potentialbc_solver

  pure subroutine fab_chargebc(lo,     hi,                  & ! dim(lo) == dim(hi) == 3
       &                     dom_lo, dom_hi,              &
       &                     data,   d_lo, d_hi, d_ncomp, & ! dim(d_lo) == dim(d_hi) == 3
       &                     ngc, dim_fill_ghost)         &
       &                     bind(C, name="fab_chargebc")

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
       if (bc_es_lo(1) .eq. 2) then ! Neumann
          do k = lo(3)-ngc_eff(3), hi(3)+ngc_eff(3)
             do j = lo(2)-ngc_eff(2), hi(2)+ngc_eff(2)
                do i = 1, ngc

                   data(lo(1)-1+i, j, k, :) = data(lo(1)-1+i, j, k, :) + data(lo(1)-i, j, k, :)

                end do
             end do
          end do
       elseif (bc_es_lo(1) .eq. 1) then !Dirichlet
          do k = lo(3)-ngc_eff(3), hi(3)+ngc_eff(3)
             do j = lo(2)-ngc_eff(2), hi(2)+ngc_eff(2)
                do i = 1, ngc

                   data(lo(1)-1+i, j, k, :) = data(lo(1)-1+i, j, k, :) - data(lo(1)-i, j, k, :)

                end do
             end do
          end do
       end if
    end if

    if (hi(1) .eq. dom_hi(1)) then ! upper bound
       if (bc_es_hi(1) .eq. 2) then !Neumann
          do k = lo(3)-ngc_eff(3), hi(3)+ngc_eff(3)
             do j = lo(2)-ngc_eff(2), hi(2)+ngc_eff(2)
                do i = 1, ngc

                   data(hi(1)+1-i, j, k, :) = data(hi(1)+1-i, j, k, :) + data(hi(1)+i, j, k, :)

                end do
             end do
          end do
       elseif (bc_es_hi(1) .eq. 1) then ! Dirichlet
          do k = lo(3)-ngc_eff(3), hi(3)+ngc_eff(3)
             do j = lo(2)-ngc_eff(2), hi(2)+ngc_eff(2)
                do i = 1, ngc

                   data(hi(1)+1-i, j, k, :) = data(hi(1)+1-i, j, k, :) - data(hi(1)+i, j, k, :)

                end do
             end do
          end do
       end if

    end if


    !____________________________________________________________________________
    ! Apply BC to Y faces

    if (lo(2) .eq. dom_lo(2)) then ! lower bound
       if (bc_es_lo(2) .eq. 2) then ! Neumann
          do k = lo(3)-ngc_eff(3), hi(3)+ngc_eff(3)
             do j = 1, ngc
                do i = lo(1)-ngc_eff(1), hi(1)+ngc_eff(1)

                   data(i, lo(2)-1+j, k, :) = data(i, lo(2)-1+j, k, :) + data(i, lo(2)-j, k, :)

                end do
             end do
          end do
       elseif (bc_es_lo(2) .eq. 1) then ! Dirichlet
          do k = lo(3)-ngc_eff(3), hi(3)+ngc_eff(3)
             do j = 1, ngc
                do i = lo(1)-ngc_eff(1), hi(1)+ngc_eff(1)

                   data(i, lo(2)-1+j, k, :) = data(i, lo(2)-1+j, k, :) - data(i, lo(2)-j, k, :)

                end do
             end do
          end do
       end if
    end if

    if (hi(2) .eq. dom_hi(2)) then ! upper bound
       if (bc_es_hi(2) .eq. 2) then ! Neumann
          do k = lo(3)-ngc_eff(3), hi(3)+ngc_eff(3)
             do j = 1, ngc
                do i = lo(1)-ngc_eff(1), hi(1)+ngc_eff(1)

                   data(i, hi(2)+1-j, k, :) = data(i, hi(2)+1-j, k, :) + data(i, hi(2)+j, k, :)

                end do
             end do
          end do
       elseif (bc_es_hi(2) .eq. 1) then ! Dirichlet
          do k = lo(3)-ngc_eff(3), hi(3)+ngc_eff(3)
             do j = 1, ngc
                do i = lo(1)-ngc_eff(1), hi(1)+ngc_eff(1)

                   data(i, hi(2)+1-j, k, :) = data(i, hi(2)+1-j, k, :) - data(i, hi(2)+j, k, :)

                end do
             end do
          end do
       end if
    end if


    !____________________________________________________________________________
    ! Apply BC to Z faces

    if (lo(3) .eq. dom_lo(3)) then ! lower bound
       if (bc_es_lo(3) .eq. 2) then ! Neumann
          do k = 1, ngc
             do j = lo(2)-ngc_eff(2), hi(2)+ngc_eff(2)
                do i = lo(1)-ngc_eff(1), hi(1)+ngc_eff(1)

                   data(i, j, lo(3)-1+k, :) = data(i, j, lo(3)-1+k, :) + data(i, j, lo(3)-k, :)

                end do
             end do
          end do
       elseif (bc_es_lo(3) .eq. 1) then ! Dirichlet
          do k = 1, ngc
             do j = lo(2)-ngc_eff(2), hi(2)+ngc_eff(2)
                do i = lo(1)-ngc_eff(1), hi(1)+ngc_eff(1)

                   data(i, j, lo(3)-1+k, :) = data(i, j, lo(3)-1+k, :) - data(i, j, lo(3)-k, :)

                end do
             end do
          end do
       end if

    end if

    if (hi(3) .eq. dom_hi(3)) then ! upper bound
       if (bc_es_hi(3) .eq. 2) then ! no slip thermal
          do k = 1, ngc
             do j = lo(2)-ngc_eff(2), hi(2)+ngc_eff(2)
                do i = lo(1)-ngc_eff(1), hi(1)+ngc_eff(1)

                   data(i, j, hi(3)+1-k, :) = data(i, j, hi(3)+1-k, :) + data(i, j, hi(3)+k, :)

                end do
             end do
          end do
       elseif (bc_es_hi(3) .eq. 1) then ! no slip thermal
          do k = 1, ngc
             do j = lo(2)-ngc_eff(2), hi(2)+ngc_eff(2)
                do i = lo(1)-ngc_eff(1), hi(1)+ngc_eff(1)

                   data(i, j, hi(3)+1-k, :) = data(i, j, hi(3)+1-k, :) - data(i, j, hi(3)+k, :)

                end do
             end do
          end do
       end if
    end if

  end subroutine fab_chargebc

#endif


#if (AMREX_SPACEDIM == 2)

  pure subroutine fab_physbc_domainvel(lo,     hi,               & ! dim(lo) == dim(hi) == 3
       &                               dom_lo, dom_hi,           &
       &                               vel, v_lo, v_hi, v_ncomp, & ! dim(v_lo) == dim(v_hi) == 3
       &                               ngc, dim_fill_ghost, dd)      &
       &                               bind(C, name="fab_physbc_domainvel")

    integer,          intent(in   ) :: lo(3), hi(3), dom_lo(3), dom_hi(3), v_lo(3), v_hi(3), dd
    integer,          intent(in   ) :: v_ncomp, dim_fill_ghost(2)
    integer, value,   intent(in   ) :: ngc
    real(amrex_real), intent(inout) :: vel(v_lo(1):v_hi(1), &
         &                                 v_lo(2):v_hi(2), v_ncomp)

    ! ** loop indices
    integer :: i,j

    ! ** number of ghost cells to fill in each dimension
    integer, dimension(2) :: ngc_eff

    ngc_eff(:) = ngc*dim_fill_ghost(:)


    ! A wee note about limits for face-centered indices: face-centred boxes will
    ! have a hi(n) = dom_hi(n)+1 (where n is the direction of the face-centred
    ! quantity) and hi(m) = dom_hi(m) for all other direction


    !____________________________________________________________________________
    ! Apply BC to X faces
    if(dd .eq.  0) then
    if(lo(1) .eq. dom_lo(1)) then ! lower bound
       if(bc_vel_lo(1) .eq. 2) then ! no slip thermal

          do j = lo(2), hi(2)

             vel(lo(1), j, :) = 0

          end do

          do j = lo(2)-ngc_eff(2), hi(2)+ngc_eff(2)
             do i = 1, ngc ! always fill the ghost cells at the bc face

                ! Normal face-centered indices are symmetric
                vel(lo(1)-i, j, :) = -vel(lo(1)+i, j, :)

             end do
          end do

       end if
    end if

    if(hi(1) .eq. (dom_hi(1)+1)) then ! upper bound (note: +1)
       if(bc_vel_hi(1) .eq. 2) then ! no slip thermal

          do j = lo(2), hi(2)

             vel(hi(1), j, :) = 0

          end do

          do j = lo(2)-ngc_eff(2), hi(2)+ngc_eff(2)
             do i = 1, ngc ! always fill the ghost cells at the bc face

                ! Normal face-centered indices are symmetric
                vel(hi(1)+i, j, :) = -vel(hi(1)-i, j, :)

             end do
          end do

       end if
    end if
    endif


    !____________________________________________________________________________
    ! Apply BC to Y faces
    if(dd .eq.  1) then
    if(lo(2) .eq. dom_lo(2)) then ! lower bound
       if(bc_vel_lo(2) .eq. 2) then ! no slip thermal

          do i = lo(1), hi(1)

             vel(i, lo(2), :) = 0

          end do

          do j = 1, ngc ! always fill the ghost cells at the bc face
             do i = lo(1)-ngc_eff(1), hi(1)+ngc_eff(1)

                ! Normal face-centered indices are symmetric
                vel(i, lo(2)-j, :) = -vel(i, lo(2)+j, :)

             end do
          end do

       end if
    end if

    if(hi(2) .eq. (dom_hi(2)+1)) then ! upper bound (note: +1)
       if(bc_vel_hi(2) .eq. 2) then ! no slip thermal

          do i = lo(1), hi(1)

             vel(i, hi(2), :) = 0

          end do

          do j = 1, ngc ! always fill the ghost cells at the bc face
             do i = lo(1)-ngc_eff(1), hi(1)+ngc_eff(1)

                ! Normal face-centered indices are symmetric
                vel(i, hi(2)+j, :) = -vel(i, hi(2)-j, :)

             end do
          end do

       end if
    end if
    endif

  end subroutine fab_physbc_domainvel

  subroutine fab_physbc_domainstress(lo,     hi,               & ! dim(lo) == dim(hi) == 3
       &                               dom_lo, dom_hi,           &
       &                               stress, v_lo, v_hi, v_ncomp, & ! dim(v_lo) == dim(v_hi) == 3
       &                               ngc, dim_fill_ghost, dd)      &
       &                               bind(C, name="fab_physbc_domainstress")

    integer,          intent(in   ) :: lo(3), hi(3), dom_lo(3), dom_hi(3), v_lo(3), v_hi(3),dd
    integer,          intent(in   ) :: v_ncomp, dim_fill_ghost(2)
    integer, value,   intent(in   ) :: ngc
    real(amrex_real), intent(inout) :: stress(v_lo(1):v_hi(1), &
         &                                 v_lo(2):v_hi(2), v_ncomp)

    ! ** loop indices
    integer :: i,j

    ! ** number of ghost cells to fill in each dimension
    integer, dimension(2) :: ngc_eff, sliplo, sliphi

    ngc_eff(:) = ngc*dim_fill_ghost(:)

    do i=1,2
      if(bc_vel_lo(i) .eq. 1) then
        sliplo(i) = 1
      elseif(bc_vel_lo(i) .eq. 2) then
        sliplo(i) = -1
      endif
    enddo

    do i=1,2
      if(bc_vel_hi(i) .eq. 1) then
        sliphi(i) = 1
      elseif(bc_vel_hi(i) .eq. 2) then
        sliphi(i) = -1
      endif
    enddo

    ! A wee note about limits for face-centered indices: face-centred boxes will
    ! have a hi(n) = dom_hi(n)+1 (where n is the direction of the face-centred
    ! quantity) and hi(m) = dom_hi(m) for all other direction

    !____________________________________________________________________________
    ! Apply BC to X faces

    if(dd .eq.  0) then
    if(lo(1) .eq. dom_lo(1)) then ! lower bound
       if(bc_vel_lo(1) .ne. -1) then ! no slip thermal


             do j = lo(2), hi(2)

                stress(lo(1), j, :) = 0.5*(1d0 + sliplo(1))*stress(lo(1), j, :)

             end do



             do j = lo(2)-ngc_eff(2), hi(2)+ngc_eff(2)
                do i = 1, ngc ! always fill the ghost cells at the bc face

                   ! Normal face-centered indices are symmetric

                  stress(lo(1)+i, j, :) = stress(lo(1)+i, j, :) + sliplo(1)*stress(lo(1)-i, j, :)                    

                end do
             end do


       end if
    end if

    if(hi(1) .eq. (dom_hi(1)+1)) then ! upper bound (note: +1)
       if(bc_vel_hi(1) .ne. -1) then ! no slip thermal

             do j = lo(2), hi(2)

                stress(hi(1), j, :) = 0.5*(1d0 + sliphi(1))*stress(hi(1), j, :)

             end do

             do j = lo(2)-ngc_eff(2), hi(2)+ngc_eff(2)
                do i = 1, ngc ! always fill the ghost cells at the bc face

                   stress(hi(1)-i, j, :) = stress(hi(1)-i, j, :) + sliphi(1)*stress(hi(1)+i, j, :)

                end do
             end do

       end if
    end if
    endif


    !____________________________________________________________________________
    ! Apply BC to Y faces
    if(dd .eq.  1) then
    if(lo(2) .eq. dom_lo(2)) then ! lower bound
       if(bc_vel_lo(2) .ne. -1) then ! no slip thermal

             do i = lo(1), hi(1)

                stress(i, lo(2), :) = 0.5*(1d0 + sliplo(2))*stress(i, lo(2), :)

             end do

             do j = 1, ngc ! always fill the ghost cells at the bc face
                do i = lo(1)-ngc_eff(1), hi(1)+ngc_eff(1)

                   ! Normal face-centered indices are symmetric
                   stress(i, lo(2)+j, :) = stress(i, lo(2)+j, :) + sliplo(2)*stress(i, lo(2)-j, :)

                end do
             end do

       end if
    end if

    if(hi(2) .eq. (dom_hi(2)+1)) then ! upper bound (note: +1)
       if(bc_vel_hi(2) .ne. -1) then ! no slip thermal


             do i = lo(1), hi(1)

                stress(i, hi(2), :) = 0.5*(1d0 + sliphi(2))*stress(i, hi(2), :)

             end do

             do j = 1, ngc ! always fill the ghost cells at the bc face
                do i = lo(1)-ngc_eff(1), hi(1)+ngc_eff(1)

                   ! Normal face-centered indices are symmetric
                   stress(i, hi(2)-j, :) = stress(i, hi(2)-j, :) + sliphi(2)*stress(i, hi(2)+j, :)

                end do
             end do

       end if
    end if
    endif

  end subroutine fab_physbc_domainstress

#elif (AMREX_SPACEDIM == 3)

  pure subroutine fab_physbc_domainvel(lo,     hi,               & ! dim(lo) == dim(hi) == 3
       &                               dom_lo, dom_hi,           &
       &                               vel, v_lo, v_hi, v_ncomp, & ! dim(v_lo) == dim(v_hi) == 3
       &                               ngc, dim_fill_ghost, dd)      &
       &                               bind(C, name="fab_physbc_domainvel")

    integer,          intent(in   ) :: lo(3), hi(3), dom_lo(3), dom_hi(3), v_lo(3), v_hi(3), dd
    integer,          intent(in   ) :: v_ncomp, dim_fill_ghost(3)
    integer, value,   intent(in   ) :: ngc
    real(amrex_real), intent(inout) :: vel(v_lo(1):v_hi(1), &
         &                                 v_lo(2):v_hi(2), &
         &                                 v_lo(3):v_hi(3), v_ncomp)

    ! ** loop indices
    integer :: i,j,k

    ! ** number of ghost cells to fill in each dimension
    integer, dimension(3) :: ngc_eff

    ngc_eff(:) = ngc*dim_fill_ghost(:)


    ! A wee note about limits for face-centered indices: face-centred boxes will
    ! have a hi(n) = dom_hi(n)+1 (where n is the direction of the face-centred
    ! quantity) and hi(m) = dom_hi(m) for all other direction


    !____________________________________________________________________________
    ! Apply BC to X faces

    if(dd .eq.  0) then
    if(lo(1) .eq. dom_lo(1)) then ! lower bound
       if(bc_vel_lo(1) .eq. 2) then ! no slip thermal

          do k = lo(3), hi(3)
             do j = lo(2), hi(2)

                vel(lo(1), j, k, :) = 0

             end do
          end do

          do k = lo(3)-ngc_eff(3), hi(3)+ngc_eff(3)
             do j = lo(2)-ngc_eff(2), hi(2)+ngc_eff(2)
                do i = 1, ngc ! always fill the ghost cells at the bc face

                   ! Normal face-centered indices are symmetric
                   vel(lo(1)-i, j, k, :) = -vel(lo(1)+i, j, k, :)

                end do
             end do
          end do

       end if
    end if

    if(hi(1) .eq. (dom_hi(1)+1)) then ! upper bound (note: +1)
       if(bc_vel_hi(1) .eq. 2) then ! no slip thermal

          do k = lo(3), hi(3)
             do j = lo(2), hi(2)

                vel(hi(1), j, k, :) = 0

             end do
          end do

          do k = lo(3)-ngc_eff(3), hi(3)+ngc_eff(3)
             do j = lo(2)-ngc_eff(2), hi(2)+ngc_eff(2)
                do i = 1, ngc ! always fill the ghost cells at the bc face

                   ! Normal face-centered indices are symmetric
                   vel(hi(1)+i, j, k, :) = -vel(hi(1)-i, j, k, :)

                end do
             end do
          end do

       end if
    end if
    endif

    !____________________________________________________________________________
    ! Apply BC to Y faces

    if(dd .eq.  1) then
    if(lo(2) .eq. dom_lo(2)) then ! lower bound
       if(bc_vel_lo(2) .eq. 2) then ! no slip thermal

          do k = lo(3), hi(3)
             do i = lo(1), hi(1)

                vel(i, lo(2), k, :) = 0

             end do
          end do

          do k = lo(3)-ngc_eff(3), hi(3)+ngc_eff(3)
             do j = 1, ngc ! always fill the ghost cells at the bc face
                do i = lo(1)-ngc_eff(1), hi(1)+ngc_eff(1)

                   ! Normal face-centered indices are symmetric
                   vel(i, lo(2)-j, k, :) = -vel(i, lo(2)+j, k, :)

                end do
             end do
          end do

       end if
    end if

    if(hi(2) .eq. (dom_hi(2)+1)) then ! upper bound (note: +1)
       if(bc_vel_hi(2) .eq. 2) then ! no slip thermal

          do k = lo(3), hi(3)
             do i = lo(1), hi(1)

                vel(i, hi(2), k, :) = 0

             end do
          end do

          do k = lo(3)-ngc_eff(3), hi(3)+ngc_eff(3)
             do j = 1, ngc ! always fill the ghost cells at the bc face
                do i = lo(1)-ngc_eff(1), hi(1)+ngc_eff(1)

                   ! Normal face-centered indices are symmetric
                   vel(i, hi(2)+j, k, :) = -vel(i, hi(2)-j, k, :)

                end do
             end do
          end do

       end if
    end if
    endif


    !____________________________________________________________________________
    ! Apply BC to Z faces
    if(dd .eq.  2) then
    if(lo(3) .eq. dom_lo(3)) then ! lower bound
       if(bc_vel_lo(3) .eq. 2) then ! no slip thermal

          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                vel(i, j, lo(3), :) = 0

             end do
          end do

          do k = 1, ngc ! always fill the ghost cells at the bc face
             do j = lo(2)-ngc_eff(2), hi(2)+ngc_eff(2)
                do i = lo(1)-ngc_eff(1), hi(1)+ngc_eff(1)

                   ! Normal face-centered indices are symmetric
                   vel(i, j, lo(3)-k, :) = -vel(i, j, lo(3)+k, :)

                end do
             end do
          end do

       end if
    end if

    if(hi(3) .eq. (dom_hi(3)+1)) then ! upper bound (note: +1)
       if(bc_vel_hi(3) .eq. 2) then ! no slip thermal

          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                vel(i, j, hi(3), :) = 0

             end do
          end do

          do k = 1, ngc ! always fill the ghost cells at the bc face
             do j = lo(2)-ngc_eff(2), hi(2)+ngc_eff(2)
                do i = lo(1)-ngc_eff(1), hi(1)+ngc_eff(1)

                   ! Normal face-centered indices are symmetric
                   vel(i, j, hi(3)+k, :) = -vel(i, j, hi(3)-k, :)

                end do
             end do
          end do

       end if
    end if
    endif

  end subroutine fab_physbc_domainvel

  subroutine fab_physbc_domainstress(lo,     hi,               & ! dim(lo) == dim(hi) == 3
       &                               dom_lo, dom_hi,           &
       &                               stress, v_lo, v_hi, v_ncomp, & ! dim(v_lo) == dim(v_hi) == 3
       &                               ngc, dim_fill_ghost, dd)      &
       &                               bind(C, name="fab_physbc_domainstress")

    integer,          intent(in   ) :: lo(3), hi(3), dom_lo(3), dom_hi(3), v_lo(3), v_hi(3),dd
    integer,          intent(in   ) :: v_ncomp, dim_fill_ghost(3)
    integer, value,   intent(in   ) :: ngc
    real(amrex_real), intent(inout) :: stress(v_lo(1):v_hi(1), &
         &                                 v_lo(2):v_hi(2), &
         &                                 v_lo(3):v_hi(3), v_ncomp)

    ! ** loop indices
    integer :: i,j,k

    ! ** number of ghost cells to fill in each dimension
    integer, dimension(3) :: ngc_eff, sliplo, sliphi

    ngc_eff(:) = ngc*dim_fill_ghost(:)

!    do i=1,3
!      if(bc_vel_lo(i) .eq. 1) then
!        sliplo(i) = 1
!      elseif(bc_vel_lo(i) .eq. 2) then
!        sliplo(i) = -1
!      endif
!    enddo

!    do i=1,3
!      if(bc_vel_hi(i) .eq. 1) then
!        sliphi(i) = 1
!      elseif(bc_vel_hi(i) .eq. 2) then
!        sliphi(i) = -1
!      endif
!    enddo


    ! A wee note about limits for face-centered indices: face-centred boxes will
    ! have a hi(n) = dom_hi(n)+1 (where n is the direction of the face-centred
    ! quantity) and hi(m) = dom_hi(m) for all other direction


    !____________________________________________________________________________
    ! Apply BC to X faces

    if(dd .eq.  0) then
    if(lo(1) .eq. dom_lo(1)) then ! lower bound
       if(bc_vel_lo(1) .ne. -1) then ! no slip thermal

          do k = lo(3)-ngc_eff(3), hi(3)+ngc_eff(3)
             do j = lo(2)-ngc_eff(2), hi(2)+ngc_eff(2)

                stress(lo(1), j, k, :) = 0

             end do
          end do

          do k = lo(3)-ngc_eff(3), hi(3)+ngc_eff(3)
             do j = lo(2)-ngc_eff(2), hi(2)+ngc_eff(2)
                do i = 1, ngc ! always fill the ghost cells at the bc face

                   ! Normal face-centered indices are symmetric

                  stress(lo(1)+i, j, k, :) = stress(lo(1)+i, j, k, :) - stress(lo(1)-i, j, k, :)

                end do
             end do
          end do

       end if
    end if

    if(hi(1) .eq. (dom_hi(1)+1)) then ! upper bound (note: +1)
       if(bc_vel_hi(1) .ne. -1) then ! no slip thermal

          do k = lo(3)-ngc_eff(3), hi(3)+ngc_eff(3)
             do j = lo(2)-ngc_eff(2), hi(2)+ngc_eff(2)

                stress(hi(1), j, k, :) = 0

             end do
          end do

          do k = lo(3)-ngc_eff(3), hi(3)+ngc_eff(3)
             do j = lo(2)-ngc_eff(2), hi(2)+ngc_eff(2)
                do i = 1, ngc ! always fill the ghost cells at the bc face

                   stress(hi(1)-i, j, k, :) = stress(hi(1)-i, j, k, :) - stress(hi(1)+i, j, k, :)

                end do
             end do
          end do

       end if
    end if
    endif

    !____________________________________________________________________________
    ! Apply BC to Y faces
    if(dd .eq.  1) then
    if(lo(2) .eq. dom_lo(2)) then ! lower bound
       if(bc_vel_lo(2) .ne. -1) then ! no slip thermal

          do k = lo(3)-ngc_eff(3), hi(3)+ngc_eff(3)
             do i = lo(1)-ngc_eff(1), hi(1)+ngc_eff(1)

                stress(i, lo(2), k, :) = 0

                !print *, "Zeroing ", i, lo(2), k

             end do
          end do

          do k = lo(3)-ngc_eff(3), hi(3)+ngc_eff(3)
             do j = 1, ngc ! always fill the ghost cells at the bc face
                do i = lo(1)-ngc_eff(1), hi(1)+ngc_eff(1)

                   ! Normal face-centered indices are symmetric

                   !print *, "FORTRAN BC!", stress(i, lo(2)+j, k, :), stress(i, lo(2)-j, k, :)
                   stress(i, lo(2)+j, k, :) = stress(i, lo(2)+j, k, :) - stress(i, lo(2)-j, k, :)

                end do
             end do
          end do

       end if
    end if

    if(hi(2) .eq. (dom_hi(2)+1)) then ! upper bound (note: +1)
       if(bc_vel_hi(2) .ne. -1) then ! no slip thermal

          do k = lo(3)-ngc_eff(3), hi(3)+ngc_eff(3)
             do i = lo(1)-ngc_eff(1), hi(1)+ngc_eff(1)

                stress(i, hi(2), k, :) = 0

             end do
          end do

          do k = lo(3)-ngc_eff(3), hi(3)+ngc_eff(3)
             do j = 1, ngc ! always fill the ghost cells at the bc face
                do i = lo(1)-ngc_eff(1), hi(1)+ngc_eff(1)

                   ! Normal face-centered indices are symmetric
                   stress(i, hi(2)-j, k, :) = stress(i, hi(2)-j, k, :) - stress(i, hi(2)+j, k, :)

                end do
             end do
          end do

       end if
    end if
    endif


    !____________________________________________________________________________
    ! Apply BC to Z faces
    if(dd .eq.  2) then
    if(lo(3) .eq. dom_lo(3)) then ! lower bound
       if(bc_vel_lo(3) .ne. -1) then ! no slip thermal

          do j = lo(2)-ngc_eff(2), hi(2)+ngc_eff(2)
             do i = lo(1)-ngc_eff(1), hi(1)+ngc_eff(1)

                stress(i, j, lo(3), :) = 0

             end do
          end do

          do k = 1, ngc ! always fill the ghost cells at the bc face
             do j = lo(2)-ngc_eff(2), hi(2)+ngc_eff(2)
                do i = lo(1)-ngc_eff(1), hi(1)+ngc_eff(1)

                   ! Normal face-centered indices are symmetric
                   stress(i, j, lo(3)+k, :) = stress(i, j, lo(3)+k, :) - stress(i, j, lo(3)-k, :)

                end do
             end do
          end do

       end if
    end if

    if(hi(3) .eq. (dom_hi(3)+1)) then ! upper bound (note: +1)
       if(bc_vel_hi(3) .ne. -1) then ! no slip thermal

          do j = lo(2)-ngc_eff(2), hi(2)+ngc_eff(2)
             do i = lo(1)-ngc_eff(1), hi(1)+ngc_eff(1)

                stress(i, j, hi(3), :) = 0

             end do
          end do

          do k = 1, ngc ! always fill the ghost cells at the bc face
             do j = lo(2)-ngc_eff(2), hi(2)+ngc_eff(2)
                do i = lo(1)-ngc_eff(1), hi(1)+ngc_eff(1)

                   ! Normal face-centered indices are symmetric
                   stress(i, j, hi(3)-k, :) = stress(i, j, hi(3)-k, :) - stress(i, j, hi(3)+k, :)

                end do
             end do
          end do

       end if
    end if
  endif

  end subroutine fab_physbc_domainstress

#endif


#if (AMREX_SPACEDIM == 2)

  pure subroutine fab_physbc_macvel(lo,     hi,               & ! dim(lo) == dim(hi) == 3
       &                            dom_lo, dom_hi,           &
       &                            vel, v_lo, v_hi, v_ncomp, & ! dim(v_lo) == dim(v_hi) == 3
       &                            ngc, dim_fill_ghost, dd)      &
       &                            bind(C, name="fab_physbc_macvel")

    integer,          intent(in   ) :: lo(3), hi(3), dom_lo(3), dom_hi(3), v_lo(3), v_hi(3), dd
    integer,          intent(in   ) :: v_ncomp, dim_fill_ghost(2)
    integer, value,   intent(in   ) :: ngc
    real(amrex_real), intent(inout) :: vel(v_lo(1):v_hi(1), &
         &                                 v_lo(2):v_hi(2), v_ncomp)

    ! ** loop indices
    integer :: i,j

    ! ** number of ghost cells to fill in each dimension
    integer, dimension(2) :: ngc_eff

    ngc_eff(:) = ngc*dim_fill_ghost(:)

    !____________________________________________________________________________
    ! Apply BC to X faces
    if(dd .ne.  0) then
    if(lo(1) .eq. dom_lo(1)) then ! lower bound
       if(bc_vel_lo(1) .eq. 2) then ! no slip thermal

          do j = lo(2)-ngc_eff(2), hi(2)+ngc_eff(2)
             do i = 1, ngc ! always fill the ghost cells at the bc face

                vel(lo(1)-i, j, :) = -vel(lo(1)-1+i, j, :)

             end do
          end do

       end if
    end if

    if(hi(1) .eq. dom_hi(1)) then ! upper bound
       if(bc_vel_hi(1) .eq. 2) then ! no slip thermal

          do j = lo(2)-ngc_eff(2), hi(2)+ngc_eff(2)
             do i = 1, ngc ! always fill the ghost cells at the bc face

                vel(hi(1)+i, j, :) = -vel(hi(1)+1-i, j, :)

             end do
          end do

       end if
    end if
    endif

    !____________________________________________________________________________
    ! Apply BC to Y faces

    if(dd .ne.  1) then
    if(lo(2) .eq. dom_lo(2)) then ! lower bound
       if(bc_vel_lo(2) .eq. 2) then ! no slip thermal

          do j = 1, ngc ! always fill the ghost cells at the bc face
             do i = lo(1)-ngc_eff(1), hi(1)+ngc_eff(1)

                vel(i, lo(2)-j, :) = -vel(i, lo(2)-1+j, :)

             end do
          end do

       end if
    end if

    if(hi(2) .eq. dom_hi(2)) then ! upper bound
       if(bc_vel_hi(2) .eq. 2) then ! no slip thermal

          do j = 1, ngc ! always fill the ghost cells at the bc face
             do i = lo(1)-ngc_eff(1), hi(1)+ngc_eff(1)

                vel(i, hi(2)+j, :) = -vel(i, hi(2)+1-j, :)

             end do
          end do

       end if
    end if
    endif

  end subroutine fab_physbc_macvel

  pure subroutine fab_physbc_macstress(lo,     hi,               & ! dim(lo) == dim(hi) == 3
       &                            dom_lo, dom_hi,           &
       &                            stress, v_lo, v_hi, v_ncomp, & ! dim(v_lo) == dim(v_hi) == 3
       &                            ngc, dim_fill_ghost, dd)      &
       &                            bind(C, name="fab_physbc_macstress")

    integer,          intent(in   ) :: lo(3), hi(3), dom_lo(3), dom_hi(3), v_lo(3), v_hi(3), dd
    integer,          intent(in   ) :: v_ncomp, dim_fill_ghost(2)
    integer, value,   intent(in   ) :: ngc
    real(amrex_real), intent(inout) :: stress(v_lo(1):v_hi(1), &
         &                                 v_lo(2):v_hi(2), v_ncomp)

    ! ** loop indices
    integer :: i,j

    ! ** number of ghost cells to fill in each dimension
    integer, dimension(2) :: ngc_eff, sliplo, sliphi

    ngc_eff(:) = ngc*dim_fill_ghost(:)

    do i=1,2
      if(bc_vel_lo(i) .eq. 1) then
        sliplo(i) = 1
      elseif(bc_vel_lo(i) .eq. 2) then
        sliplo(i) = -1
      endif
    enddo

    do i=1,2
      if(bc_vel_hi(i) .eq. 1) then
        sliphi(i) = 1
      elseif(bc_vel_hi(i) .eq. 2) then
        sliphi(i) = -1
      endif
    enddo

    !____________________________________________________________________________
    ! Apply BC to X faces

    if(dd .ne.  0) then
    if(lo(1) .eq. dom_lo(1)) then ! lower bound
       if(bc_vel_lo(1) .ne. -1) then ! no slip thermal

             do j = lo(2)-ngc_eff(2), hi(2)+ngc_eff(2)
                do i = 1, ngc ! always fill the ghost cells at the bc face

                   stress(lo(1)-1+i, j, :) = stress(lo(1)-1+i, j, :) + sliplo(1)*stress(lo(1)-i, j, :)

                end do
             end do

       end if
    end if

    if(hi(1) .eq. dom_hi(1)) then ! upper bound
       if(bc_vel_hi(1) .ne. -1) then ! no slip thermal

             do j = lo(2)-ngc_eff(2), hi(2)+ngc_eff(2)
                do i = 1, ngc ! always fill the ghost cells at the bc face

                   stress(hi(1)+1-i, j, :) = stress(hi(1)+1-i, j, :) + sliphi(1)*stress(hi(1)+i, j, :)

                end do
             end do

       end if
    end if
    endif


    !____________________________________________________________________________
    ! Apply BC to Y faces
    if(dd .ne.  1) then
    if(lo(2) .eq. dom_lo(2)) then ! lower bound
       if(bc_vel_lo(2) .ne. -1) then ! no slip thermal


             do j = 1, ngc ! always fill the ghost cells at the bc face
                do i = lo(1)-ngc_eff(1), hi(1)+ngc_eff(1)

                   stress(i, lo(2)-1+j, :) = stress(i, lo(2)-1+j, :) + sliplo(2)*stress(i, lo(2)-j, :)

                end do
             end do


       end if
    end if

    if(hi(2) .eq. dom_hi(2)) then ! upper bound
       if(bc_vel_hi(2) .ne. -1) then ! no slip thermal

             do j = 1, ngc ! always fill the ghost cells at the bc face
                do i = lo(1)-ngc_eff(1), hi(1)+ngc_eff(1)

                   stress(i, hi(2)+1-j, :) = stress(i, hi(2)+1-j, :) + sliphi(2)*stress(i, hi(2)+j, :)

                end do
             end do

       end if
    end if
    endif

    !endif

  end subroutine fab_physbc_macstress

#elif (AMREX_SPACEDIM == 3)

  pure subroutine fab_physbc_macvel(lo,     hi,               & ! dim(lo) == dim(hi) == 3
       &                            dom_lo, dom_hi,           &
       &                            vel, v_lo, v_hi, v_ncomp, & ! dim(v_lo) == dim(v_hi) == 3
       &                            ngc, dim_fill_ghost,dd)      &
       &                            bind(C, name="fab_physbc_macvel")

    integer,          intent(in   ) :: lo(3), hi(3), dom_lo(3), dom_hi(3), v_lo(3), v_hi(3),dd
    integer,          intent(in   ) :: v_ncomp, dim_fill_ghost(3)
    integer, value,   intent(in   ) :: ngc
    real(amrex_real), intent(inout) :: vel(v_lo(1):v_hi(1), &
         &                                 v_lo(2):v_hi(2), &
         &                                 v_lo(3):v_hi(3), v_ncomp)

    ! ** loop indices
    integer :: i,j,k

    ! ** number of ghost cells to fill in each dimension
    integer, dimension(3) :: ngc_eff

    ngc_eff(:) = ngc*dim_fill_ghost(:)

    !____________________________________________________________________________
    ! Apply BC to X faces

    if(dd .ne.  0) then
    if(lo(1) .eq. dom_lo(1)) then ! lower bound
       if(bc_vel_lo(1) .eq. 2) then ! no slip thermal

          do k = lo(3)-ngc_eff(3), hi(3)+ngc_eff(3)
             do j = lo(2)-ngc_eff(2), hi(2)+ngc_eff(2)
                do i = 1, ngc ! always fill the ghost cells at the bc face

                   vel(lo(1)-i, j, k, :) = -vel(lo(1)-1+i, j, k, :)

                end do
             end do
          end do

       end if
    end if

    if(hi(1) .eq. dom_hi(1)) then ! upper bound
       if(bc_vel_hi(1) .eq. 2) then ! no slip thermal

          do k = lo(3)-ngc_eff(3), hi(3)+ngc_eff(3)
             do j = lo(2)-ngc_eff(2), hi(2)+ngc_eff(2)
                do i = 1, ngc ! always fill the ghost cells at the bc face

                   vel(hi(1)+i, j, k, :) = -vel(hi(1)+1-i, j, k, :)

                end do
             end do
          end do

       end if
    end if
    endif


    !____________________________________________________________________________
    ! Apply BC to Y faces

    if(dd .ne.  1) then
    if(lo(2) .eq. dom_lo(2)) then ! lower bound
       if(bc_vel_lo(2) .eq. 2) then ! no slip thermal

          do k = lo(3)-ngc_eff(3), hi(3)+ngc_eff(3)
             do j = 1, ngc ! always fill the ghost cells at the bc face
                do i = lo(1)-ngc_eff(1), hi(1)+ngc_eff(1)

                   vel(i, lo(2)-j, k, :) = -vel(i, lo(2)-1+j, k, :)

                end do
             end do
          end do

       end if
    end if

    if(hi(2) .eq. dom_hi(2)) then ! upper bound
       if(bc_vel_hi(2) .eq. 2) then ! no slip thermal

          do k = lo(3)-ngc_eff(3), hi(3)+ngc_eff(3)
             do j = 1, ngc ! always fill the ghost cells at the bc face
                do i = lo(1)-ngc_eff(1), hi(1)+ngc_eff(1)

                   vel(i, hi(2)+j, k, :) = -vel(i, hi(2)+1-j, k, :)

                end do
             end do
          end do

       end if
    end if
    endif

    !____________________________________________________________________________
    ! Apply BC to Z faces
    if(dd .ne.  2) then
    if(lo(3) .eq. dom_lo(3)) then ! lower bound
       if(bc_vel_lo(3) .eq. 2) then ! no slip thermal

          do k = 1, ngc ! always fill the ghost cells at the bc face
             do j = lo(2)-ngc_eff(2), hi(2)+ngc_eff(2)
                do i = lo(1)-ngc_eff(1), hi(1)+ngc_eff(1)

                   vel(i, j, lo(3)-k, :) = -vel(i, j, lo(3)-1+k, :)

                end do
             end do
          end do

       end if
    end if

    if(hi(3) .eq. dom_hi(3)) then ! upper bound
       if(bc_vel_hi(3) .eq. 2) then ! no slip thermal

          do k = 1, ngc ! always fill the ghost cells at the bc face
             do j = lo(2)-ngc_eff(2), hi(2)+ngc_eff(2)
                do i = lo(1)-ngc_eff(1), hi(1)+ngc_eff(1)

                   vel(i, j, hi(3)+k, :) = -vel(i, j, hi(3)+1-k, :)

                end do
             end do
          end do

       end if
    end if
    endif

  end subroutine fab_physbc_macvel

  pure subroutine fab_physbc_macstress(lo,     hi,               & ! dim(lo) == dim(hi) == 3
       &                            dom_lo, dom_hi,           &
       &                            stress, v_lo, v_hi, v_ncomp, & ! dim(v_lo) == dim(v_hi) == 3
       &                            ngc, dim_fill_ghost, dd)      &
       &                            bind(C, name="fab_physbc_macstress")

    integer,          intent(in   ) :: lo(3), hi(3), dom_lo(3), dom_hi(3), v_lo(3), v_hi(3), dd
    integer,          intent(in   ) :: v_ncomp, dim_fill_ghost(3)
    integer, value,   intent(in   ) :: ngc
    real(amrex_real), intent(inout) :: stress(v_lo(1):v_hi(1), &
         &                                 v_lo(2):v_hi(2), &
         &                                 v_lo(3):v_hi(3), v_ncomp)

    ! ** loop indices
    integer :: i,j,k

    ! ** number of ghost cells to fill in each dimension
    integer, dimension(3) :: ngc_eff, sliplo, sliphi

    ngc_eff(:) = ngc*dim_fill_ghost(:)

    do i=1,3
      if(bc_vel_lo(i) .eq. 1) then
        sliplo(i) = 1
      elseif(bc_vel_lo(i) .eq. 2) then
        sliplo(i) = -1
      endif
    enddo

    do i=1,2
      if(bc_vel_hi(i) .eq. 1) then
        sliphi(i) = 1
      elseif(bc_vel_hi(i) .eq. 2) then
        sliphi(i) = -1
      endif
    enddo

    !____________________________________________________________________________
    ! Apply BC to X faces

    if(dd .ne.  0) then
    if(lo(1) .eq. dom_lo(1)) then ! lower bound
       if(bc_vel_lo(1) .ne. -1) then ! no slip thermal

          do k = lo(3)-ngc_eff(3), hi(3)+ngc_eff(3)
             do j = lo(2)-ngc_eff(2), hi(2)+ngc_eff(2)
                do i = 1, ngc ! always fill the ghost cells at the bc face

                   stress(lo(1)-1+i, j, k, :) = stress(lo(1)-1+i, j, k, :) + sliplo(1)*stress(lo(1)-i, j, k, :)

                end do
             end do
          end do

       end if
    end if

    if(hi(1) .eq. dom_hi(1)) then ! upper bound
       if(bc_vel_hi(1) .ne. -1) then ! no slip thermal

          do k = lo(3)-ngc_eff(3), hi(3)+ngc_eff(3)
             do j = lo(2)-ngc_eff(2), hi(2)+ngc_eff(2)
                do i = 1, ngc ! always fill the ghost cells at the bc face

                   stress(hi(1)+1-i, j, k, :) = stress(hi(1)+1-i, j, k, :) + sliphi(1)*stress(hi(1)+i, j, k, :)

                end do
             end do
          end do

       end if
    end if
    endif


    !____________________________________________________________________________
    ! Apply BC to Y faces
    if(dd .ne.  1) then
    if(lo(2) .eq. dom_lo(2)) then ! lower bound
       if(bc_vel_lo(2) .ne. -1) then ! no slip thermal

          do k = lo(3)-ngc_eff(3), hi(3)+ngc_eff(3)
             do j = 1, ngc ! always fill the ghost cells at the bc face
                do i = lo(1)-ngc_eff(1), hi(1)+ngc_eff(1)

                   stress(i, lo(2)-1+j, k, :) = stress(i, lo(2)-1+j, k, :) + sliplo(2)*stress(i, lo(2)-j, k, :)

                end do
             end do
          end do

       end if
    end if

    if(hi(2) .eq. dom_hi(2)) then ! upper bound
       if(bc_vel_hi(2) .ne. -1) then ! no slip thermal

          do k = lo(3)-ngc_eff(3), hi(3)+ngc_eff(3)
             do j = 1, ngc ! always fill the ghost cells at the bc face
                do i = lo(1)-ngc_eff(1), hi(1)+ngc_eff(1)

                   stress(i, hi(2)+1-j, k, :) = stress(i, hi(2)+1-j, k, :) + sliphi(2)*stress(i, hi(2)+j, k, :)

                end do
             end do
          end do

       end if
    end if
    endif

    !____________________________________________________________________________
    ! Apply BC to Z faces
    if(dd .ne.  2) then
    if(lo(3) .eq. dom_lo(3)) then ! lower bound
       if(bc_vel_lo(3) .ne. -1) then ! no slip thermal

          do k = 1, ngc ! always fill the ghost cells at the bc face
             do j = lo(2)-ngc_eff(2), hi(2)+ngc_eff(2)
                do i = lo(1)-ngc_eff(1), hi(1)+ngc_eff(1)

                   stress(i, j, lo(3)-1+k, :) = stress(i, j, lo(3)-1+k, :) + sliplo(3)*stress(i, j, lo(3)-k, :)

                end do
             end do
          end do

       end if
    end if

    if(hi(3) .eq. dom_hi(3)) then ! upper bound
       if(bc_vel_hi(3) .ne. -1) then ! no slip thermal

          do k = 1, ngc ! always fill the ghost cells at the bc face
             do j = lo(2)-ngc_eff(2), hi(2)+ngc_eff(2)
                do i = lo(1)-ngc_eff(1), hi(1)+ngc_eff(1)

                   stress(i, j, hi(3)+1-k, :) = stress(i, j, hi(3)+1-k, :) + sliphi(3)*stress(i, j, hi(3)+k, :)

                end do
             end do
          end do

       end if
    end if
    endif

  end subroutine fab_physbc_macstress

#endif

end module multifab_physbc_module
