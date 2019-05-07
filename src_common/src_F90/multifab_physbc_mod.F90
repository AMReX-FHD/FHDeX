module multifab_physbc_module

  use amrex_fort_module,      only : amrex_real
  use common_namelist_module, only : bc_lo, bc_hi

  implicit none

  private

  public :: fab_physbc
  public :: fab_electricbc
  public :: fab_physbc_domainvel
  public :: fab_physbc_macvel
  public :: fab_physbc_domainstress
  public :: fab_physbc_macstress


contains

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
       if(bc_lo(1) .eq. 2) then ! no slip thermal

          do j = lo(2)-ngc_eff(2), hi(2)+ngc_eff(2)
             do i = 1, ngc ! always fill the ghost cells at the bc face

                data(lo(1)-i, j, :) = data(lo(1)-1+i, j, :)

             end do
          end do

       end if
    end if

    if(hi(1) .eq. dom_hi(1)) then ! upper bound
       if(bc_hi(1) .eq. 2) then ! no slip thermal

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
       if(bc_lo(2) .eq. 2) then ! no slip thermal

          do j = 1, ngc ! always fill the ghost cells at the bc face
             do i = lo(1)-ngc_eff(1), hi(1)+ngc_eff(1)

                data(i, lo(2)-j, :) = data(i, lo(2)-1+j, :)

             end do
          end do

       end if
    end if

    if(hi(2) .eq. dom_hi(2)) then ! upper bound
       if(bc_hi(2) .eq. 2) then ! no slip thermal

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
       if(bc_lo(1) .eq. 2) then ! no slip thermal

          do j = lo(2)-ngc_eff(2), hi(2)+ngc_eff(2)
             do i = 1, ngc ! always fill the ghost cells at the bc face

                data(lo(1)-i, j, :) = -data(lo(1)-1+i, j, :)

             end do
          end do

       end if
    end if

    if(hi(1) .eq. dom_hi(1)) then ! upper bound
       if(bc_hi(1) .eq. 2) then ! no slip thermal

          do j = lo(2)-ngc_eff(2), hi(2)+ngc_eff(2)
             do i = 1, ngc ! always fill the ghost cells at the bc face

                data(hi(1)+i, j, :) = -data(hi(1)+1-i, j, :)

             end do
          end do

       end if
    end if


    !____________________________________________________________________________
    ! Apply BC to Y faces

    if(lo(2) .eq. dom_lo(2)) then ! lower bound
       if(bc_lo(2) .eq. 2) then ! no slip thermal

          do j = 1, ngc ! always fill the ghost cells at the bc face
             do i = lo(1)-ngc_eff(1), hi(1)+ngc_eff(1)

                data(i, lo(2)-j, :) = -data(i, lo(2)-1+j, :)

             end do
          end do

       end if
    end if

    if(hi(2) .eq. dom_hi(2)) then ! upper bound
       if(bc_hi(2) .eq. 2) then ! no slip thermal

          do j = 1, ngc ! always fill the ghost cells at the bc face
             do i = lo(1)-ngc_eff(1), hi(1)+ngc_eff(1)

                data(i, hi(2)+j, :) = -data(i, hi(2)+1-j, :)

             end do
          end do

       end if
    end if

  end subroutine fab_electricbc

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
       if (bc_lo(1) .eq. 2) then ! no slip thermal

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
       if (bc_hi(1) .eq. 2) then ! no slip thermal

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
       if (bc_lo(2) .eq. 2) then ! no slip thermal

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
       if (bc_hi(2) .eq. 2) then ! no slip thermal

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
       if (bc_lo(3) .eq. 2) then ! no slip thermal

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
       if (bc_hi(3) .eq. 2) then ! no slip thermal

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
       if (bc_lo(1) .eq. 2) then ! no slip thermal

          do k = lo(3)-ngc_eff(3), hi(3)+ngc_eff(3)
             do j = lo(2)-ngc_eff(2), hi(2)+ngc_eff(2)
                do i = 1, ngc ! always fill the ghost cells at the bc face

                   data(lo(1)-i, j, k, :) = -data(lo(1)-1+i, j, k, :)

                end do
             end do
          end do

       end if
    end if

    if (hi(1) .eq. dom_hi(1)) then ! upper bound
       if (bc_hi(1) .eq. 2) then ! no slip thermal

          do k = lo(3)-ngc_eff(3), hi(3)+ngc_eff(3)
             do j = lo(2)-ngc_eff(2), hi(2)+ngc_eff(2)
                do i = 1, ngc ! always fill the ghost cells at the bc face

                   data(hi(1)+i, j, k, :) = -data(hi(1)+1-i, j, k, :)

                end do
             end do
          end do

       end if
    end if


    !____________________________________________________________________________
    ! Apply BC to Y faces

    if (lo(2) .eq. dom_lo(2)) then ! lower bound
       if (bc_lo(2) .eq. 2) then ! no slip thermal

          do k = lo(3)-ngc_eff(3), hi(3)+ngc_eff(3)
             do j = 1, ngc ! always fill the ghost cells at the bc face
                do i = lo(1)-ngc_eff(1), hi(1)+ngc_eff(1)

                   data(i, lo(2)-j, k, :) = -data(i, lo(2)-1+j, k, :)

                end do
             end do
          end do

       end if
    end if

    if (hi(2) .eq. dom_hi(2)) then ! upper bound
       if (bc_hi(2) .eq. 2) then ! no slip thermal

          do k = lo(3)-ngc_eff(3), hi(3)+ngc_eff(3)
             do j = 1, ngc ! always fill the ghost cells at the bc face
                do i = lo(1)-ngc_eff(1), hi(1)+ngc_eff(1)

                   data(i, hi(2)+j, k, :) = -data(i, hi(2)+1-j, k, :)

                end do
             end do
          end do

       end if
    end if


    !____________________________________________________________________________
    ! Apply BC to Z faces

    if (lo(3) .eq. dom_lo(3)) then ! lower bound
       if (bc_lo(3) .eq. 2) then ! no slip thermal

          do k = 1, ngc ! always fill the ghost cells at the bc face
             do j = lo(2)-ngc_eff(2), hi(2)+ngc_eff(2)
                do i = lo(1)-ngc_eff(1), hi(1)+ngc_eff(1)

                   data(i, j, lo(3)-k, :) = -data(i, j, lo(3)-1+k, :)

                end do
             end do
          end do

       end if
    end if

    if (hi(3) .eq. dom_hi(3)) then ! upper bound
       if (bc_hi(3) .eq. 2) then ! no slip thermal

          do k = 1, ngc ! always fill the ghost cells at the bc face
             do j = lo(2)-ngc_eff(2), hi(2)+ngc_eff(2)
                do i = lo(1)-ngc_eff(1), hi(1)+ngc_eff(1)

                   data(i, j, hi(3)+k, :) = -data(i, j, hi(3)+1-k, :)

                end do
             end do
          end do

       end if
    end if

  end subroutine fab_electricbc

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
       if(bc_lo(1) .eq. 2) then ! no slip thermal

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
       if(bc_hi(1) .eq. 2) then ! no slip thermal

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
       if(bc_lo(2) .eq. 2) then ! no slip thermal

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
       if(bc_hi(2) .eq. 2) then ! no slip thermal

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
       if(bc_lo(1) .eq. 2) then ! no slip thermal

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
       if(bc_hi(1) .eq. 2) then ! no slip thermal

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
       if(bc_lo(2) .eq. 2) then ! no slip thermal

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
       if(bc_hi(2) .eq. 2) then ! no slip thermal

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
       if(bc_lo(3) .eq. 2) then ! no slip thermal

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
       if(bc_hi(3) .eq. 2) then ! no slip thermal

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
       if(bc_lo(1) .eq. 2) then ! no slip thermal

          do j = lo(2)-ngc_eff(2), hi(2)+ngc_eff(2)
             do i = 1, ngc ! always fill the ghost cells at the bc face

                vel(lo(1)-i, j, :) = -vel(lo(1)-1+i, j, :)

             end do
          end do

       end if
    end if

    if(hi(1) .eq. dom_hi(1)) then ! upper bound
       if(bc_hi(1) .eq. 2) then ! no slip thermal

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
       if(bc_lo(2) .eq. 2) then ! no slip thermal

          do j = 1, ngc ! always fill the ghost cells at the bc face
             do i = lo(1)-ngc_eff(1), hi(1)+ngc_eff(1)

                vel(i, lo(2)-j, :) = -vel(i, lo(2)-1+j, :)

             end do
          end do

       end if
    end if

    if(hi(2) .eq. dom_hi(2)) then ! upper bound
       if(bc_hi(2) .eq. 2) then ! no slip thermal

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
    integer :: i,j,k

    ! ** number of ghost cells to fill in each dimension
    integer, dimension(3) :: ngc_eff

    ngc_eff(:) = ngc*dim_fill_ghost(:)

    !____________________________________________________________________________
    ! Apply BC to X faces

    if(dd .ne.  0) then
    if(lo(1) .eq. dom_lo(1)) then ! lower bound
       if(bc_lo(1) .eq. 2) then ! no slip thermal

          do k = lo(3)-ngc_eff(3), hi(3)+ngc_eff(3)
             do j = lo(2)-ngc_eff(2), hi(2)+ngc_eff(2)
                do i = 1, ngc ! always fill the ghost cells at the bc face

                   stress(lo(1)-1+i, j, k, :) = stress(lo(1)-1+i, j, k, :) - stress(lo(1)-i, j, k, :)

                end do
             end do
          end do

       end if
    end if

    if(hi(1) .eq. dom_hi(1)) then ! upper bound
       if(bc_hi(1) .eq. 2) then ! no slip thermal

          do k = lo(3)-ngc_eff(3), hi(3)+ngc_eff(3)
             do j = lo(2)-ngc_eff(2), hi(2)+ngc_eff(2)
                do i = 1, ngc ! always fill the ghost cells at the bc face

                   stress(hi(1)+1-i, j, k, :) = stress(hi(1)+1-i, j, k, :) - stress(hi(1)+i, j, k, :)

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
       if(bc_lo(2) .eq. 2) then ! no slip thermal

          do k = lo(3)-ngc_eff(3), hi(3)+ngc_eff(3)
             do j = 1, ngc ! always fill the ghost cells at the bc face
                do i = lo(1)-ngc_eff(1), hi(1)+ngc_eff(1)

                   stress(i, lo(2)-1+j, k, :) = stress(i, lo(2)-1+j, k, :) - stress(i, lo(2)-j, k, :)

                end do
             end do
          end do

       end if
    end if

    if(hi(2) .eq. dom_hi(2)) then ! upper bound
       if(bc_hi(2) .eq. 2) then ! no slip thermal

          do k = lo(3)-ngc_eff(3), hi(3)+ngc_eff(3)
             do j = 1, ngc ! always fill the ghost cells at the bc face
                do i = lo(1)-ngc_eff(1), hi(1)+ngc_eff(1)

                   stress(i, hi(2)+1-j, k, :) = stress(i, hi(2)+1-j, k, :) - stress(i, hi(2)+j, k, :)

                end do
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
       if(bc_lo(1) .eq. 2) then ! no slip thermal

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
       if(bc_hi(1) .eq. 2) then ! no slip thermal

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
       if(bc_lo(2) .eq. 2) then ! no slip thermal

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
       if(bc_hi(2) .eq. 2) then ! no slip thermal

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
       if(bc_lo(3) .eq. 2) then ! no slip thermal

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
       if(bc_hi(3) .eq. 2) then ! no slip thermal

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
    integer, dimension(3) :: ngc_eff

    ngc_eff(:) = ngc*dim_fill_ghost(:)

    !____________________________________________________________________________
    ! Apply BC to X faces

    if(dd .ne.  0) then
    if(lo(1) .eq. dom_lo(1)) then ! lower bound
       if(bc_lo(1) .eq. 2) then ! no slip thermal

          do k = lo(3)-ngc_eff(3), hi(3)+ngc_eff(3)
             do j = lo(2)-ngc_eff(2), hi(2)+ngc_eff(2)
                do i = 1, ngc ! always fill the ghost cells at the bc face

                   stress(lo(1)-1+i, j, k, :) = stress(lo(1)-1+i, j, k, :) - stress(lo(1)-i, j, k, :)

                end do
             end do
          end do

       end if
    end if

    if(hi(1) .eq. dom_hi(1)) then ! upper bound
       if(bc_hi(1) .eq. 2) then ! no slip thermal

          do k = lo(3)-ngc_eff(3), hi(3)+ngc_eff(3)
             do j = lo(2)-ngc_eff(2), hi(2)+ngc_eff(2)
                do i = 1, ngc ! always fill the ghost cells at the bc face

                   stress(hi(1)+1-i, j, k, :) = stress(hi(1)+1-i, j, k, :) - stress(hi(1)+i, j, k, :)

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
       if(bc_lo(2) .eq. 2) then ! no slip thermal

          do k = lo(3)-ngc_eff(3), hi(3)+ngc_eff(3)
             do j = 1, ngc ! always fill the ghost cells at the bc face
                do i = lo(1)-ngc_eff(1), hi(1)+ngc_eff(1)

                   stress(i, lo(2)-1+j, k, :) = stress(i, lo(2)-1+j, k, :) - stress(i, lo(2)-j, k, :)

                end do
             end do
          end do

       end if
    end if

    if(hi(2) .eq. dom_hi(2)) then ! upper bound
       if(bc_hi(2) .eq. 2) then ! no slip thermal

          do k = lo(3)-ngc_eff(3), hi(3)+ngc_eff(3)
             do j = 1, ngc ! always fill the ghost cells at the bc face
                do i = lo(1)-ngc_eff(1), hi(1)+ngc_eff(1)

                   stress(i, hi(2)+1-j, k, :) = stress(i, hi(2)+1-j, k, :) - stress(i, hi(2)+j, k, :)

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
       if(bc_lo(3) .eq. 2) then ! no slip thermal

          do k = 1, ngc ! always fill the ghost cells at the bc face
             do j = lo(2)-ngc_eff(2), hi(2)+ngc_eff(2)
                do i = lo(1)-ngc_eff(1), hi(1)+ngc_eff(1)

                   stress(i, j, lo(3)-1+k, :) = stress(i, j, lo(3)-1+k, :) - stress(i, j, lo(3)-k, :)

                end do
             end do
          end do

       end if
    end if

    if(hi(3) .eq. dom_hi(3)) then ! upper bound
       if(bc_hi(3) .eq. 2) then ! no slip thermal

          do k = 1, ngc ! always fill the ghost cells at the bc face
             do j = lo(2)-ngc_eff(2), hi(2)+ngc_eff(2)
                do i = lo(1)-ngc_eff(1), hi(1)+ngc_eff(1)

                   stress(i, j, hi(3)+1-k, :) = stress(i, j, hi(3)+1-k, :) - stress(i, j, hi(3)+k, :)

                end do
             end do
          end do

       end if
    end if
    endif

  end subroutine fab_physbc_macstress

#endif

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
    integer, dimension(3) :: ngc_eff

    ngc_eff(:) = ngc*dim_fill_ghost(:)


    ! A wee note about limits for face-centered indices: face-centred boxes will
    ! have a hi(n) = dom_hi(n)+1 (where n is the direction of the face-centred
    ! quantity) and hi(m) = dom_hi(m) for all other direction


    !____________________________________________________________________________
    ! Apply BC to X faces
    
    if(dd .eq.  0) then
    if(lo(1) .eq. dom_lo(1)) then ! lower bound
       if(bc_lo(1) .eq. 2) then ! no slip thermal

          do k = lo(3), hi(3)
             do j = lo(2), hi(2)

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
       if(bc_hi(1) .eq. 2) then ! no slip thermal

          do k = lo(3), hi(3)
             do j = lo(2), hi(2)

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
       if(bc_lo(2) .eq. 2) then ! no slip thermal

          do k = lo(3), hi(3)
             do i = lo(1), hi(1)

                stress(i, lo(2), k, :) = 0

             end do
          end do

          do k = lo(3)-ngc_eff(3), hi(3)+ngc_eff(3)
             do j = 1, ngc ! always fill the ghost cells at the bc face
                do i = lo(1)-ngc_eff(1), hi(1)+ngc_eff(1)

                   ! Normal face-centered indices are symmetric
                   stress(i, lo(2)+j, k, :) = stress(i, lo(2)+j, k, :) - stress(i, lo(2)-j, k, :)

                end do
             end do
          end do

       end if
    end if

    if(hi(2) .eq. (dom_hi(2)+1)) then ! upper bound (note: +1)
       if(bc_hi(2) .eq. 2) then ! no slip thermal

          do k = lo(3), hi(3)
             do i = lo(1), hi(1)

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
       if(bc_lo(3) .eq. 2) then ! no slip thermal

          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

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
       if(bc_hi(3) .eq. 2) then ! no slip thermal

          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

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



end module multifab_physbc_module
