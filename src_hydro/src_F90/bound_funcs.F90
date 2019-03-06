module bound_module

  use amrex_fort_module, only : amrex_real
  use common_namelist_module, only : bc_lo, bc_hi

  implicit none

  private

  ! public :: set_bc
  public :: fab_physbc
  public :: fab_physbc_domainvel
  public :: fab_physbc_macvel

contains

# if (AMREX_SPACEDIM == 2)

  subroutine fab_physbc(lo,     hi,                     & ! dim(lo) == dim(hi) == 3
       &                dom_lo, dom_hi,                 &
       &                pressure, p_lo, p_hi, p_ncomp,  & ! dim(p_lo) == dim(p_hi) == 3
       &                ngc, dim_fill_ghost)            &
       &                bind(C, name="fab_physbc")

    integer,          intent(in   ) :: lo(3), hi(3), dom_lo(3), dom_hi(3), &
         &                             p_lo(3), p_hi(3), p_ncomp
    integer,          intent(in   ) :: dim_fill_ghost(2)
    integer, value,   intent(in   ) :: ngc
    real(amrex_real), intent(inout) :: pressure(p_lo(1):p_hi(1), &
         &                                      p_lo(2):p_hi(2), p_ncomp)

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

                pressure(lo(1)-i, j, :) = pressure(lo(1)-1+i, j, :)

             end do
          end do

       end if
    end if

    if(hi(1) .eq. dom_hi(1)) then ! upper bound
       if(bc_hi(1) .eq. 2) then ! no slip thermal

          do j = lo(2)-ngc_eff(2), hi(2)+ngc_eff(2)
             do i = 1, ngc ! always fill the ghost cells at the bc face

                pressure(hi(1)+i, j, :) = pressure(hi(1)+1-i, j, :)

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

                pressure(i, lo(2)-j, :) = pressure(i, lo(2)-1+j, :)

             end do
          end do

       end if
    end if

    if(hi(2) .eq. dom_hi(2)) then ! upper bound
       if(bc_hi(2) .eq. 2) then ! no slip thermal

          do j = 1, ngc ! always fill the ghost cells at the bc face
             do i = lo(1)-ngc_eff(1), hi(1)+ngc_eff(1)

                pressure(i, hi(2)+j, :) = pressure(i, hi(2)+1-j, :)

             end do
          end do

       end if
    end if

  end subroutine fab_physbc

#elif (AMREX_SPACEDIM == 3)

  subroutine fab_physbc(lo,     hi,                     & ! dim(lo) == dim(hi) == 3
       &                dom_lo, dom_hi,                 &
       &                pressure, p_lo, p_hi, p_ncomp,  & ! dim(p_lo) == dim(p_hi) == 3
       &                ngc, dim_fill_ghost)            &
       &                bind(C, name="fab_physbc")

    integer,          intent(in   ) :: lo(3), hi(3), dom_lo(3), dom_hi(3), &
         &                             p_lo(3), p_hi(3), p_ncomp
    integer,          intent(in   ) :: dim_fill_ghost(3)
    integer, value,   intent(in   ) :: ngc
    real(amrex_real), intent(inout) :: pressure(p_lo(1):p_hi(1), &
         &                                      p_lo(2):p_hi(2), &
         &                                      p_lo(3):p_hi(3), p_ncomp)

    ! ** loop indices
    integer :: i,j,k

    ! ** number of ghost cells to fill in each dimension
    integer, dimension(3) :: ngc_eff

    ngc_eff(:) = ngc*dim_fill_ghost(:)


    !____________________________________________________________________________
    ! Apply BC to X faces

    if (lo(1) .eq. dom_lo(1)) then ! lower bound
       if(bc_lo(1) .eq. 2) then ! no slip thermal

          do k = lo(3)-ngc_eff(3), hi(3)+ngc_eff(3)
             do j = lo(2)-ngc_eff(2), hi(2)+ngc_eff(2)
                do i = 1, ngc ! always fill the ghost cells at the bc face

                   pressure(lo(1)-i, j, k, :) = pressure(lo(1)-1+i, j, k, :)

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

                   pressure(hi(1)+i, j, k, :) = pressure(hi(1)+1-i, j, k, :)

                end do
             end do
          end do

       end if
    end if


    !____________________________________________________________________________
    ! Apply BC to Y faces

    if(lo(2) .eq. dom_lo(2)) then ! lower bound
       if(bc_lo(2) .eq. 2) then ! no slip thermal

          do k = lo(3)-ngc_eff(3), hi(3)+ngc_eff(3)
             do j = 1, ngc ! always fill the ghost cells at the bc face
                do i = lo(1)-ngc_eff(1), hi(1)+ngc_eff(1)

                   pressure(i, lo(2)-j, k, :) = pressure(i, lo(2)-1+j, k, :)

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

                   pressure(i, hi(2)+j, k, :) = pressure(i, hi(2)+1-j, k, :)

                end do
             end do
          end do

       end if
    end if


    !____________________________________________________________________________
    ! Apply BC to Z faces

    if(lo(3) .eq. dom_lo(3)) then ! lower bound
       if(bc_lo(3) .eq. 2) then ! no slip thermal

          do k = 1, ngc ! always fill the ghost cells at the bc face
             do j = lo(2)-ngc_eff(2), hi(2)+ngc_eff(2)
                do i = lo(1)-ngc_eff(1), hi(1)+ngc_eff(1)

                   pressure(i, j, lo(3)-k, :) = pressure(i, j, lo(3)-1+k, :)

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

                   pressure(i, j, hi(3)+k, :) = pressure(i, j, hi(3)+1-k, :)

                end do
             end do
          end do

       end if
    end if

  end subroutine fab_physbc
#endif


#if (AMREX_SPACEDIM == 2)

  subroutine fab_physbc_domainvel(lo,     hi,               & ! dim(lo) == dim(hi) == 3
       &                          dom_lo, dom_hi,           &
       &                          vel, v_lo, v_hi, v_ncomp, & ! dim(v_lo) == dim(v_hi) == 3
       &                          ngc, dim_fill_ghost)      &
       &                          bind(C, name="fab_physbc_domainvel")

    integer,          intent(in   ) :: lo(3), hi(3), dom_lo(3), dom_hi(3), v_lo(3), v_hi(3)
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


    !____________________________________________________________________________
    ! Apply BC to Y faces

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

  end subroutine fab_physbc_domainvel

#elif (AMREX_SPACEDIM == 3)

  subroutine fab_physbc_domainvel(lo,     hi,               & ! dim(lo) == dim(hi) == 3
       &                          dom_lo, dom_hi,           &
       &                          vel, v_lo, v_hi, v_ncomp, & ! dim(v_lo) == dim(v_hi) == 3
       &                          ngc, dim_fill_ghost)      &
       &                          bind(C, name="fab_physbc_domainvel")

    integer,          intent(in   ) :: lo(3), hi(3), dom_lo(3), dom_hi(3), v_lo(3), v_hi(3)
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


    !____________________________________________________________________________
    ! Apply BC to Y faces

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

                   ! if ( -vel(i, lo(2)+j, k, 1) /= -vel(i, lo(2)+j, k, 1) ) then
                   !    write(*,*) lo, hi, ngc
                   !    write(*,*) "domainvel lo", i, lo(2)+j, k, -vel(i, lo(2)+j, k, :)

                   !    stop 0
                   ! end if

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

                   ! if ( -vel(i, hi(2)-j, k, 1) /= -vel(i, hi(2)-j, k, 1) ) then
                   !    write(*,*) lo, hi
                   !    write(*,*) "domainvel hi", i, hi(2)-j, k, -vel(i, hi(2)-j, k, :)

                   !    stop 0
                   ! end if


                end do
             end do
          end do

       end if
    end if


    !____________________________________________________________________________
    ! Apply BC to Z faces

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

  end subroutine fab_physbc_domainvel
#endif


#if (AMREX_SPACEDIM == 2)

  subroutine fab_physbc_macvel(lo,     hi,               & ! dim(lo) == dim(hi) == 3
       &                       dom_lo, dom_hi,           &
       &                       vel, v_lo, v_hi, v_ncomp, & ! dim(v_lo) == dim(v_hi) == 3
       &                       ngc, dim_fill_ghost)      &
       &                       bind(C, name="fab_physbc_macvel")

    integer,          intent(in   ) :: lo(3), hi(3), dom_lo(3), dom_hi(3), v_lo(3), v_hi(3)
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


    !____________________________________________________________________________
    ! Apply BC to Y faces

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

  end subroutine fab_physbc_macvel

#elif (AMREX_SPACEDIM == 3)

  subroutine fab_physbc_macvel(lo,     hi,               & ! dim(lo) == dim(hi) == 3
       &                       dom_lo, dom_hi,           &
       &                       vel, v_lo, v_hi, v_ncomp, & ! dim(v_lo) == dim(v_hi) == 3
       &                       ngc, dim_fill_ghost)      &
       &                       bind(C, name="fab_physbc_macvel")

    integer,          intent(in   ) :: lo(3), hi(3), dom_lo(3), dom_hi(3), v_lo(3), v_hi(3)
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


    !____________________________________________________________________________
    ! Apply BC to Y faces

    if(lo(2) .eq. dom_lo(2)) then ! lower bound
       if(bc_lo(2) .eq. 2) then ! no slip thermal

          do k = lo(3)-ngc_eff(3), hi(3)+ngc_eff(3)
             do j = 1, ngc ! always fill the ghost cells at the bc face
                do i = lo(1)-ngc_eff(1), hi(1)+ngc_eff(1)

                   vel(i, lo(2)-j, k, :) = -vel(i, lo(2)-1+j, k, :)

                   ! if ( -vel(i, lo(2)-1+j, k, 1) /= -vel(i, lo(2)-1+j, k, 1) ) then
                   !    write(*,*) lo, hi
                   !    write(*,*) "macvel", i, lo(2)-1+j, k, -vel(i, lo(2)-1+j, k, :)

                   !    stop 0
                   ! end if

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

                   ! if ( -vel(i, hi(2)+1-j, k, 1) /= -vel(i, hi(2)+1-j, k, 1) ) then
                   !    write(*,*) lo, hi
                   !    write(*,*) "macvel", i, hi(2)+1-j, k, -vel(i, hi(2)+1-j, k, :)

                   !    stop 0
                   ! end if


                end do
             end do
          end do

       end if
    end if


    !____________________________________________________________________________
    ! Apply BC to Z faces

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

  end subroutine fab_physbc_macvel

#endif


  ! subroutine set_bc(lo, hi, u_mac, v_mac,  w_mac ) bind(C, name="set_bc")

  !   integer,          intent(in   ) :: lo(3), hi(3)

  !   real(amrex_real), intent(inout) :: u_mac(lo(1)-ngc(1):hi(1)+ngc(1), &
  !        &                                   lo(2)-ngc(2):hi(2)+ngc(2), &
  !        &                                   lo(3)-ngc(3):hi(3)+ngc(3))
  !   real(amrex_real), intent(inout) :: v_mac(lo(1)-ngc(1):hi(1)+ngc(1), &
  !        &                                   lo(2)-ngc(2):hi(2)+ngc(2), &
  !        &                                   lo(3)-ngc(3):hi(3)+ngc(3))
  !   real(amrex_real), intent(inout) :: w_mac(lo(1)-ngc(1):hi(1)+ngc(1), &
  !        &                                   lo(2)-ngc(2):hi(2)+ngc(2), &
  !        &                                   lo(3)-ngc(3):hi(3)+ngc(3))


  !   ! ** loop indices
  !   integer :: i,j,k


  !   !____________________________________________________________________________
  !   ! Apply BC to X faces

  !   if(lo(1) .eq. 0) then ! lower bound
  !      if(bc_lo(1) .eq. 2) then ! no slip thermal

  !         do k = lo(3), hi(3)
  !            do j = lo(2), hi(2)

  !               u_mac(lo(1), j, k) = 0

  !            end do
  !         end do

  !         do k = lo(3)-ngc(3), hi(3)+ngc(3)
  !            do j = lo(2)-ngc(2), hi(2)+ngc(2)
  !               do i = 1, ngc(1)

  !                  ! Indices about x-face are symmetric
  !                  u_mac(lo(1)-i, j, k) = -u_mac(lo(1)+i, j, k)

  !                  v_mac(lo(1)-i, j, k) = -v_mac(lo(1)-1+i, j, k)
  !                  w_mac(lo(1)-i, j, k) = -w_mac(lo(1)-1+i, j, k)

  !               end do
  !            end do
  !         end do

  !      end if
  !   end if

  !   if(hi(1) .eq. (n_cells(1)-1)) then ! upper bound
  !      if(bc_hi(1) .eq. 2) then ! no slip thermal

  !         do k = lo(3), hi(3)
  !            do j = lo(2), hi(2)

  !               u_mac(hi(1), j, k) = 0

  !            end do
  !         end do

  !         do k = lo(3)-ngc(3), hi(3)+ngc(3)
  !            do j = lo(2)-ngc(2), hi(2)+ngc(2)
  !               do i = 1, ngc(1)

  !                  ! Indices about x-face are symmetric
  !                  u_mac(hi(1)+i, j, k) = -u_mac(hi(1)-i, j, k)

  !                  v_mac(hi(1)+i, j, k) = -v_mac(hi(1)+1-i, j, k)
  !                  w_mac(hi(1)+i, j, k) = -w_mac(hi(1)+1-i, j, k)

  !               end do
  !            end do
  !         end do

  !      end if
  !   end if


  !   !____________________________________________________________________________
  !   ! Apply BC to Y faces

  !   if(lo(2) .eq. 0) then ! lower bound
  !      if(bc_lo(2) .eq. 2) then ! no slip thermal

  !         do k = lo(3), hi(3)
  !            do i = lo(1), hi(1)

  !               v_mac(i, lo(2), k) = 0

  !            end do
  !         end do

  !         do k = lo(3)-ngc(3), hi(3)+ngc(3)
  !            do j = 1, ngc(2)
  !               do i = lo(1)-ngc(1), hi(1)+ngc(1)

  !                  ! Indices about y-face are symmetric
  !                  v_mac(i, lo(2)-j, k) = -v_mac(i, lo(2)+j, k)

  !                  u_mac(i, lo(2)-j, k) = -u_mac(i, lo(2)-1+j, k)
  !                  w_mac(i, lo(2)-j, k) = -w_mac(i, lo(2)-1+j, k)

  !                  if ( -w_mac(i, lo(2)-1+j, k) /= -w_mac(i, lo(2)-1+j, k) ) then
  !                     write(*,*) lo, hi
  !                     write(*,*) "vlo", i, lo(2)+j, k, -v_mac(i, lo(2)+j, k)
  !                     write(*,*) "ulo", i, lo(2)-1+j, k, -u_mac(i, lo(2)-1+j, k)
  !                     write(*,*) "wlo", i, lo(2)-1+j, k, -w_mac(i, lo(2)-1+j, k)
  !                  end if


  !               end do
  !            end do
  !         end do

  !      end if
  !   end if

  !   if(hi(2) .eq. (n_cells(2)-1)) then ! upper bound
  !      if(bc_hi(2) .eq. 2) then ! no slip thermal

  !         do k = lo(3), hi(3)
  !            do i = lo(1), hi(1)

  !               v_mac(i, hi(2), k) = 0

  !            end do
  !         end do

  !         do k = lo(3)-ngc(3), hi(3)+ngc(3)
  !            do j = 1, ngc(2)
  !               do i = lo(1)-ngc(1), hi(1)+ngc(1)

  !                  ! Indices about y-face are symmetric
  !                  v_mac(i, hi(2)+j, k) = -v_mac(i, hi(2)-j, k)

  !                  u_mac(i, hi(2)+j, k) = -u_mac(i, hi(2)+1-j, k)
  !                  w_mac(i, hi(2)+j, k) = -w_mac(i, hi(2)+1-j, k)

  !                  if ( -w_mac(i, hi(2)+1-j, k) /= -w_mac(i, hi(2)+1-j, k) ) then
  !                     write(*,*) lo, hi
  !                     write(*,*) "vhi", i, hi(2)-j, k, -v_mac(i, hi(2)-j, k)
  !                     write(*,*) "uhi", i, hi(2)+1-j, k, -u_mac(i, hi(2)+1-j, k)
  !                     write(*,*) "whi", i, hi(2)+1-j, k, -w_mac(i, hi(2)+1-j, k)
  !                  end if


  !               end do
  !            end do
  !         end do

  !      end if
  !   end if


  !   !____________________________________________________________________________
  !   ! Apply BC to Z faces

  !   if(lo(3) .eq. 0) then ! lower bound
  !      if(bc_lo(3) .eq. 2) then ! no slip thermal

  !         do j = lo(2), hi(2)
  !            do i = lo(1), hi(1)

  !               w_mac(i, j, lo(3)) = 0

  !            end do
  !         end do

  !         do k = 1, ngc(3)
  !            do j = lo(2)-ngc(2), hi(2)+ngc(2)
  !               do i = lo(1)-ngc(1), hi(1)+ngc(1)

  !                  ! Indices about z-face are symmetric
  !                  w_mac(i, j, lo(3)-k) = -w_mac(i, j, lo(3)+k)

  !                  u_mac(i, j, lo(3)-k) = -u_mac(i, j, lo(3)-1+k)
  !                  v_mac(i, j, lo(3)-k) = -v_mac(i, j, lo(3)-1+k)

  !               end do
  !            end do
  !         end do

  !      end if
  !   end if

  !   if(hi(3) .eq. (n_cells(3)-1)) then ! upper bound
  !      if(bc_hi(3) .eq. 2) then ! no slip thermal

  !         do j = lo(2), hi(2)
  !            do i = lo(1), hi(1)

  !               w_mac(i, j, hi(3)) = 0

  !            end do
  !         end do

  !         do k = 1, ngc(3)
  !            do j = lo(2)-ngc(2), hi(2)+ngc(2)
  !               do i = lo(1)-ngc(1), hi(1)+ngc(1)

  !                  ! Indices about z-face are symmetric
  !                  w_mac(i, j, hi(3)+k) = -w_mac(i, j, hi(3)-k)

  !                  u_mac(i, j, hi(3)+k) = -u_mac(i, j, hi(3)+1-k)
  !                  v_mac(i, j, hi(3)+k) = -v_mac(i, j, hi(3)+1-k)

  !               end do
  !            end do
  !         end do

  !      end if
  !   end if

  ! end subroutine set_bc

end module bound_module
