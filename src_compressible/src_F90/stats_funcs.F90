module stats_module

  use amrex_fort_module, only : amrex_real
  use common_namelist_module, only : ngc, nvars, nprimvars, nspecies, cell_depth, k_b, &
       n_cells, hcv, cross_cell
  use conv_module
  implicit none

  private

  public :: multifab_yzav

contains

  subroutine multifab_yzav(lo, hi, fabin, fabout, comps) bind(c,name='multifab_yzav')

      implicit none

      integer,          intent(in      ) :: lo(3), hi(3), comps

      double precision, intent(in      ) :: fabin(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3), comps)
      double precision, intent(inout   ) :: fabout(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3), comps)

      integer i,j,k,l, counter
      double precision holder

      do l = 1, comps

      do i = lo(1), hi(1)

        holder = 0;
        counter = 0;

        do k = lo(3), hi(3)
          do j = lo(2), hi(2)

           holder = holder + fabin(i,j,k,l)
           counter = counter + 1

          enddo
        enddo

        holder = holder/counter

        do k = lo(3), hi(3)
          do j = lo(2), hi(2)

           fabout(i,j,k,l) = holder

          enddo
        enddo

      enddo
      enddo

    end subroutine multifab_yzav


end module stats_module
