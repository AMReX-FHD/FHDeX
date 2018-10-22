module membrane_module

  use amrex_fort_module, only : amrex_real
  use common_namelist_module, only : ngc, nvars, nprimvars, nspecies, cell_depth, k_b, bc_lo, bc_hi, n_cells, hcv
  use conv_module
  implicit none

  private

  public :: do_ssa

contains

  subroutine do_ssa(lo, hi, cu, prim, dt) bind(c,name='do_ssa')

      implicit none

      integer,          intent(in      ) :: steps, lo(3), hi(3)
      double precision, intent(inout   ) :: del1, del2, del3

      double precision, intent(inout   ) :: cu(lo(1)-ngc:hi(1)+ngc,lo(2)-ngc:hi(2)+ngc,lo(3)-ngc:hi(3)+ngc, nvars)
      double precision, intent(inout   ) :: cumeans(lo(1)-ngc:hi(1)+ngc,lo(2)-ngc:hi(2)+ngc,lo(3)-ngc:hi(3)+ngc, nvars)

      double precision, intent(inout   ) :: prim(lo(1)-ngc:hi(1)+ngc,lo(2)-ngc:hi(2)+ngc,lo(3)-ngc:hi(3)+ngc, nprimvars)
      double precision, intent(inout   ) :: primmeans(lo(1)-ngc:hi(1)+ngc,lo(2)-ngc:hi(2)+ngc,lo(3)-ngc:hi(3)+ngc, nprimvars)

      !double precision fac1, fac2, fac3, test, pairfrac
      integer i,j,k,l, cells, ti
      double precision stepsminusone, stepsinv, densitymeaninv, fracvec(nspecies), massvec(nspecies)

      stepsminusone = steps - 1
      stepsinv = 1d0/steps

      del1 = 0
     ! del2 = 0
      del3 = 0

      

      do k = lo(3), hi(3)
        do j = lo(2), hi(2)
          do i = lo(1), hi(1)

            do l=1,nvars
              cumeans(i,j,k,l) = (cumeans(i,j,k,l)*stepsminusone + cu(i,j,k,l))*stepsinv
            enddo

            fracvec = cumeans(i,j,k,6:nvars)
            massvec = fracvec*cumeans(i,j,k,1)

            densitymeaninv = 1.0/cumeans(i,j,k,1)

            primmeans(i,j,k,1) = cumeans(i,j,k,1)
            primmeans(i,j,k,2) = cumeans(i,j,k,2)*densitymeaninv
            primmeans(i,j,k,3) = cumeans(i,j,k,3)*densitymeaninv
            primmeans(i,j,k,4) = cumeans(i,j,k,4)*densitymeaninv

            !print *, "Density: ", cu(i,j,k,1)

            call get_temperature(cumeans(i,j,k,5), massvec, primmeans(i,j,k,5))
            call get_pressure_gas(primmeans(i,j,k,6), fracvec, cumeans(i,j,k,1),cumeans(i,j,k,5))

            del1 = del1 + cu(i,j,k,1)
            !del2 = del2 + cu(i,j,k,2)
            del3 = del3 + cu(i,j,k,5)

          enddo
        enddo
      enddo

      !print *, hi    

      !if(hi(2) .lt. 25) then

          !print *, "steps: ", stepsminusone, stepsinv
          !print *, "density: ", cu(1,1,1,1), ", mean: ", cumeans(1,1,1,1), ", density: ", cu(0,1,0,1), ", mean: ", cumeans(0,1,0,1)
      !endif

       ti = 19
  !    tj = 0
  !    tk = 0

      if((ti .ge. lo(1)) .and. (ti .le. hi(1))) then
      
        del2 = cu(ti,0,0,2)-cumeans(ti,0,0,2)

      else
    
       del2 = 0

      endif

      !print *, "del2 in: ", del2
          
    end subroutine do_ssa

end module membrane_module
