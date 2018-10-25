module membrane_module

  use amrex_fort_module, only : amrex_real
  use common_namelist_module, only : ngc, nvars, nprimvars, nspecies, cell_depth, k_b, bc_lo, bc_hi, n_cells, hcv, membrane_cell, molmass,transmission
  use conv_module
  use rng_functions_module
  implicit none

  private

  public :: do_ssa

contains

  subroutine do_ssa(lo, hi, cu, prim, xflux, dx, dt) bind(c,name='do_ssa')

      implicit none

      integer,          intent(in      ) :: lo(3), hi(3)
      double precision, intent(in      ) :: dx(3), dt

      double precision, intent(in      ) :: cu(lo(1)-ngc:hi(1)+ngc,lo(2)-ngc:hi(2)+ngc,lo(3)-ngc:hi(3)+ngc, nvars)
      double precision, intent(in      ) :: prim(lo(1)-ngc:hi(1)+ngc,lo(2)-ngc:hi(2)+ngc,lo(3)-ngc:hi(3)+ngc, nprimvars)

      double precision, intent(inout   ) :: xflux(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3), 6)


      !double precision fac1, fac2, fac3, test, pairfrac
      integer i,j,k,l, fluxcount
      double precision vol, area, mm, rr, hole, tm, massflux, energyflux, taul, taur, rn1, rn2, rn3, rn4, ratel, rater, nl, nr, theta, energy


      if((lo(1) .eq. membrane_cell) .or. (hi(2)+1 .eq. membrane_cell)) then

#if (AMREX_SPACEDIM == 3)
        vol = dx(1)*dx(2)*dx(3)
        area = dx(2)*dx(3)
#else
        vol = dx(1)*dx(2)*cell_depth
        area = dx(2)*cell_depth
#endif
        mm = molmass(1)/(6.02d23)
        rr = k_B/mm
        hole = transmission*area

        do  k=lo(3),hi(3)
          do  j=lo(2),hi(2)

             tm = 0
             massflux = 0
             energyflux = 0
             fluxcount = 0

             nl = cu(membrane_cell-1,j,k,1)/mm
             nr = cu(membrane_cell,j,k,1)/mm

             taul = (1d0/(hole*nl))*sqrt(2d0*3.142/(rr*prim(membrane_cell-1,j,k,5)))
             taur = (1d0/(hole*nr))*sqrt(2d0*3.142/(rr*prim(membrane_cell,j,k,5)))

             ratel = 1d0/taul
             rater = 1d0/taur

             rn1 = get_uniform_func()
             theta = -log(rn1)/(ratel + rater)
             tm = tm + theta

             do while (tm .lt. dt)

              rn1 = get_uniform_func()
              rn2 = get_uniform_func()
              rn3 = get_uniform_func()
              rn4 = get_uniform_func()

              if(ratel/(ratel+rater) .gt. rn3) then
                !left to right
                energy = -k_B*prim(membrane_cell-1,j,k,5)*log(rn1*rn2)

                massflux = massflux + mm
                energyflux = energyflux + energy
                
              else
                !right to left 

                energy = -k_B*prim(membrane_cell,j,k,5)*log(rn1*rn2)

                massflux = massflux - mm
                energyflux = energyflux - energy

              endif

              theta = -log(rn4)/(ratel + rater)

              tm = tm + theta
              fluxcount = fluxcount + 1

             enddo
             xflux(membrane_cell,j,k,1) = massflux/vol
             xflux(membrane_cell,j,k,5) = energyflux/vol
          
          enddo
        enddo

      endif
        
    end subroutine do_ssa


  subroutine apply_effusion(lo, hi, cu, xflux, dx, dt) bind(c,name='apply_effusion')

      implicit none

      integer,          intent(in      ) :: lo(3), hi(3)
      double precision, intent(in      ) :: dx(3), dt

      double precision, intent(inout   ) :: cu(lo(1)-ngc:hi(1)+ngc,lo(2)-ngc:hi(2)+ngc,lo(3)-ngc:hi(3)+ngc, nvars)

      double precision, intent(in      ) :: xflux(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3), 6)

      integer i,j,k

      if(lo(1) .eq. membrane_cell) then

        do  k=lo(3),hi(3)
          do  j=lo(2),hi(2)

            cu(membrane_cell,j,k,1) = cu(membrane_cell,j,k,1) + xflux(membrane_cell,j,k,1)
            cu(membrane_cell,j,k,5) = cu(membrane_cell,j,k,5) + xflux(membrane_cell,j,k,5)

                

              !print *, "Effusing!"

            if((cu(membrane_cell,j,k,1) .lt. 0) .or. (cu(membrane_cell,j,k,5) .lt. 0)) then
    
              print *, "Negative effusion removed"

              cu(membrane_cell,j,k,1) = cu(membrane_cell,j,k,1) - xflux(membrane_cell,j,k,1)
              cu(membrane_cell,j,k,5) = cu(membrane_cell,j,k,5) - xflux(membrane_cell,j,k,5)

            endif


          enddo
        enddo

      endif


      if(hi(1) .eq. membrane_cell-1) then

        do  k=lo(3),hi(3)
          do  j=lo(2),hi(2)

            cu(membrane_cell-1,j,k,1) = cu(membrane_cell-1,j,k,1) - xflux(membrane_cell,j,k,1)
            cu(membrane_cell-1,j,k,5) = cu(membrane_cell-1,j,k,5) - xflux(membrane_cell,j,k,5)

              !print *, "Effusing!"

            if((cu(membrane_cell-1,j,k,1) .lt. 0) .or. (cu(membrane_cell-1,j,k,5) .lt. 0)) then
    
              print *, "Negative effusion removed"

              cu(membrane_cell-1,j,k,1) = cu(membrane_cell-1,j,k,1) + xflux(membrane_cell,j,k,1)
              cu(membrane_cell-1,j,k,5) = cu(membrane_cell-1,j,k,5) + xflux(membrane_cell,j,k,5)

            endif


          enddo
        enddo

      endif
        
    end subroutine apply_effusion

end module membrane_module
