module membrane_module

  use amrex_fort_module, only : amrex_real
  use common_namelist_module, only : ngc, nvars, nprimvars, nspecies, cell_depth, k_b, &
       n_cells, hcv, membrane_cell, molmass,transmission
  use conv_module

  implicit none

  private

  public :: do_ssa

contains

  subroutine do_ssa(lo, hi, cu, prim, xflux, dx, dt) bind(c,name='do_ssa')

      implicit none

      integer,          intent(in      ) :: lo(3), hi(3)
      double precision, intent(in      ) :: dx(3), dt

      double precision, intent(in      ) :: cu(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3), nvars)
      double precision, intent(in      ) :: prim(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3), nprimvars)

      double precision, intent(inout   ) :: xflux(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3), 6)


      !double precision fac1, fac2, fac3, test, pairfrac
      integer i,j,k,l, fluxcount
      double precision vol, area, mm, rr, hole, tm, massflux, energyflux, taul, taur, rn1, rn2, rn3, rn4, ratel, rater, nl, nr, theta, energy

      if((lo(1) .eq. membrane_cell) .or. (hi(1)+1 .eq. membrane_cell)) then

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

             !print *, "tm: ", tm, " dt: ", dt

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

              fluxcount = fluxcount + 1
                
              else
                !right to left 

                energy = -k_B*prim(membrane_cell,j,k,5)*log(rn1*rn2)

                massflux = massflux - mm
                energyflux = energyflux - energy

              !fluxcount = fluxcount + 1

              endif

              theta = -log(rn4)/(ratel + rater)

              tm = tm + theta

             enddo
             xflux(membrane_cell,j,k,1) = massflux/vol
             xflux(membrane_cell,j,k,5) = energyflux/vol

             !print *, "ssa rho: ", xflux(membrane_cell,j,k,1)
             !print *, "ssa E: ", xflux(membrane_cell,j,k,5)

              !print *, "rho flux: ", massflux/vol
              !print *, "E flux: ", energyflux/vol
!                print *, "lo(1): ", lo(1), ", hi(1): ", hi(1)
!                print *, "fluxcount: ", fluxcount
!                print *, "saving: ", xflux(membrane_cell,j,k,5)
             print *, "Fluxcount: :", fluxcount
          
          enddo
        enddo

      endif
        
    end subroutine do_ssa

  subroutine do_langevin(lo, hi, cu, prim, xflux, dx, dt) bind(c,name='do_langevin')

      implicit none

      integer,          intent(in      ) :: lo(3), hi(3)
      double precision, intent(in      ) :: dx(3), dt

      double precision, intent(in      ) :: cu(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3), nvars)
      double precision, intent(in      ) :: prim(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3), nprimvars)

      double precision, intent(inout   ) :: xflux(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3), 6)


      !double precision fac1, fac2, fac3, test, pairfrac
      integer i,j,k,l, fluxcount
      double precision vol, area, mm, tl, tr, rhol, rhor, fac1,fac3, fac5, sqrttl, sqrttr, um, nm, uv, nv, cross, corr, rn1, rn2, rn3

      mm = molmass(1)/(6.02d23)

        !print *, mm

      fac5 = transmission*(k_b**2.5)*6.0/sqrt(2*mm*3.142)      
      fac3 = transmission*(k_b**1.5)*2.0/sqrt(2*mm*3.142)
      fac1 = transmission*(k_b**0.5)/sqrt(2*mm*3.142)

      if((lo(1) .eq. membrane_cell) .or. (hi(1)+1 .eq. membrane_cell)) then

#if (AMREX_SPACEDIM == 3)
        vol = dx(1)*dx(2)*dx(3)
        area = dx(2)*dx(3)
#else
        vol = dx(1)*dx(2)*cell_depth
        area = dx(2)*cell_depth
#endif

        do  k=lo(3),hi(3)
          do  j=lo(2),hi(2)

             rhol = cu(membrane_cell-1,j,k,1)
             rhor = cu(membrane_cell,j,k,1)

             tl = prim(membrane_cell-1,j,k,5)
             tr = prim(membrane_cell,j,k,5)

             sqrttl = sqrt(tl)
             sqrttr = sqrt(tr)

             um = fac3*(sqrttl*tl*rhol-sqrttr*tr*rhor)
             nm = fac1*(sqrttl*rhol-sqrttr*rhor)

             uv = fac5*(sqrttl*tl*tl*rhol+sqrttr*tr*tr*rhor)
             nv = fac1*(sqrttl*rhol+sqrttr*rhor)

             cross = fac3*(sqrttl*tl*rhol+sqrttr*tr*rhor)

             corr = cross/(sqrt(uv)*sqrt(nv))

                !print *, corr

             rn1 = get_fhd_normal_func()
             rn2 = get_fhd_normal_func()
             rn3 = rn1*corr + sqrt(1-corr**2)*rn2
            
             !xflux(membrane_cell,j,k,1) = (dt*area*nm + sqrt(dt*area*mm*nv)*rn1)/vol
             !xflux(membrane_cell,j,k,5) = (dt*area*um + sqrt(dt*area*mm*uv)*rn3)/(vol*mm)



             xflux(membrane_cell,j,k,1) = (dt*area*nm + sqrt(dt*area*mm*nv)*rn1)/vol
             xflux(membrane_cell,j,k,5) = (dt*area*um + sqrt(dt*area*mm*uv)*rn3)/(vol*mm)


             !print *, dt*area*nm, dt*area*nv

             !xflux(membrane_cell,j,k,1) = mm*(nm )/vol
             !xflux(membrane_cell,j,k,5) = (um)/vol

           ! print *, "rho flux: ", dt*area*(nm)/vol
           ! print *, "E flux: ", dt*area*(um)/(vol*mm)

             !print *, "langevin rho: ",  (-k_b*mm*l11/nav)*nd + (1.5*k_b*mm*l11/tav - mm*l12/(tav**2))*td
             !print *, "langevin E: ", (-k_b*l21/nav)*nd + (1.5*k_b*l21/tav - l12/(tav**2))*td

          enddo
        enddo

      endif


!      mm = molmass(1)/(6.02d23)

!      bigA = transmission*sqrt(mm*k_B/(2.0*3.142))

!      if((lo(1) .eq. membrane_cell) .or. (hi(1)+1 .eq. membrane_cell)) then

!#if (AMREX_SPACEDIM == 3)
!        vol = dx(1)*dx(2)*dx(3)
!        area = dx(2)*dx(3)
!#else
!        vol = dx(1)*dx(2)*cell_depth
!        area = dx(2)*cell_depth
!#endif

!        do  k=lo(3),hi(3)
!          do  j=lo(2),hi(2)

!             nl = cu(membrane_cell-1,j,k,1)/mm
!             nr = cu(membrane_cell,j,k,1)/mm

!             tl = prim(membrane_cell-1,j,k,5)
!             tr = prim(membrane_cell,j,k,5)

!             nav = (nl+nr)/2
!             tav = (tl+tr)/2

!             nd = nr-nl
!             td = tr-tl

!             l11 = bigA*nav*sqrt(tav)/k_b
!             l12 = 2*bigA*nav*tav*sqrt(tav)
!             l21 = l12
!             l22 = 4.5*k_b*bigA*nav*tav**2.5             

!             !rn1 = get_uniform_func()
!             
!             xflux(membrane_cell,j,k,1) = 0.1*((-k_b*mm*l11/nav)*nd + (1.5*k_b*mm*l11/tav - mm*l12/(tav**2))*td)/vol
!             xflux(membrane_cell,j,k,5) = 0.1*((-k_b*l21/nav)*nd + (1.5*k_b*l21/tav - l22/(tav**2))*td)/vol

!             !print *, "langevin rho: ",  (-k_b*mm*l11/nav)*nd + (1.5*k_b*mm*l11/tav - mm*l12/(tav**2))*td
!             !print *, "langevin E: ", (-k_b*l21/nav)*nd + (1.5*k_b*l21/tav - l12/(tav**2))*td

!          enddo
!        enddo

!      endif
        
        
    end subroutine do_langevin


  subroutine apply_effusion(lo, hi, cu, xflux, dx, dt) bind(c,name='apply_effusion')

      implicit none

      integer,          intent(in      ) :: lo(3), hi(3)
      double precision, intent(in      ) :: dx(3), dt

      double precision, intent(inout   ) :: cu(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3), nvars)

      double precision, intent(in      ) :: xflux(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3), 6)

      integer i,j,k

      if(lo(1) .eq. membrane_cell) then

        do  k=lo(3),hi(3)
          do  j=lo(2),hi(2)

            !print *, "old energy: ", cu(membrane_cell,j,k,5)
            !print *, "adding: ", xflux(membrane_cell,j,k,5)

            cu(membrane_cell,j,k,1) = cu(membrane_cell,j,k,1) + xflux(membrane_cell,j,k,1)
            cu(membrane_cell,j,k,5) = cu(membrane_cell,j,k,5) + xflux(membrane_cell,j,k,5)

                
            !print *, "new energy: ", cu(membrane_cell,j,k,5)
            ! print *, "Effusing!"

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
