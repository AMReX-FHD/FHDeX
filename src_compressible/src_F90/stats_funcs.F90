module stats_module

  use amrex_fort_module, only : amrex_real
  use common_namelist_module, only : ngc, nvars, nprimvars, nspecies, cell_depth, k_b, bc_lo, bc_hi, n_cells, hcv
  use conv_module
  implicit none

  private

  public :: evaluate_means, evaluate_corrs

contains

  subroutine evaluate_means(lo, hi, cu, cumeans, prim, primmeans, steps, del1, del2) bind(c,name='evaluate_means')

      implicit none

      integer,          intent(in      ) :: steps, lo(3), hi(3)
      double precision, intent(inout   ) :: del1, del2

      double precision, intent(inout   ) :: cu(lo(1)-ngc:hi(1)+ngc,lo(2)-ngc:hi(2)+ngc,lo(3)-ngc:hi(3)+ngc, nvars)
      double precision, intent(inout   ) :: cumeans(lo(1)-ngc:hi(1)+ngc,lo(2)-ngc:hi(2)+ngc,lo(3)-ngc:hi(3)+ngc, nvars)

      double precision, intent(inout   ) :: prim(lo(1)-ngc:hi(1)+ngc,lo(2)-ngc:hi(2)+ngc,lo(3)-ngc:hi(3)+ngc, nprimvars)
      double precision, intent(inout   ) :: primmeans(lo(1)-ngc:hi(1)+ngc,lo(2)-ngc:hi(2)+ngc,lo(3)-ngc:hi(3)+ngc, nprimvars)

      !double precision fac1, fac2, fac3, test, pairfrac
      integer i,j,k,l
      double precision stepsminusone, stepsinv, densitymeaninv, fracvec(nspecies), massvec(nspecies)

      stepsminusone = steps - 1
      stepsinv = 1d0/steps

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

            call get_temperature(cumeans(i,j,k,5), massvec, primmeans(i,j,k,5))
            call get_pressure_gas(primmeans(i,j,k,6), fracvec, cumeans(i,j,k,1),cumeans(i,j,k,5))

          enddo
        enddo
      enddo

      !print *, hi    

      !if(hi(2) .lt. 25) then

          !print *, "steps: ", stepsminusone, stepsinv
          !print *, "density: ", cu(1,1,1,1), ", mean: ", cumeans(1,1,1,1), ", density: ", cu(0,1,0,1), ", mean: ", cumeans(0,1,0,1)
      !endif

  !    ti = 19
  !    tj = 0
  !    tk = 0

  !    if((ti .ge. mlo(1)) .and. (ti .le. mhi(1))) then
  !    
  !      del2 = pxmean(ti,tj,tk)

  !    else
  !  
  !      del2 = 0

  !    endif

      !print *, "del1 in: ", del1
          
    end subroutine evaluate_means

  subroutine evaluate_corrs(lo, hi, cu, cumeans, cuvars, prim, primmeans, primvars, steps, del1, del2) bind(c,name='evaluate_corrs')

      !use iso_c_binding, only: c_ptr, c_int, c_f_pointer

      implicit none

      integer,          intent(in      ) :: steps, lo(3), hi(3)
      double precision, intent(inout   ) :: del1, del2

      double precision, intent(inout   ) :: cu(lo(1)-ngc:hi(1)+ngc,lo(2)-ngc:hi(2)+ngc,lo(3)-ngc:hi(3)+ngc, nvars)
      double precision, intent(inout   ) :: cumeans(lo(1)-ngc:hi(1)+ngc,lo(2)-ngc:hi(2)+ngc,lo(3)-ngc:hi(3)+ngc, nvars)
      double precision, intent(inout   ) :: cuvars(lo(1)-ngc:hi(1)+ngc,lo(2)-ngc:hi(2)+ngc,lo(3)-ngc:hi(3)+ngc, nvars)

      double precision, intent(inout   ) :: prim(lo(1)-ngc:hi(1)+ngc,lo(2)-ngc:hi(2)+ngc,lo(3)-ngc:hi(3)+ngc, nprimvars)
      double precision, intent(inout   ) :: primmeans(lo(1)-ngc:hi(1)+ngc,lo(2)-ngc:hi(2)+ngc,lo(3)-ngc:hi(3)+ngc, nprimvars)
      double precision, intent(inout   ) :: primvars(lo(1)-ngc:hi(1)+ngc,lo(2)-ngc:hi(2)+ngc,lo(3)-ngc:hi(3)+ngc, nprimvars + 4)

      !double precision fac1, fac2, fac3, test, pairfrac
      integer i,j,k,l
      double precision stepsminusone, stepsinv, cv, cvinv, rr, delg, qmean, delpx, delpy, delpz, delrho, delvelx, delvely, delvelz, delenergy, densitymeaninv

      stepsminusone = steps - 1
      stepsinv = 1d0/steps

      rr = 1d0
      cvinv = 2.0/(3.0*rr)
      cv = 1.0/cvinv

      do k = lo(3), hi(3)
        do j = lo(2), hi(2)
          do i = lo(1), hi(1)

            !Vars

            delrho = cu(i,j,k,1) - cumeans(i,j,k,1)

            cuvars(i,j,k,1) = (cuvars(i,j,k,1)*stepsminusone + delrho**2)*stepsinv
            !vars(i,j,k,1) = delrho**2

          enddo
        enddo
      enddo


      do k = lo(3), hi(3)
        do j = lo(2), hi(2)
          do i = lo(1), hi(1)

          cv = 0
          do l=1,nspecies
            cv = cv + hcv(l)*cumeans(i,j,k,5+l)
          enddo

          cvinv = 1.0/cv

          !Vars
          qmean = cv*primmeans(i,j,k,5)-0.5*(primmeans(i,j,k,2)**2 + primmeans(i,j,k,3)**2 + primmeans(i,j,k,4)**2)

          densitymeaninv = 1.0/cumeans(i,j,k,1)

          delrho = cu(i,j,k,1) - cumeans(i,j,k,1)
          delpx = cu(i,j,k,2) - cumeans(i,j,k,2)
          delpy = cu(i,j,k,3) - cumeans(i,j,k,3)
          delpz = cu(i,j,k,4) - cumeans(i,j,k,4)
          delenergy = cu(i,j,k,5) - cumeans(i,j,k,5)

          cuvars(i,j,k,1) = (cuvars(i,j,k,1)*stepsminusone + delrho**2)*stepsinv
          cuvars(i,j,k,2) = (cuvars(i,j,k,2)*stepsminusone + delpx**2)*stepsinv
          cuvars(i,j,k,3) = (cuvars(i,j,k,3)*stepsminusone + delpy**2)*stepsinv
          cuvars(i,j,k,4) = (cuvars(i,j,k,4)*stepsminusone + delpz**2)*stepsinv
          cuvars(i,j,k,5) = (cuvars(i,j,k,5)*stepsminusone + delenergy**2)*stepsinv      

          delvelx = (delpx - primmeans(i,j,k,2)*delrho)*densitymeaninv
          delvely = (delpy - primmeans(i,j,k,3)*delrho)*densitymeaninv
          delvelz = (delpz - primmeans(i,j,k,4)*delrho)*densitymeaninv

          primvars(i,j,k,1) = cuvars(i,j,k,1)
          primvars(i,j,k,2) = (primvars(i,j,k,2)*stepsminusone + delvelx**2)*stepsinv
          primvars(i,j,k,3) = (primvars(i,j,k,3)*stepsminusone + delvely**2)*stepsinv
          primvars(i,j,k,4) = (primvars(i,j,k,4)*stepsminusone + delvelz**2)*stepsinv

         

          delg = primvars(i,j,k,2)*delpx + primvars(i,j,k,3)*delpy + primvars(i,j,k,4)*delpz

          primvars(i,j,k,nprimvars+1) = (primvars(i,j,k,nprimvars+1)*stepsminusone + delg**2)*stepsinv  !gvar

          primvars(i,j,k,nprimvars+2) = (primvars(i,j,k,nprimvars+2)*stepsminusone + delg*delenergy)*stepsinv   !kgcross
          primvars(i,j,k,nprimvars+3) = (primvars(i,j,k,nprimvars+3)*stepsminusone + delrho*delenergy)*stepsinv !krcross
          primvars(i,j,k,nprimvars+4) = (primvars(i,j,k,nprimvars+4)*stepsminusone + delrho*delg)*stepsinv      !rgcross

          primvars(i,j,k,5) = (primvars(i,j,k,5)*stepsminusone + cvinv*cvinv*densitymeaninv*densitymeaninv*(cuvars(i,j,k,5) + primvars(i,j,k,nprimvars+1) - 2*primvars(i,j,k,nprimvars+2) &
                            + qmean*(qmean*cuvars(i,j,k,1) - 2*primvars(i,j,k,nprimvars+3) + 2*primvars(i,j,k,nprimvars+4))))*stepsinv

          !deltemp = (delenergy - delg - qmean*delrho)*cvinv*densitymeaninv

          !spatialcross1(i,j,k) = (spatialcross1(i,j,k)*stepsminusone + del1*density(i,j,k))*stepsinv

          !spatialcross2(i,j,k) = spatialcross1(i,j,k) - densitymean(i,j,k)*del2

          !spatialcross1(i,j,k) = (spatialcross1(i,j,k)*stepsminusone + del1*density(i,j,k) - del2*densitymean(i,j,k))*stepsinv
          !spatialcross2(i,j,k) = (spatialcross2(i,j,k)*stepsminusone + del2*deltemp)*stepsinv

          !print *, spatialcross1(i,j,k)

        enddo
      enddo
    enddo
          
    end subroutine evaluate_corrs

end module stats_module
