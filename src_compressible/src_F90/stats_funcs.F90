module stats_module

  use amrex_fort_module, only : amrex_real
  use common_namelist_module, only : ngc, nvars, nprimvars, nspecies, cell_depth, k_b, bc_lo, bc_hi, n_cells, hcv, cross_cell
  use conv_module
  implicit none

  private

  public :: evaluate_means, evaluate_corrs

contains

  subroutine evaluate_means(lo, hi, cu, cumeans, prim, primmeans, steps, totalmass) bind(c,name='evaluate_means')
 
      implicit none

      integer,          intent(in      ) :: steps, lo(3), hi(3)
      double precision, intent(inout   ) :: totalmass

      double precision, intent(inout   ) :: cu(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3), nvars)
      double precision, intent(inout   ) :: cumeans(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3), nvars)

      double precision, intent(inout   ) :: prim(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3), nprimvars)
      double precision, intent(inout   ) :: primmeans(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3), nprimvars)

      !double precision fac1, fac2, fac3, test, pairfrac
      integer i,j,k,l, cells
      double precision stepsminusone, stepsinv, densitymeaninv, fracvec(nspecies), massvec(nspecies)
 
      stepsminusone = steps - 1
      stepsinv = 1d0/steps

      totalmass = 0  

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

            totalmass = totalmass + cu(i,j,k,1)


          enddo
        enddo
      enddo
          
    end subroutine evaluate_means

  subroutine evaluate_corrs(lo, hi, cu, cumeans, cuvars, prim, primmeans, primvars, steps) bind(c,name='evaluate_corrs')

      !use iso_c_binding, only: c_ptr, c_int, c_f_pointer

      implicit none

      integer,          intent(in      ) :: steps, lo(3), hi(3)

      double precision, intent(inout   ) :: cu(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3), nvars)
      double precision, intent(inout   ) :: cumeans(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3), nvars)
      double precision, intent(inout   ) :: cuvars(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3), nvars)

      double precision, intent(inout   ) :: prim(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3), nprimvars)
      double precision, intent(inout   ) :: primmeans(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3), nprimvars)
      double precision, intent(inout   ) :: primvars(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3), nprimvars + 5)

      integer i,j,k,l
      double precision stepsminusone, stepsinv, cv, cvinv,delg, qmean, delpx, delpy, delpz, delrho, delvelx, delvely, delvelz, delenergy, densitymeaninv, deltemp

      stepsminusone = steps - 1
      stepsinv = 1d0/steps

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

      
          delg = primmeans(i,j,k,2)*delpx + primmeans(i,j,k,3)*delpy + primmeans(i,j,k,4)*delpz

          primvars(i,j,k,nprimvars+1) = (primvars(i,j,k,nprimvars+1)*stepsminusone + delg**2)*stepsinv  !gvar

          primvars(i,j,k,nprimvars+2) = (primvars(i,j,k,nprimvars+2)*stepsminusone + delg*delenergy)*stepsinv   !kgcross
          primvars(i,j,k,nprimvars+3) = (primvars(i,j,k,nprimvars+3)*stepsminusone + delrho*delenergy)*stepsinv !krcross
          primvars(i,j,k,nprimvars+4) = (primvars(i,j,k,nprimvars+4)*stepsminusone + delrho*delg)*stepsinv      !rgcross

          primvars(i,j,k,5) = (primvars(i,j,k,5)*stepsminusone + cvinv*cvinv*densitymeaninv*densitymeaninv*(cuvars(i,j,k,5) + primvars(i,j,k,nprimvars+1) - 2*primvars(i,j,k,nprimvars+2) &
                            + qmean*(qmean*cuvars(i,j,k,1) - 2*primvars(i,j,k,nprimvars+3) + 2*primvars(i,j,k,nprimvars+4))))*stepsinv

          ! deltemp = (delenergy - delg - qmean*delrho)*cvinv*densitymeaninv

        enddo
      enddo
    enddo
          
    end subroutine evaluate_corrs

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
