module stats_module

  use amrex_fort_module, only : amrex_real
  use common_namelist_module, only : ngc, nvars, nprimvars, nspecies, cell_depth, k_b, &
       n_cells, hcv, cross_cell
  use conv_module
  implicit none

  private

  public :: evaluate_means, evaluate_corrs, multifab_yzav

contains

  subroutine evaluate_means(lo, hi, cu, cumeans, prim, primmeans, steps, miscstats, miscvals, totalmass) bind(c,name='evaluate_means')
 
      implicit none

      integer,          intent(in   ) :: steps, lo(3), hi(3)
      double precision, intent(inout   ) :: totalmass, miscvals(10)
      double precision, intent(in   ) :: cu       (lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3), nvars)
      double precision, intent(inout) :: cumeans  (lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3), nvars)
      double precision, intent(in   ) :: prim     (lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3), nprimvars)
      double precision, intent(inout) :: primmeans(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3), nprimvars)
      double precision, intent(inout) :: miscstats(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3), 10)

      integer i,j,k,l, cells, ti, jc, kc, counter
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

!              mixedstats(i,j,k,1) = cumeans(i,j,k,l)

!            do l=1,nprimvars
!              primmeans(i,j,k,l) = (primmeans(i,j,k,l)*stepsminusone + prim(i,j,k,l))*stepsinv
!            enddo

                !print *, "means: ", cumeans(i,j,k,l), stepsminusone


            fracvec = cumeans(i,j,k,6:nvars)/cumeans(i,j,k,1)
            massvec = cumeans(i,j,k,6:nvars)

            densitymeaninv = 1.0/cumeans(i,j,k,1)

            primmeans(i,j,k,1) = cumeans(i,j,k,1)
            primmeans(i,j,k,2) = cumeans(i,j,k,2)*densitymeaninv
            primmeans(i,j,k,3) = cumeans(i,j,k,3)*densitymeaninv
            primmeans(i,j,k,4) = cumeans(i,j,k,4)*densitymeaninv

            call get_temperature(cumeans(i,j,k,5), massvec, primmeans(i,j,k,5))
            call get_pressure_gas(primmeans(i,j,k,6), fracvec, cumeans(i,j,k,1),cumeans(i,j,k,5))

            totalmass = totalmass + cu(i,j,k,1)


          enddo
        enddo
      enddo

    do i = 1,10
      miscVals(i) = 0
    enddo

    counter = 0

    if((cross_cell .ge. lo(1)) .and. (cross_cell .le. hi(1))) then
      
      do k = lo(3), hi(3)
        do j = lo(2), hi(2)
          
          miscvals(1) = miscvals(1) + cumeans(cross_cell,j,k,2) !slice average of mean x momentum
          miscvals(2) = miscvals(2) + cu(cross_cell,j,k,2) !slice average of instant x momentum
          miscvals(3) = miscvals(3) + primmeans(cross_cell,j,k,2) !slice average of mean x velocity
          miscvals(4) = miscvals(4) + cumeans(cross_cell,j,k,1) !slice average of mean rho
          miscvals(5) = miscvals(5) + cu(cross_cell,j,k,1) !slice average of instant rho
          miscvals(6) = miscvals(6) + prim(cross_cell,j,k,2) !slice average of instant x velocity

          counter = counter + 1

         enddo
       enddo

       miscvals(1) = miscvals(1)/counter !slice average of mean x momentum
       miscvals(2) = miscvals(2)/counter !slice average of instant x momentum
       miscvals(3) = miscvals(3)/counter !slice average of mean x velocity
       miscvals(4) = miscvals(4)/counter !slice average of mean rho
       miscvals(5) = miscvals(5)/counter !slice average of instant rho
       miscvals(6) = miscvals(6)/counter !slice average of instant x velocity

     endif
          
  end subroutine evaluate_means

  subroutine evaluate_corrs(lo, hi, cu, cumeans, cuvars, prim, primmeans, primvars, spatialcross, steps, miscstats, miscvals) bind(c,name='evaluate_corrs')

      implicit none

      integer,          intent(in   ) :: steps, lo(3), hi(3)
      double precision, intent(inout   ) :: miscvals(10)
      double precision, intent(in   ) :: cu       (lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3), nvars)
      double precision, intent(in   ) :: cumeans  (lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3), nvars)
      double precision, intent(inout) :: cuvars   (lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3), nvars)
      double precision, intent(in   ) :: prim     (lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3), nprimvars)
      double precision, intent(in   ) :: primmeans(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3), nprimvars)
      double precision, intent(inout) :: primvars (lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3), nprimvars + 5)
      double precision, intent(inout) :: spatialcross(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3), 6)
      double precision, intent(inout) :: miscstats(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3), 10)

      double precision slices(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3), 10)

      integer i,j,k,l, counter
      double precision stepsminusone, stepsinv, cv, cvinv,delg, qmean, delpx, delpy, delpz, delrho, delvelx, delvely, delvelz, delenergy, densitymeaninv, deltemp, delpdelrho, delrhoS, delrhostarS

      stepsminusone = steps - 1
      stepsinv = 1d0/steps

      counter = 0

     do k = lo(3), hi(3)
        do j = lo(2), hi(2)
 
        counter = counter +1
      enddo
    enddo


      do l = 1,10
      do k = lo(3), hi(3)
        do j = lo(2), hi(2)
          do i = lo(1), hi(1)

            slices(i,j,k,l) = 0

        enddo
      enddo
    enddo
    enddo

 
      do k = lo(3), hi(3)
        do j = lo(2), hi(2)
          do i = lo(1), hi(1)

            slices(i,1,1,1) = slices(i,1,1,1) + cu(i,j,k,1) !rho instant slices
            slices(i,1,1,2) = slices(i,1,1,2) + cumeans(i,j,k,1) !rho mean slices

        enddo
      enddo
    enddo

    do i = lo(1), hi(1)

        slices(i,1,1,1) = slices(i,1,1,1)/counter
        slices(i,1,1,2) = slices(i,1,1,2)/counter

    enddo

      do k = lo(3), hi(3)
        do j = lo(2), hi(2)
          do i = lo(1), hi(1)

            slices(i,j,k,1) = slices(i,1,1,1) !rho instant slices
            slices(i,j,k,2) = slices(i,1,1,2) !rho mean slices

        enddo
      enddo
    enddo

      do k = lo(3), hi(3)
        do j = lo(2), hi(2)
          do i = lo(1), hi(1)

          cv = 0
          do l=1,nspecies
            cv = cv + hcv(l)*cumeans(i,j,k,5+l)/cumeans(i,j,k,1)
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

          deltemp = (delenergy - delg - qmean*delrho)*cvinv*densitymeaninv

          miscstats(i,j,k,1) = (miscstats(i,j,k,1)*stepsminusone + miscvals(2)*slices(i,j,k,1))*stepsinv  ! <p(x*)rho(x)>, sliced

          delpdelrho = miscstats(i,j,k,1) - miscvals(1)*cumeans(i,j,k,1) !<p(x*)rho(x)> - <p(x*)><rho(x)>, sliced

          delrhoS = slices(i,j,k,1) - slices(i,j,k,2) !rho(x) - <rho(x)>, sliced
          delrhostarS = miscvals(5) - miscvals(4) !rho(x*) - <rho(x*)>, sliced

          miscstats(i,j,k,2) = (miscstats(i,j,k,2)*stepsminusone + delrhoS*delrhostarS)*stepsinv  ! <(rho(x*)-<rho(x*)>)(rho(x)-<rho(x)>)>, sliced

          !print *, (n_cells(2))*(k) + (j+1), k, j

!          spatialcross(i,j,k,1) = (spatialcross(i,j,k,1)*stepsminusone + delrho*delholder1((n_cells(2))*(k) + (j+1)))*stepsinv
!          spatialcross(i,j,k,2) = (spatialcross(i,j,k,2)*stepsminusone + delenergy*delholder2((n_cells(2))*(k) + (j+1)))*stepsinv
!          spatialcross(i,j,k,3) = (spatialcross(i,j,k,3)*stepsminusone + delrho*delholder3((n_cells(2))*(k) + (j+1)))*stepsinv
!          spatialcross(i,j,k,4) = (spatialcross(i,j,k,4)*stepsminusone + deltemp*delholder4((n_cells(2))*(k) + (j+1)))*stepsinv
          spatialcross(i,j,k,5) = (spatialcross(i,j,k,5)*stepsminusone + delrhoS*(miscvals(6)-miscvals(3)))*stepsinv  !<(u(x*)-<u(x*)>)(rho(x)-<rho(x)>)> !sliced
          spatialcross(i,j,k,6) = (delpdelrho - miscvals(3)*miscstats(i,j,k,2))/miscvals(4)

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
