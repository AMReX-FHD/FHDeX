module stats_module

  use amrex_fort_module, only : amrex_real
  use common_namelist_module, only : ngc, nvars, nprimvars, nspecies, cell_depth, k_b, &
       n_cells, hcv, cross_cell
  use conv_module
  implicit none

  private

  public :: evaluate_corrs, multifab_yzav

contains

  subroutine evaluate_corrs(lo, hi, cu, cumeans, cuvars, prim, primmeans, primvars, spatialcross, steps, miscstats, miscvals) bind(c,name='evaluate_corrs')

      implicit none

      integer,          intent(in   ) :: steps, lo(3), hi(3)
      double precision, intent(inout   ) :: miscvals(20)
      double precision, intent(in   ) :: cu       (lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3), nvars)
      double precision, intent(in   ) :: cumeans  (lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3), nvars)
      double precision, intent(inout) :: cuvars   (lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3), nvars)
      double precision, intent(in   ) :: prim     (lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3), nprimvars)
      double precision, intent(in   ) :: primmeans(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3), nprimvars)
      double precision, intent(inout) :: primvars (lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3), nprimvars + 5)
      double precision, intent(inout) :: spatialcross(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3), 6)
      double precision, intent(inout) :: miscstats(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3), 10)

      double precision slices(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3), 20)

      integer i,j,k,l, counter
      double precision stepsminusone, stepsinv, cv, cvinv,delg, qmean, delpx, delpy, delpz, delrho, delvelx, delvely, delvelz, delenergy, densitymeaninv
      double precision deltemp, delpdelrho, delrhoS, delrhoSstar,delpxS, delpyS, delpzS, cvinvS, cvinvSstar, delES, delESstar, qmeanS, qmeanSstar
      double precision delgS, delgSstar, delpxSstar, delpySstar, delpzSstar, deltempS, deltempSstar, densitymeaninvS, densitymeaninvSstar

      stepsminusone = steps - 1.d0
      stepsinv = 1.d0/steps


      !print *, "STATS!"
      counter = 0

     do k = lo(3), hi(3)
        do j = lo(2), hi(2)
 
        counter = counter +1
      enddo
    enddo


     slices = 0

 
      do k = lo(3), hi(3)
        do j = lo(2), hi(2)
          do i = lo(1), hi(1)

            slices(i,lo(2),lo(3),1) = slices(i,lo(2),lo(3),1) + cu(i,j,k,1) !rho instant slices
            slices(i,lo(2),lo(3),2) = slices(i,lo(2),lo(3),2) + cumeans(i,j,k,1) !rho mean slices
            slices(i,lo(2),lo(3),3) = slices(i,lo(2),lo(3),3) + cu(i,j,k,5) !energy instant slices
            slices(i,lo(2),lo(3),4) = slices(i,lo(2),lo(3),4) + cumeans(i,j,k,5) !energy mean slices

            slices(i,lo(2),lo(3),5) = slices(i,lo(2),lo(3),5) + cu(i,j,k,2) !x momentum instant slices
            slices(i,lo(2),lo(3),6) = slices(i,lo(2),lo(3),6) + cumeans(i,j,k,2) !x momentum mean slices

            slices(i,lo(2),lo(3),7) = slices(i,lo(2),lo(3),7) + cu(i,j,k,3) !y momentum instant slices
            slices(i,lo(2),lo(3),8) = slices(i,lo(2),lo(3),8) + cumeans(i,j,k,3) !y momentum mean slices


            slices(i,lo(2),lo(3),9) = slices(i,lo(2),lo(3),9) + cu(i,j,k,4) !z momentum instant slices
            slices(i,lo(2),lo(3),10) = slices(i,lo(2),lo(3),10) + cumeans(i,j,k,4) !z momentum mean slices

            slices(i,lo(2),lo(3),11) = slices(i,lo(2),lo(3),11) + prim(i,j,k,2) !x vel instant slices
            slices(i,lo(2),lo(3),12) = slices(i,lo(2),lo(3),12) + primmeans(i,j,k,2) !x vel mean slices

            slices(i,lo(2),lo(3),13) = slices(i,lo(2),lo(3),13) + prim(i,j,k,3) !y vel instant slices
            slices(i,lo(2),lo(3),14) = slices(i,lo(2),lo(3),14) + primmeans(i,j,k,3) !y vel mean slices


            slices(i,lo(2),lo(3),15) = slices(i,lo(2),lo(3),15) + prim(i,j,k,4) !z vel instant slices
            slices(i,lo(2),lo(3),16) = slices(i,lo(2),lo(3),16) + primmeans(i,j,k,4) !z vel mean slices

            cv = 0
            do l=1,nspecies
              cv = cv + hcv(l)*cumeans(i,j,k,5+l)/cumeans(i,j,k,1)
            enddo

            slices(i,lo(2),lo(3),17) = slices(i,lo(2),lo(3),17) + cv !cv mean slices

            slices(i,lo(2),lo(3),18) = slices(i,lo(2),lo(3),18) + prim(i,j,k,5) !temperature instant slices
            slices(i,lo(2),lo(3),19) = slices(i,lo(2),lo(3),19) + primmeans(i,j,k,5) !temperature mean slices

        enddo
      enddo
    enddo

    slices = slices/counter

      do k = lo(3), hi(3)
        do j = lo(2), hi(2)
          do i = lo(1), hi(1)

            slices(i,j,k,1) = slices(i,lo(2),lo(3),1) 
            slices(i,j,k,2) = slices(i,lo(2),lo(3),2) 
            slices(i,j,k,3) = slices(i,lo(2),lo(3),3)
            slices(i,j,k,4) = slices(i,lo(2),lo(3),4) 
            slices(i,j,k,5) = slices(i,lo(2),lo(3),5) 
            slices(i,j,k,6) = slices(i,lo(2),lo(3),6)
            slices(i,j,k,7) = slices(i,lo(2),lo(3),7)
            slices(i,j,k,8) = slices(i,lo(2),lo(3),8) 
            slices(i,j,k,9) = slices(i,lo(2),lo(3),9) 
            slices(i,j,k,10) = slices(i,lo(2),lo(3),10)
            slices(i,j,k,11) = slices(i,lo(2),lo(3),11) 
            slices(i,j,k,12) = slices(i,lo(2),lo(3),12) 
            slices(i,j,k,13) = slices(i,lo(2),lo(3),13)
            slices(i,j,k,14) = slices(i,lo(2),lo(3),14) 
            slices(i,j,k,15) = slices(i,lo(2),lo(3),15) 
            slices(i,j,k,16) = slices(i,lo(2),lo(3),16)
            slices(i,j,k,17) = slices(i,lo(2),lo(3),17)
            slices(i,j,k,18) = slices(i,lo(2),lo(3),18) 
            slices(i,j,k,19) = slices(i,lo(2),lo(3),19) 
            slices(i,j,k,20) = slices(i,lo(2),lo(3),20) 

        enddo
      enddo
    enddo

      do k = lo(3), hi(3)
        do j = lo(2), hi(2)
          do i = lo(1), hi(1)

          cv = 0.d0
          do l=1,nspecies
            cv = cv + hcv(l)*cumeans(i,j,k,5+l)/cumeans(i,j,k,1)
          enddo

          qmean = cv*primmeans(i,j,k,5)-0.5d0*(primmeans(i,j,k,2)**2 + primmeans(i,j,k,3)**2 + primmeans(i,j,k,4)**2)
          cvinv = 1.0d0/cv
          cvinvS = 1.0d0/slices(i,j,k,17)
          cvinvSstar = 1.0d0/miscvals(13)

          !Vars
          qmean = cv*primmeans(i,j,k,5)-0.5*(primmeans(i,j,k,2)**2 + primmeans(i,j,k,3)**2 + primmeans(i,j,k,4)**2)
          qmeanS = slices(i,j,k,17)*slices(i,j,k,19)-0.5*(slices(i,j,k,12)**2 + slices(i,j,k,14)**2 + slices(i,j,k,16)**2)
          qmeanSstar =  miscvals(13)*slices(i,lo(2),lo(3),19)-0.5*(miscvals(3)**2 + miscvals(15)**2 + miscvals(16)**2)

          densitymeaninv = 1.0/cumeans(i,j,k,1)
          densitymeaninvS = 1.0/slices(i,j,k,2)
          densitymeaninvSstar = 1.0/miscvals(4)

          delrho = cu(i,j,k,1) - cumeans(i,j,k,1)
          delpx = cu(i,j,k,2) - cumeans(i,j,k,2)
          delpy = cu(i,j,k,3) - cumeans(i,j,k,3)
          delpz = cu(i,j,k,4) - cumeans(i,j,k,4)
          delenergy = cu(i,j,k,5) - cumeans(i,j,k,5)

          delrhoS = slices(i,j,k,1) - slices(i,j,k,2) !rho(x) - <rho(x)>, sliced
          delES = slices(i,j,k,4) - slices(i,j,k,3) !E(x) - <E(x)>, sliced
          delpxS = slices(i,j,k,6) - slices(i,j,k,5)
          delpyS = slices(i,j,k,8) - slices(i,j,k,7)
          delpzS = slices(i,j,k,10) - slices(i,j,k,9)

          delrhoSstar = miscvals(5) - miscvals(4)
          delESstar = miscvals(7) - miscvals(8)
          delpxSstar = miscvals(2) - miscvals(1)
          delpySstar = miscvals(9) - miscvals(10)
          delpzSstar = miscvals(11) - miscvals(12)


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

          delgS = slices(i,j,k,12)*delpxS + slices(i,j,k,14)*delpyS + slices(i,j,k,16)*delpzS

          delgSstar = miscvals(3)*delpxSstar + miscvals(15)*delpySstar + miscvals(16)*delpzSstar

          primvars(i,j,k,nprimvars+1) = (primvars(i,j,k,nprimvars+1)*stepsminusone + delg**2)*stepsinv  !gvar

          primvars(i,j,k,nprimvars+2) = (primvars(i,j,k,nprimvars+2)*stepsminusone + delg*delenergy)*stepsinv   !kgcross
          primvars(i,j,k,nprimvars+3) = (primvars(i,j,k,nprimvars+3)*stepsminusone + delrho*delenergy)*stepsinv !krcross
          primvars(i,j,k,nprimvars+4) = (primvars(i,j,k,nprimvars+4)*stepsminusone + delrho*delg)*stepsinv      !rgcross

          primvars(i,j,k,5) = (primvars(i,j,k,5)*stepsminusone + cvinv*cvinv*densitymeaninv*densitymeaninv*(cuvars(i,j,k,5) + primvars(i,j,k,nprimvars+1) - 2*primvars(i,j,k,nprimvars+2) &
                            + qmean*(qmean*cuvars(i,j,k,1) - 2*primvars(i,j,k,nprimvars+3) + 2*primvars(i,j,k,nprimvars+4))))*stepsinv

          deltemp = (delenergy - delg - qmean*delrho)*cvinv*densitymeaninv

          deltempS = (delES - delgS - qmeanS*delrhoS)*cvinvS*densitymeaninvS

          deltempSstar = (delESstar - delgSstar - qmeanSstar*delrhoSstar)*cvinvSstar*densitymeaninvSstar

          miscstats(i,j,k,1) = (miscstats(i,j,k,1)*stepsminusone + miscvals(2)*slices(i,j,k,1))*stepsinv  ! <p(x*)rho(x)>, sliced

          delpdelrho = miscstats(i,j,k,1) - miscvals(1)*cumeans(i,j,k,1) !<p(x*)rho(x)> - <p(x*)><rho(x)>, sliced


          miscstats(i,j,k,2) = (miscstats(i,j,k,2)*stepsminusone + delrhoS*delrhoSstar)*stepsinv  ! <(rho(x*)-<rho(x*)>)(rho(x)-<rho(x)>)>, sliced

          miscstats(i,j,k,3) = (miscstats(i,j,k,3)*stepsminusone + miscvals(17)*slices(i,j,k,18))*stepsinv  ! <(T(x*)T(x))>
          miscstats(i,j,k,4) = (miscstats(i,j,k,4)*stepsminusone + miscvals(17)*slices(i,j,k,1))*stepsinv  ! <(T(x*)rho(x))>

          !print *, slices(i,j,k,1)

!          spatialcross(i,j,k,1) = (spatialcross(i,j,k,1)*stepsminusone + delES*delrhoSstar)*stepsinv
!          spatialcross(i,j,k,2) = (spatialcross(i,j,k,2)*stepsminusone + delES*delESstar))*stepsinv
!          spatialcross(i,j,k,3) = (spatialcross(i,j,k,3)*stepsminusone + delES*delESstar)*stepsinv

          spatialcross(i,j,k,1) = miscvals(14)
          spatialcross(i,j,k,2) = slices(i,j,k,19)
          spatialcross(i,j,k,3) = miscstats(i,j,k,3)

          spatialcross(i,j,k,4) = miscstats(i,j,k,3) - slices(i,j,k,19)*miscvals(14)
          spatialcross(i,j,k,5) = miscstats(i,j,k,4) - slices(i,j,k,2)*miscvals(14)
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
