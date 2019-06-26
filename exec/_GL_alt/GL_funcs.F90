module time_step_module

  use amrex_fort_module, only : amrex_real
  use common_namelist_module, only : ngc, nvars, nprimvars, cfl, prob_lo, prob_hi, n_cells
  use GL_namelist_module
  implicit none

  private

  public :: rk2_stage1, rk2_stage2, initphi, integrate, setdt, Stat_Quant, inc_phi0_Adapt, Umbrella_Adjust

contains

  subroutine rk2_stage1(lo,hi, phi, phin, rannums, integral, energy, teng, dx, dt,phi_avg,phi_sq_avg) bind(C,name="rk2_stage1")

      integer         , intent(in   ) :: lo(2),hi(2)
      real(amrex_real), intent(in   ) :: dx(2), dt
      real(amrex_real), intent(inout) :: energy, teng
      real(amrex_real), intent(inout) :: phi_avg, phi_sq_avg


      real(amrex_real), intent(inout) :: phi(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2))
      real(amrex_real), intent(inout) :: phin(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2))

      real(amrex_real), intent(in)    :: rannums(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2))
      real(amrex_real), intent(in)    :: integral

      integer :: i,j,k,l
      real(amrex_real) :: func,factor, dele

       factor = sqrt(2.d0*noise_coef/(dt*dx(1)*dx(2)))

         do  j=lo(2),hi(2)
           do  i=lo(1),hi(1)

             func =  acoef+2.d0*bcoef*phi(i,j)+3.d0*ccoef*phi(i,j)**2+4.d0*dcoef*phi(i,j)**3

             phin(i,j) = phi(i,j) + dt*diff_coef*(phi(i+1,j)-2.d0*phi(i,j)+phi(i-1,j))/dx(1)**2  &
              + dt*diff_coef*(phi(i,j+1)-2.d0*phi(i,j)+phi(i,j-1))/dx(2)**2  &
              -dt*func-dt*umbrella*integral + dt*factor * rannums(i,j)

              dele =  phi(i,j)*(acoef+bcoef*phi(i,j)+ccoef*phi(i,j)**2+dcoef*phi(i,j)**3) &
                + 0.5d0*diff_coef*((phi(i,j)-phi(i-1,j))**2+(phi(i,j)-phi(i,j-1))**2)

              energy = energy + dele

!             energy = energy + phi(i,j)*(acoef+bcoef*phi(i,j)+ccoef*phi(i,j)**2+dcoef*phi(i,j)**3) &
!               + diff_coef/2.d0*((phi(i,j)-phi(i-1,j))**2+(phi(i,j)-phi(i,j-1))**2)

              teng = teng + dele + 0.5d0*umbrella*integral

          enddo
        enddo

            energy = energy *dx(1)*dx(2)
            teng = teng *dx(1)*dx(2)

         do  j=lo(2),hi(2)
           do  i=lo(1),hi(1)

               phi(i,j)= phin(i,j)

          enddo
        enddo
        phi_avg=0.0d0
        phi_sq_avg=0.0d0
      call Stat_Quant(lo,hi, phi,phi_avg,phi_sq_avg)
  end subroutine rk2_stage1

  subroutine rk2_stage2(lo,hi, phi, phin, rannums, integral, dx, dt) bind(C,name="rk2_stage2")

      integer         , intent(in   ) :: lo(2),hi(2)
      real(amrex_real), intent(in   ) :: dx(2), dt

      real(amrex_real), intent(inout) :: phi(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2))
      real(amrex_real), intent(inout) :: phin(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2))

      real(amrex_real), intent(in)    :: rannums(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2))
      real(amrex_real), intent(in)    :: integral

      integer :: i,j,k,l
      real(amrex_real) :: func,factor

       factor = sqrt(2.d0*noise_coef/(dt*dx(1)*dx(2)))

         do  j=lo(2),hi(2)
           do  i=lo(1),hi(1)


             func =  acoef+2.d0*bcoef*phi(i,j)+3.d0*ccoef*phi(i,j)**2+4.d0*dcoef*phi(i,j)**3

             phin(i,j) = phi(i,j) + dt*diff_coef*(phi(i+1,j)-2.d0*phi(i,j)+phi(i-1,j))/dx(1)**2  &
              + dt*diff_coef*(phi(i,j+1)-2.d0*phi(i,j)+phi(i,j-1))/dx(2)**2  &
              -dt*func-dt*umbrella*integral + dt*factor * rannums(i,j)

          enddo
        enddo

  end subroutine rk2_stage2

  subroutine initphi(lo,hi, phi, dx) bind(C,name="initphi")

      integer         , intent(in   ) :: lo(2),hi(2)
      real(amrex_real), intent(in   ) :: dx(2)

      real(amrex_real), intent(inout) :: phi(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2))

      real(amrex_real) :: xloc,yloc, xcen, ycen
      integer :: i,j

         xcen = 0.5d0*(prob_lo(1) + prob_hi(1))
         ycen = 0.5d0*(prob_lo(2) + prob_hi(2))

         do  j=lo(2),hi(2)
           do  i=lo(1),hi(1)

              xloc = dx(1)*dfloat(i-1)+prob_lo(1)
              yloc = dx(2)*dfloat(j-1)+prob_lo(2)

           !  phi(i,j) =  (sin(2.d0*pi*xloc)*sin(2.d0*pi*yloc)) **2
           !  if( (xloc-.5d0)**2 + (yloc-.5d0)**2 .lt. rad**2)then
              if( (xloc-xcen)**2 + (yloc-xcen)**2 .lt. rad**2)then
                 phi(i,j) = 1.d0
              else
                 phi(i,j) = 0.d0
              endif


          enddo
        enddo

  end subroutine initphi

  subroutine integrate(lo,hi, phi, dx, integral) bind(C,name="integrate")

      integer         , intent(in   ) :: lo(2),hi(2)
      real(amrex_real), intent(in   ) :: dx(2)
      real(amrex_real), intent(out  ) :: integral

      real(amrex_real), intent(inout) :: phi(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2))

      real(amrex_real) :: xloc,yloc
      integer :: i,j

         integral = 0.d0

         do  j=lo(2),hi(2)
           do  i=lo(1),hi(1)

              integral = integral + phi(i,j) - phi0

          enddo
        enddo

        integral = integral*dx(1)*dx(2)

  end subroutine integrate





  subroutine inc_phi0 ( step ) bind(C,name="inc_phi0")

  integer, intent(in) :: step

  if(mod(step,n_inc_phi) .eq. 0)then

     phi0 = phi0 + phi_inc
     write(6,*)" phi0 changed to ",phi0 ," at step ", step

  endif

     phi0 = min(phi0, 1.d0)

  end subroutine inc_phi0



  subroutine inc_phi0_Adapt ( Expec,MAD,r1,Shift_Flag ) bind(C,name="inc_phi0_Adapt")

    real(amrex_real), intent(inout  ) :: Expec,MAD,r1
    integer, intent(inout ) :: Shift_Flag
    write(6,*)" phi0 is",phi0 
    if(Shift_Flag .NE. 1) then
      phi0 =(Expec+r1*MAD)
    else 
      phi0=phi0+r1*MAD
    end if
    write(6,*)" phi0 changed to ",phi0 


  end subroutine inc_phi0_Adapt

  subroutine Umbrella_Adjust ( sucessful_iter,alpha,umbrella_size,sucessful_iter_prev) bind(C,name="Umbrella_Adjust")

    real(amrex_real), intent(in  ) :: alpha
    integer, intent(inout ) :: sucessful_iter,umbrella_size,sucessful_iter_prev

    if (umbrella <=50) then
      umbrella_size=1
    end if
    !! check previous step, current step, and umbrella_size. This section must account for all positbiltities
    write(6,*)"umbrella is",umbrella 

    if (sucessful_iter_prev==1 .AND. sucessful_iter==1 .AND. umbrella_size==0) then 
         umbrella=umbrella/alpha
    else if(sucessful_iter_prev==1 .AND. sucessful_iter==1 .AND. umbrella_size==1) then
         umbrella=umbrella
    else if(sucessful_iter_prev==0 .AND. sucessful_iter==1) then
         umbrella=umbrella
    else if(sucessful_iter_prev==0 .AND. sucessful_iter==0) then
         umbrella=umbrella*alpha
    else if(sucessful_iter_prev==1 .AND. sucessful_iter==0) then
         umbrella=umbrella*alpha
         umbrella_size=1
    end if

    write(6,*)" umbrella changed to ",umbrella 

  end subroutine Umbrella_Adjust


  subroutine setdt (  dx, dt ) bind(C,name="setdt")

  double precision::  dt, dx(2)

  if(cfl.gt.0.d0)then

      dt = .25d0*cfl*min(dx(1)**2, dx(2)**2)/diff_coef
    
  endif

  end subroutine setdt




  subroutine Stat_Quant(lo,hi, phi,phi_avg,phi_sq_avg) bind(C,name="Stat_Quant")

    integer         , intent(in   ) :: lo(2),hi(2)
    real(amrex_real), intent(in   ) :: phi(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2))
    real(amrex_real), intent(inout  ) :: phi_avg
    real(amrex_real), intent(inout  ) :: phi_sq_avg
    integer :: i,j
  
      phi_avg = 0.d0
       do  j=lo(2),hi(2)
         do  i=lo(1),hi(1)
  
          phi_avg = phi_avg + phi(i,j)
  
        enddo
      enddo
      phi_avg = phi_avg/(n_cells(2)*n_cells(1))
  
  
      phi_sq_avg = 0.d0
      do  j=lo(2),hi(2)
        do  i=lo(1),hi(1)
  
          phi_sq_avg = phi_sq_avg + phi(i,j)*phi(i,j)
  
       enddo
     enddo
  
     phi_sq_avg = phi_sq_avg/(n_cells(2)*n_cells(1))
  
  end subroutine Stat_Quant




end module time_step_module



