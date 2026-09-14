module time_step_module

  use amrex_fort_module, only : amrex_real
  use common_namelist_module, only : ngc, nvars, nprimvars, cfl, prob_lo, prob_hi, n_cells
  use GL_namelist_module
  implicit none

  private

  public :: inc_phi0_Adapt, Umbrella_Adjust, Param_Output,umbrella_reset,set_inputs,fixed_inc_phi0
  public :: Stat_Quant_1D, integrate_1D,  initphi_1D, rk2_stage1_1D,setdt_1D

contains

  subroutine set_inputs(Plot_Skip_out,Number_of_Samples_out,Equil_out,alpha_out,r1_out,r2_out, adaptive_param,reverse_param) bind(C,name="set_inputs")
    integer,    intent(inout)   :: Plot_Skip_out, Number_of_Samples_out,Equil_out,adaptive_param,reverse_param
    real(amrex_real), intent(inout) :: alpha_out,r1_out,r2_out
    ! This subroutine takes inputs read in the "GL_namelist.F90" file, but not in the common name list, and sets outputs
    ! of the subroutine to these values. This is primarily used to send input file quantities to the c++ functions

    ! NOTE: Explanations of these quantities are in the input file
    Plot_Skip_out=Plot_Skip
    Number_of_Samples_out=Number_of_Samples
    Equil_out=Equil
    alpha_out=alpha
    r1_out=r1
    r2_out=r2
    adaptive_param=adaptive
    reverse_param=Reverse
  end subroutine set_inputs



  subroutine inc_phi0 ( step ) bind(C,name="inc_phi0")
    ! This subroutine increments the umbrella center by a value specified in the input value (phi_inc)
    ! This subroutine increments \phi whenever a muliple of "n_inc_phi"(an input parameter is reached).
    !INPUT:
    ! step -- an integer corresponsing to the index of the time step


    ! Update: Since 06/2019, both the adaptive and fixed increment approaches do not use this function anymore
  integer, intent(in) :: step

  if(mod(step,n_inc_phi) .eq. 0)then

     phi0 = phi0 + phi_inc
     write(6,*)"phi0 changed to ",phi0 ," at step ", step

            !  open(7,file="Console_output_Fortran.txt",status="old",position="append",action="write")
            !  write(7,*)"phi0 changed to ",phi0 ," at step ", step
            !  close(7)

  endif

     phi0 = min(phi0, 1.d0)

  end subroutine inc_phi0



  subroutine fixed_inc_phi0 ( forward_input ) bind(C,name="fixed_inc_phi0")
    ! This subroutines increments \phi by "phi_inc" whenever it is called.
    ! This subroutine takes an integer 0 or 1 input that determines whether or not the steps increments
    ! by adding or subtracting. This step is then written to a file and the console ouput.

    ! INPUT:
    ! forward_input -- flag that indicates direction to increment phi_0
    ! forward_input==1 means    phi0 = phi0 + phi_inc (forward direction)
    ! forward_input==0 means    phi0 = phi0 - phi_inc (backward direction)
    integer, intent(in) :: forward_input

    if(forward_input .eq. 1)then
       phi0 = phi0 + phi_inc      !increment in forward direction
       write(6,*)"phi0 changed to ",phi0

              !  open(7,file="Console_output_Fortran.txt",status="old",position="append",action="write")
              !  write(7,*)"phi0 changed to ",phi0
              !  close(7)
    endif

    if(forward_input .eq. 0)then
      phi0 = phi0 - phi_inc     !increment in backward direction
      write(6,*)"phi0 changed to ",phi0

              ! open(7,file="Console_output_Fortran.txt",status="old",position="append",action="write")
              ! write(7,*)"phi0 changed to ",phi0
              ! close(7)
   endif

    end subroutine fixed_inc_phi0



  subroutine inc_phi0_Adapt ( Expec,MAD,r1,Shift_Flag ) bind(C,name="inc_phi0_Adapt")
    !This subroutine computes and sets \phi_0 according to criteria determined by the adaptive algorithm and writes this to file and console.

    ! The subroutine gives the next \phi_0 using one of two methods

    ! First approach: The next \phi_0 is determined as the spatial average of the current \phi field plus some
    ! parameter "r1" times some measure of the spread of the data "MAD" (The median average deviation has been used here).
    ! This is usually the method used before high energy barriers are encountered.

    ! Second approach: The next \phi_0 is determined as the previous \phi_0 plus some
    ! parameter "r1" times some measure of the spread of the data "MAD" (The median average deviation has been used here).
    ! This is usually the method used when high energy barriers are encountered and spatial averages of \phi deviate from \phi_0

    ! INPUT:
    ! Expec- A measure of center of an umbrella window (e.g average or median). Set as average at present
    ! MAD - A measure of the spread of the data in umbrella window (e.g standard deviation, median absolute deviation, etc). Set as  median absolute deviation.
    ! r1 - A scaling factor for the spread of the data used when determining the "next" phi_0. See overleaf for details.
    ! Shift_Flag- A flag that specifies which "approach" is being used

    real(amrex_real), intent(inout  ) :: Expec,MAD,r1
    integer, intent(inout ) :: Shift_Flag
            ! open(7,file="Console_output_Fortran.txt",status="old",position="append",action="write")

    write(6,*)" phi0 is",phi0
            ! write(7,*)" phi0 is",phi0

    if(Shift_Flag .NE. 1) then    ! First approach
      phi0 =(Expec+r1*MAD)
    else
      phi0=phi0+r1*MAD     ! Second approach
    end if
    write(6,*)" phi0 changed to ",phi0
            ! write(7,*)" phi0 changed to ",phi0
            ! close(7)

  end subroutine inc_phi0_Adapt





  subroutine Umbrella_Adjust ( sucessful_iter,alpha,umbrella_size,sucessful_iter_prev) bind(C,name="Umbrella_Adjust")
    ! Subroutine for adjusting the umbrella for some common scenarious that occur.
    ! Precise details of the adaptive algorithm are in the overleaf notes.

    ! INPUT:
    ! sucessful_iter -- Flag that states whether the current step had sufficient overlap(i.e a successful step). 1= yes, 0=no.
    ! alpha -- Parameter for scaling the spring constant
    ! umbrella_size -- Flag ehose value depends on whether the spring constant is too large or small. (1=make no smaller, 2=make no larger)
    ! sucessful_iter_prev -- Flag that states whether the previous step was successful. 1= yes, 0=no.

    !NOTE: Sometimes, the method fails if the minimum spring constant is too weak. Some trail and error is needed to determine
    ! good choices of the input paramters
    real(amrex_real), intent(in  ) :: alpha
    integer, intent(inout ) :: sucessful_iter,umbrella_size,sucessful_iter_prev

            ! open(7,file="Console_output_Fortran.txt",status="old",position="append",action="write")

    if (umbrella <=umbrella_min) then
      umbrella_size=1       ! This value serves as a flag when the spring constant cannot be made smaller
    else if (umbrella .GE. umbrella_max .AND. sucessful_iter==0 ) then
      write(6,*)"Umbrella is large, and we have no overlap. Umbrella is ",umbrella
              ! write(7,*)"Umbrella is large, and we have no overlap. Umbrella is ",umbrella
      write(6,*)"Shifting down phi_0 "
              ! write(7,*)"Shifting down phi_0 "
      umbrella_size=2        !This value serves as a flag when the spring constant cannot be made stronger
    end if
    !! check previous step, current step, and umbrella_size. This section must account for all positbiltities
    write(6,*)"umbrella is",umbrella

      if (sucessful_iter_prev==1 .AND. sucessful_iter==1 .AND. umbrella_size==0) then
          umbrella=umbrella/alpha ! Scale down spring constant if both current and previous steps are succesful and the umbrella can be made smaller
      else if(sucessful_iter_prev==1 .AND. sucessful_iter==1 .AND. umbrella_size==1) then
          umbrella=umbrella ! Do not change umbrella when both previous and current steps are successful and the umbrella should not be made smaller
      else if(sucessful_iter_prev==0 .AND. sucessful_iter==1) then
          umbrella=umbrella ! Do not change spring constant if the current step was sucessful and the previous step was not successful
      else if(sucessful_iter_prev==0 .AND. sucessful_iter==0) then
          umbrella=umbrella*alpha ! Scale up the spring constant when neither the current or previous steps were successful
      else if(sucessful_iter_prev==1 .AND. sucessful_iter==0) then
          umbrella=umbrella*alpha ! If current step was not successful, but the previous step was, then the umbrella is set to the previous value.
          umbrella_size=1 ! The spring constant should no longer be changed in this case
      end if

      if (umbrella .LE. umbrella_max .AND. sucessful_iter .NE. 0 ) then
        write(6,*)"Umbrella changed to ",umbrella
                !write(7,*)"Umbrella changed to ",umbrella
      end if
              !close(7)
  end subroutine Umbrella_Adjust





  subroutine Param_Output ( umbrella_output,phi0_output) bind(C,name="Param_Output")
    ! This subroutine outputs the current spring constant value and phi_0 value.
    ! This is usually used when writing quantities to file

    ! INPUT/OUTPUT:
    ! umbrella_output -- a double that is given the value of the current spring constant as an output
    ! phi0_output -- a  double that is given the value of the current phi_0 as an output
    real(amrex_real), intent(out  ) :: umbrella_output,phi0_output

    umbrella_output=umbrella
    phi0_output=phi0

  end subroutine Param_Output



  subroutine umbrella_reset ( umbrella_input) bind(C,name="umbrella_reset")
    ! This subroutine sets the umbrella spring constant when is it called with the value provided in the input
    ! If the adaptive feature is used, the new value should be between the minimum and maximum admissible values.

    !INPUT
    ! umbrella_input -- A new value of the spring constant. The spring constant is set as this value.
    real(amrex_real), intent(in  ) :: umbrella_input

            ! open(7,file="Console_output_Fortran.txt",status="old",position="append",action="write")

    umbrella=umbrella_input

    write(6,*)"Umbrella is set to",umbrella
            ! write(7,*)"Umbrella is set to",umbrella
            ! close(7)
  end subroutine umbrella_reset

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! subroutines changed to 1D
  subroutine Stat_Quant_1D(lo,hi, phi,phi_avg) bind(C,name="Stat_Quant_1D")
    !Computes the spatial average of the phi field provided as an input

    !INPUT:
    ! lo,hi -- multifab ends
    ! phi -- the current phi field multifab

    !OUTPUT:
    ! phi_avg  spatial average of field
    integer         , intent(in   ) :: lo(1),hi(1)
    real(amrex_real), intent(in   ) :: phi(lo(1)-ngc(1):hi(1)+ngc(1))
    real(amrex_real), intent(inout  ) :: phi_avg
    integer :: i

         do  i=lo(1),hi(1)
          phi_avg = phi_avg + phi(i)
        enddo

  end subroutine Stat_Quant_1D



  subroutine integrate_1D(lo,hi, phi, dx, integral) bind(C,name="integrate_1D")
    ! Subroutine that computes the integral term of the umbrella contraint term provided the current phi field.

    !INPUT:
    ! lo,hi -- multifab ends
    ! phi -- the current phi field multifab
    ! dx -- spatial grid spacing array

    !OUTPUT:
    ! integral -- value of the integral umbrella term for the current time-step.

    ! NOTE: phi_0 is specified in the GL namespace list
    integer         , intent(in   ) :: lo(1),hi(1)
    real(amrex_real), intent(in   ) :: dx(1)
    real(amrex_real), intent(out  ) :: integral

    real(amrex_real), intent(inout) :: phi(lo(1)-ngc(1):hi(1)+ngc(1))

    real(amrex_real) :: xloc
    integer :: i

         do  i=lo(1),hi(1)

            integral = integral +  (phi(i) - phi0)*dx(1)  !use 1D mid-point rule to compute integral

        enddo

end subroutine integrate_1D


subroutine initphi_1D(lo,hi, phi, dx) bind(C,name="initphi_1D")
    ! subroutine that sets the initial conition
    ! Set to be a  zero

    !INPUT:
    ! lo,hi -- multifab ends
    ! phi -- the current phi field multifab
    ! dx -- spatial grid spacing array
    ! dt -- time step

    !OUTPUT:
    ! phi -- phi is set as the specified initial condition
  integer         , intent(in   ) :: lo(1),hi(1)
  real(amrex_real), intent(in   ) :: dx(1)

  real(amrex_real), intent(inout) :: phi(lo(1)-ngc(1):hi(1)+ngc(1))

  integer :: i


      do  i=lo(1),hi(1)
             phi(i) = 0.d0
      enddo

          ! open(7,file="Console_output_Fortran.txt") !create output file for fortran output from run
          ! close(7)

end subroutine initphi_1D



subroutine rk2_stage1_1D(lo,hi, phi, phin, rannums, integral, energy, teng, H1_semi, dx, dt,phi_avg) bind(C,name="rk2_stage1_1D")
   ! This subroutine uses an explicit forward euler step with central finite differences to obtain the solution for the
    ! next time step.

    !INPUT:
    ! lo,hi -- multifab ends
    ! phi -- the current phi field multifab
    ! rannums -- the random number multifab
    ! integral -- value of the integral umbrella term for the current time-step.
    ! dx -- spatial grid spacing array
    ! dt -- time step

    !OUTPUT:
    ! energy --  G-L free energy functional with NO umbrella contribution
    ! teng --  G-L free energy functional WITH  umbrella contribution
    ! phi_avg  spatial average of field
    ! phin -- the next time step phi field multifab
    ! phi -- phi is set as the current phi field multifab


    ! NOTE: the integral term is computed with "integral" subroutine and is provided as an input for this functions.
    ! NOTE: that "phi" is set from the previous field value to the current time-step field
    ! NOTE: We assume \mu=1
    ! NOTE : This subroutine calls the Stat_Quant subroutine to compute spatial averages of \phi
  integer         , intent(in   ) :: lo(1),hi(1)
  real(amrex_real), intent(in   ) :: dx(1), dt
  real(amrex_real), intent(inout) :: energy, teng, H1_semi
  real(amrex_real), intent(inout) :: phi_avg


  real(amrex_real), intent(inout) :: phi(lo(1)-ngc(1):hi(1)+ngc(1))
  real(amrex_real), intent(inout) :: phin(lo(1)-ngc(1):hi(1)+ngc(1))

  real(amrex_real), intent(in)    :: rannums(lo(1)-ngc(1):hi(1)+ngc(1))
  real(amrex_real), intent(in)    :: integral

  integer :: i,k,l
  real(amrex_real) :: func,factor, dele

   factor = sqrt(2.d0*noise_coef/(dt*dx(1)))  !noise_coef=k_B T

       do  i=lo(1),hi(1)

         func =  acoef+2.d0*bcoef*phi(i)+3.d0*ccoef*phi(i)**2+4.d0*dcoef*phi(i)**3  !polynomial part of energy functional

         phin(i) = phi(i) + dt*diff_coef*(phi(i+1)-2.d0*phi(i)+phi(i-1))/dx(1)**2  &
          -dt*func-dt*umbrella*integral + dt*factor * rannums(i)   ! Forward euler time-step with central spatial discretization

          dele =  phi(i)*(acoef+bcoef*phi(i)+ccoef*phi(i)**2+dcoef*phi(i)**3) &
            - 0.5d0*diff_coef*phi(i)*(phi(i+1)-2.d0*phi(i)+phi(i-1))/dx(1)**2

          energy = energy + dele*dx(1) !G-L free energy functional with NO umbrella contribution
          teng = teng + ( dele + 0.5d0*umbrella*integral**2 )*dx(1) !G-L free energy functional WITH  umbrella contribution
          H1_semi = H1_semi -0.5d0*diff_coef*phi(i)*(phi(i+1)-2.d0*phi(i)+phi(i-1))/dx(1)

      enddo

       do  i=lo(1),hi(1)

           phi(i)= phin(i)  !rename old field as the new field

      enddo
    call Stat_Quant_1D(lo,hi, phi,phi_avg)
end subroutine rk2_stage1_1D


subroutine setdt_1D (  dx, dt ) bind(C,name="setdt_1D")
  ! Subroutine for determining the time-step for the explicit method.

  ! NOTE: This time-step estimate does not depend on the umbrella parameters. Certain choices of those
  ! parameters may violate stabilitity requirements of the expicit method. If this is the case, either the time-step should
  ! be made smaller, or the umbrella parameters should be set to smaller values.

    !INPUT:
    ! dx -- spatial grid spacing array

    !OUTPUT:
    ! dt -- time step
  double precision::  dt, dx(1)

  if(cfl.gt.0.d0)then

      dt = .25d0*cfl*(dx(1)**2)/diff_coef

  endif

  end subroutine setdt_1D

end module time_step_module



