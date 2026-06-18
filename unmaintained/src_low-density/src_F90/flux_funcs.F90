module flux_module

  use amrex_fort_module, only : amrex_real
  use common_namelist_module, only : ngc, nvars, nprimvars, nspecies, molmass, cell_depth, k_b, Runiv, bc_vel_lo, bc_vel_hi, n_cells, membrane_cell, visc_type, algorithm_type
  use rng_functions_module

  implicit none

  private

  public :: stoch_flux, stoch_flux_BC, rejection_sampler, eval_density, gaussian_density

contains


  subroutine stoch_flux(lo,hi, cons, fluxx, fluxy, &
#if (AMREX_SPACEDIM == 3)
       fluxz, &
#endif
       ranfluxx, ranfluxy, &
#if (AMREX_SPACEDIM == 3)
       ranfluxz, &
#endif
       dx, dt) bind(C,name="stoch_flux")

    integer         , intent(in   ) :: lo(3),hi(3)
    real(amrex_real), intent(in   ) :: dx(3), dt
    real(amrex_real), intent(inout) :: fluxx(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3), nvars)
    real(amrex_real), intent(inout) :: fluxy(lo(1):hi(1),lo(2):hi(2)+1,lo(3):hi(3), nvars)
#if (AMREX_SPACEDIM == 3)
    real(amrex_real), intent(inout) :: fluxz(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)+1, nvars)
#endif

    real(amrex_real), intent(inout) :: ranfluxx(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3), nvars)
    real(amrex_real), intent(inout) :: ranfluxy(lo(1):hi(1),lo(2):hi(2)+1,lo(3):hi(3), nvars)
#if (AMREX_SPACEDIM == 3)
    real(amrex_real), intent(inout) :: ranfluxz(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)+1, nvars)
#endif

    real(amrex_real), intent(in   ) :: cons(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3), nvars)


    real(amrex_real) :: volinv, dtinv, vol, A, B, P, birth, death, x

    real(amrex_real) :: rFracs(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3), nvars)
    integer :: i,j,k,l
    integer :: ll, ns

    birth = 1.0
    death = 1.0

    dtinv = 1d0/dt
#if (AMREX_SPACEDIM == 3)
    vol = dx(1)*dx(2)*dx(3)
    volinv = 1d0/(vol)
#endif

#if (AMREX_SPACEDIM == 2)
    vol = dx(1)*dx(2)*cell_depth
    volinv = 1d0/(vol)
#endif

    if (abs(visc_type) .gt. 1) then

        !get uniform random fraction of the density that goes right
        CALL RANDOM_NUMBER(rFracs)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!! JB's tensor form !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!! x-flux !!!!!!!!!!!!!!!!!!!

! currently storing all fractions going right. optimize later

       do k = lo(3),hi(3)
          do j = lo(2),hi(2)
            !fix ghost cell values - hack for now
            rFracs(lo(1)-1,j,k,1) = rFracs(hi(1),j,k,1)
            rFracs(hi(1)+1,j,k,1) = rFracs(lo(1),j,k,1)
            do i = lo(1),hi(1)+1

                !effective number of particles on left and right
                A = cons(i-1,j,k,1)*vol!*rFracs(i-1,j,k,1)
                B = cons(i,j,k,1)*vol!*(1-rFracs(i,j,k,1))
                !A = cons(i-1,j,k,1)*rFracs(i-1,j,k,1)
                !B = cons(i,j,k,1)*(1-rFracs(i,j,k,1))
                P = A+B

                !draw from distribution
                !print *, "Particles and volume:",A,B, P, vol
                CALL rejection_sampler(birth,death,B,P,x)
                !print *, "dx = ", dx(1)
                ranfluxx(i,j,k,1) = x*volinv
                fluxx(i,j,k,1) = x*volinv

            end do
            ranfluxx(lo(1),j,k,1) = 0
            fluxx(lo(1),j,k,1) = 0
            ranfluxx(hi(1)+1,j,k,1) = 0
            fluxx(hi(1)+1,j,k,1) = 0

          end do
       end do


       ! print*, "Hack: got here (end) stochflux"

!!!!!!!!!!!!!!!!!!! y-flux !!!!!!!!!!!!!!!!!!!

       do k = lo(3),hi(3)
          do j = lo(2),hi(2)+1
             do i = lo(1),hi(1)


             end do
          end do
       end do

!!!!!!!!!!!!!!!!!!! z-flux !!!!!!!!!!!!!!!!!!!

       do k = lo(3),hi(3)+1
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)



             end do
          end do
       end do

    else

    endif

!!!!!!!!!!!!! Enforce flux boundary conditions !!!!!!!!!!!!!

    call stoch_flux_BC(lo,hi, fluxx, fluxy, &
#if (AMREX_SPACEDIM == 3)
         fluxz &
#endif
         )

  end subroutine stoch_flux

  subroutine stoch_flux_BC(lo,hi, xflux, yflux, &
#if (AMREX_SPACEDIM == 3)
       zflux &
#endif
       ) bind(C,name="stoch_flux_bounds")

    integer         , intent(in   ) :: lo(3),hi(3)
    real(amrex_real), intent(inout) :: xflux(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3), nvars)
    real(amrex_real), intent(inout) :: yflux(lo(1):hi(1),lo(2):hi(2)+1,lo(3):hi(3), nvars)
#if (AMREX_SPACEDIM == 3)
    real(amrex_real), intent(inout) :: zflux(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)+1, nvars)
#endif

    real(amrex_real) :: sqrtTwo

    integer :: i,j,k,l

    sqrtTwo = sqrt(2.0)

!!!!!!!!!!!!!! x-flux BCs !!!!!!!!!!!!!!

    !if on lower bound and specular
    if((lo(1) .eq. 0) .and. (bc_vel_lo(1) .eq. 1)) then
       do l = 2,nvars
          do k = lo(3),hi(3)
             do j = lo(2),hi(2)

                xflux(0,j,k,l) = 0

             end do
          end do
       end do
    endif

    !if on upper bound and specular
    if((hi(1) .eq. n_cells(1)-1) .and. (bc_vel_hi(1) .eq. 1)) then
       do l = 2,nvars
          do k = lo(3),hi(3)
             do j = lo(2),hi(2)

                xflux(hi(1)+1,j,k,l) = 0

             end do
          end do
       end do
    endif

    !if on lower bound and diff
    if((lo(1) .eq. 0) .and. (bc_vel_lo(1) .eq. 2)) then
       do l = 2,nvars
          do k = lo(3),hi(3)
             do j = lo(2),hi(2)

                xflux(0,j,k,l) = sqrtTwo*xflux(0,j,k,l)

             end do
          end do
       end do
    endif

    !if on upper bound and diff
    if((hi(1) .eq. n_cells(1)-1) .and. (bc_vel_hi(1) .eq. 2)) then
       do l = 2,nvars
          do k = lo(3),hi(3)
             do j = lo(2),hi(2)

                xflux(hi(1)+1,j,k,l) = sqrtTwo*xflux(hi(1)+1,j,k,l)

             end do
          end do
       end do
    endif

!!!!!!!!!!!!!! y-flux BCs !!!!!!!!!!!!!!

    !if on lower bound and specular
    if((lo(2) .eq. 0) .and. (bc_vel_lo(2) .eq. 1)) then
       do l = 2,nvars
          do k = lo(3),hi(3)
             do i = lo(1),hi(1)

                yflux(i,0,k,l) = 0

             end do
          end do
       end do
    endif

    !if on upper bound and specular
    if((hi(2) .eq. n_cells(2)-1) .and. (bc_vel_hi(2) .eq. 1)) then
       do l = 2,nvars
          do k = lo(3),hi(3)
             do i = lo(1),hi(1)

                yflux(i,hi(2)+1,k,l) = 0

             end do
          end do
       end do
    endif

    !if on lower bound and diff
    if((lo(2) .eq. 0) .and. (bc_vel_lo(2) .eq. 2)) then
       do l = 2,nvars
          do k = lo(3),hi(3)
             do i = lo(1),hi(1)

                yflux(i,0,k,l) = sqrtTwo*yflux(i,0,k,l)

             end do
          end do
       end do
    endif

    !if on upper bound and diff
    if((hi(2) .eq. n_cells(2)-1) .and. (bc_vel_hi(2) .eq. 2)) then
       do l = 2,nvars
          do k = lo(3),hi(3)
             do i = lo(1),hi(1)

                yflux(i,hi(2)+1,k,l) = sqrtTwo*yflux(i,hi(2)+1,k,l)

             end do
          end do
       end do
    endif

!!!!!!!!!!!!!! z-flux BCs !!!!!!!!!!!!!!

    !if on lower bound and specular
    if((lo(3) .eq. 0) .and. (bc_vel_lo(3) .eq. 1)) then
       do l = 2,nvars
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)

                zflux(i,j,0,l) = 0

             end do
          end do
       end do
    endif

    !if on upper bound and specular
    if((hi(3) .eq. n_cells(3)-1) .and. (bc_vel_hi(3) .eq. 1)) then
       do l = 2,nvars
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)

                zflux(i,j,hi(3)+1,l) = 0

             end do
          end do
       end do
    endif

    !if on lower bound and diff
    if((lo(3) .eq. 0) .and. (bc_vel_lo(3) .eq. 2)) then
       do l = 2,nvars
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)

                zflux(i,j,0,l) = sqrtTwo*zflux(i,j,0,l)

             end do
          end do
       end do
    endif

    !if on upper bound and diff
    if((hi(3) .eq. n_cells(3)-1) .and. (bc_vel_hi(3) .eq. 2)) then
       do l = 2,nvars
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)

                zflux(i,j,hi(3)+1,l) = sqrtTwo*zflux(i,j,hi(3)+1,l)

             end do
          end do
       end do
    endif


  end subroutine stoch_flux_BC

  subroutine rejection_sampler(b,d,N,P,x)

    real(amrex_real), intent(in   ) :: b,d,N,P
    real(amrex_real), intent(inout) :: x

    real(amrex_real) :: mean, variance, M, U, proposal, F,G
    integer :: Gtrials, Utrials, trialCount

    Gtrials = 500
    Utrials = 500
    trialCount = 0
    M = 1  !might be big/small, experiment with this

    !get the mean and variance of a proposal Gaussian
    !print *, b,d,N,P
    mean = (b*P)/(b+d)-N
    !print *, mean
    variance = b*d*P/(b+d)**2
    proposal = 0
    F = 0
    G = 0

    !sample from the proposal Gaussian distribution - accept/rejection
    do while (trialCount < Gtrials)
        trialCount = trialCount + 1

        !!if using uniform proposal
        CALL RANDOM_NUMBER(proposal)
        proposal = P*proposal - N
        G = 1

        !!if using gaussian proposal
        !proposal = get_fhd_normal_func() !make sure this is actually N(0,1)
        !proposal = sqrt(variance)*proposal+mean
        !CALL gaussian_density(mean, variance, proposal, G)

        !print *, "Proposal with bounds", proposal, -N, P-N

        CALL RANDOM_NUMBER(U)
        CALL eval_density(b,d,N,P,proposal,F)

        !print *, 'ratio', F/G, F, G

        if (U < F/(M*G)) then
            x = proposal
            EXIT
        endif

    end do



  end subroutine rejection_sampler

  subroutine eval_density(b,d,N,P,x,f)

    real(amrex_real), intent(in   ) :: b,d,N,P,x
    real(amrex_real), intent(inout) :: f

    if (x .gt. -N .and. x .lt. P-N) then
        f = (x+N)*log(d/b*(x+N)/(P-(x+N)))
        !print *, f
        f = f +P*log(P-(x+N))
        !print *, f
        f = f +P*log((b+d)/(d*P))
        !print *,f
        f = exp(-f)
        !print *,f
    else
        f = 0
        !print *, 'outside'
    endif

  end subroutine eval_density

  subroutine gaussian_density(mean,var,x,f)

    real(amrex_real), intent(in   ) :: mean, var, x
    real(amrex_real), intent(inout) :: f

    real(amrex_real) :: pi

    pi=4.D0*DATAN(1.D0)
    !print *,pi, x-mean, var

    !f = 1.0/sqrt(2*pi*var)*exp(-1.0/(2*var)*(x-mean)**2)
    f = exp(-1.0/(2*var)*(x-mean)**2)

  end subroutine gaussian_density

end module flux_module


































