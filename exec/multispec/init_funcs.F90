  ! IMPORTANT
  ! In the low-Mach code (init_rho_and_umac), c_init specifies initial MASS FRACTIONS
  ! (should sum to unity!... but we overwrite the final concentration so sum(c_i)=1 before computing rho)
  ! The density follows from the EOS in the LM case so it cannot be specified
  ! Same applies to boundary conditions

  ! prob_types codes for init_lowmach:

  !=============================================================
  ! case 1:

  !=========================================================
  ! case 2:


#if(AMREX_SPACEDIM == 2)

subroutine init_rho_and_umac(lo,hi, &
                             rho,rlo,rhi,nspecies, &
                             u,ulo,uhi, v,vlo,vhi, &
                             dx, prob_lo, prob_hi) bind(C, name="init_rho_and_umac")

  use amrex_error_module
  use amrex_fort_module, only : amrex_real
  use common_namelist_module, only: prob_type, algorithm_type, T_init, &
                                    rho0, rhobar, initial_variance_mass, &
                                    smoothing_width
  use multispec_namelist_module, only: Dbar, Dtherm, charge_per_mass, &
                                       c_init_1, c_init_2, temp_type, is_ideal_mixture

  implicit none

  integer         , intent(in   ) :: lo(2), hi(2), nspecies
  integer         , intent(in   ) :: rlo(2),rhi(2)
  integer         , intent(in   ) :: ulo(2),uhi(2), vlo(2),vhi(2)
  double precision, intent(in   ) :: prob_lo(2)
  double precision, intent(in   ) :: prob_hi(2)
  double precision, intent(in   ) :: dx(2)
  double precision, intent(inout) :: rho(rlo(1):rhi(1),rlo(2):rhi(2),nspecies)
  double precision, intent(inout) ::   u(ulo(1):uhi(1),ulo(2):uhi(2))
  double precision, intent(inout) ::   v(vlo(1):vhi(1),vlo(2):vhi(2))

  ! local concentration
  double precision :: c(rlo(1):rhi(1),rlo(2):rhi(2),nspecies)

  ! local variables
  integer          :: i,j,n
  double precision :: half,x,y,rad,L(2),sumtot,c_loc,r,r1,r2,m_e,l1,l2

  double precision :: random

  double precision :: rho_total, sinx, sinz

  double precision :: x1, x2, y1, y2, coeff

  half = 0.5d0

  L(1:2) = prob_hi(1:2)-prob_lo(1:2) ! Domain length

  select case (abs(prob_type))

  case (1) 

     !=============================================================
     ! bubble with radius = 1/4 of domain in x
     ! c=c_init_1(:) inside, c=c_init_2(:) outside
     ! can be discontinous or smooth depending on smoothing_width
     !=============================================================

     u = 0.d0
     v = 0.d0

     rad = L(1)/4.d0

     !$omp parallel do private(i,j,k,x,y,z,r)
     do j=lo(2),hi(2)
        y = prob_lo(2) + (dble(j)+half)*dx(2) - half*(prob_lo(2)+prob_hi(2))
        do i=lo(1),hi(1)
           x = prob_lo(1) + (dble(i)+half)*dx(1) - half*(prob_lo(1)+prob_hi(1))

           r = sqrt(x**2 + y**2)

           if (smoothing_width .eq. 0) then

              ! discontinuous interface
              if (r .lt. rad) then
                 c(i,j,1:nspecies) = c_init_1(1:nspecies)
              else
                 c(i,j,1:nspecies) = c_init_2(1:nspecies)
              end if

           else

              ! smooth interface
              c(i,j,1:nspecies-1) = c_init_1(1:nspecies-1) + &
                   (c_init_2(1:nspecies-1) - c_init_1(1:nspecies-1))* &
                   0.5d0*(1.d0 + tanh((r-rad)/(smoothing_width*dx(1))))

           end if

        end do
     end do
     !$omp end parallel do

  case (4) 

     !=============================================================
     ! bubble with radius = 1/4 of domain in x
     ! c=c_init_1(:) inside, c=c_init_2(:) outside
     ! can be discontinous or smooth depending on smoothing_width
     !=============================================================

     u = 0.d0
     v = 0.d0

     l1 = L(2)/3
     l2 = 2*l1

     !$omp parallel do private(i,j,k,x,y,z,r)
     do j=lo(2),hi(2)
        y = prob_lo(2) + (dble(j)+half)*dx(2)

        do i=lo(1),hi(1)
           x = prob_lo(1) + (dble(i)+half)*dx(1)

           if (smoothing_width .eq. 0) then

              ! discontinuous interface
              if (y .lt. l1) then
                 c(i,j,1:nspecies) = c_init_1(1:nspecies)
              elseif(y .lt. l2) then
                 c(i,j,1:nspecies) = c_init_2(1:nspecies)
              else
                 c(i,j,1:nspecies) = c_init_1(1:nspecies)
              endif

           else

              ! smooth interface
              c(i,j,1:nspecies-1) = c_init_1(1:nspecies-1) + &
                   (c_init_2(1:nspecies-1) - c_init_1(1:nspecies-1))* &
                   (1/(1+Exp(-smoothing_width*(y-l1))) - 1/(1+Exp(-smoothing_width*(y-l2))))

           end if           

        end do
     end do

  case (2) 

     !=========================================================
     ! constant concentration gradient along y
     ! c=c_init_1(:) on bottom, c=c_init_2(:) on top
     !=========================================================

     u = 0.d0
     v = 0.d0

     !$omp parallel do private(i,j,k,x,y,z)
     do j=lo(2),hi(2)
        y = prob_lo(2) + (dble(j)+half)*dx(2)
        do i=lo(1),hi(1)
           x = prob_lo(1) + (dble(i)+half)*dx(1)

           c(i,j,1:nspecies) = c_init_1(1:nspecies) + &
                (c_init_2(1:nspecies) - c_init_1(1:nspecies))*(y-prob_lo(2))/L(2)

        end do
     end do
     !$omp end parallel do

  case (12)

       !=============================================================
       ! Gaussian bubble centered in domain
       ! c=c_init_1(:) inside; c=c_init_2(:) outside
       ! lo- and hi-y walls move with prescribed velocity,
       ! see inhomogeneous_bc_val.f90
       !=============================================================

       u = 0.d0
       v = 0.d0

       do j=lo(2),hi(2)
          y = prob_lo(2) + dx(2) * (dble(j)+0.5d0) - 0.5d0*(prob_lo(2)+prob_hi(2))
          do i=lo(1),hi(1)
             x = prob_lo(1) + dx(1) * (dble(i)+0.5d0) - 0.5d0*(prob_lo(1)+prob_hi(1))

             r = sqrt (x**2 + y**2)

             ! set c using Gaussian bump
             c(i,j,1:nspecies-1) = c_init_1(1:nspecies-1)*exp(-75.d0*r**2)

          enddo
       enddo

  case (15)

  !=========================================================
  ! case +/-15: mostly for testing electrodiffusion
  ! Discontinuous band in central 1/2 (case 15)
  ! c=c_init_1(:) inside; c=c_init_2(:) outside
  !=========================================================

     u = 0.d0
     v = 0.d0
     
     ! first quarter of domain
     y1 = (3*prob_lo(2) + prob_hi(2)) / 4.d0
     x1 = (3*prob_lo(1) + prob_hi(1)) / 4.d0
        
     ! last quarter of domain
     y2 = (prob_lo(2) + 3*prob_hi(2)) / 4.d0
     x2 = (prob_lo(1) + 3*prob_hi(1)) / 4.d0

     if (smoothing_width .gt. 0.d0) then

          ! smoothed version
          do j=lo(2),hi(2)
             y = prob_lo(2) + dx(2)*(dble(j)+0.5d0) - y1
             do i=lo(1),hi(1)
                x = prob_lo(1) + (dble(i)+0.5d0)*dx(1) - x1
                do n=1,nspecies

                   ! tanh smoothing in y
                   coeff=0.5d0*(tanh(y/(smoothing_width*dx(2)))+1.d0)*0.5d0*(tanh((-y+y2-y1)/(smoothing_width*dx(2)))+1.d0)

                   ! Donev: prob_type = -15: a special case for doing ternary diffusion NaCl + KCl
                   ! Here the last two species have a tanh profile in both x and y (species are Na+,Cl-,K+,water)
                   ! Aadd a tanh smoothing for central 50% of domain in x for second-to-last species
                   if( (prob_type==-15) .and. (n==nspecies-1) ) then 
                      coeff=0.5d0*(tanh(x/(smoothing_width*dx(1)))+1.d0)*0.5d0*(tanh((-x+x2-x1)/(smoothing_width*dx(1)))+1.d0)*coeff
                   end if  
 
                   ! smooth between c_init_1(:) and c_init_2(:)
                   c_loc = c_init_2(n) + (c_init_1(n)-c_init_2(n))*coeff
                   c(i,j,n) = c_loc

                   ! for 4-species test, need to add Cl to central square to balance the K
                   if ( (prob_type==-15) .and. (nspecies .eq. 4) .and. (n .eq. nspecies-1) ) then
                      c(i,j,2) = c(i,j,2) - charge_per_mass(3)/charge_per_mass(2)*c_loc
                   end if

                end do
             end do
          end do

       end if
     
  case default

     call amrex_error("Desired prob_type not supported in 3D")

  end select

  do j=lo(2),hi(2)
     do i=lo(1),hi(1)

        ! set final c_i such that sumtot(c_i) = 1 to within roundoff
        sumtot = 0.d0
        do n=1,nspecies-1
           sumtot = sumtot + c(i,j,n)
        end do
        c(i,j,nspecies) = 1.d0 - sumtot

        ! calculate rho_total from eos
        if (algorithm_type .eq. 6) then
           rho_total = rho0
        else
           sumtot = 0.d0
           do n=1,nspecies
              ! sumtot represents rhoinv
              sumtot = sumtot + c(i,j,n)/rhobar(n)
           end do
           rho_total = 1.d0/sumtot
        end if

        ! calculate rho_i
        rho(i,j,1:nspecies) = rho_total*c(i,j,1:nspecies)

     end do
  end do

end subroutine init_rho_and_umac

#endif


#if(AMREX_SPACEDIM == 3)

subroutine init_rho_and_umac(lo,hi, &
                             rho,rlo,rhi,nspecies, &
                             u,ulo,uhi, v,vlo,vhi, w,wlo,whi, &
                             dx, prob_lo, prob_hi) bind(C, name="init_rho_and_umac")

  use amrex_error_module
  use amrex_fort_module, only : amrex_real
  use common_namelist_module, only: prob_type, algorithm_type, T_init, &
                                    rho0, rhobar, initial_variance_mass, &
                                    smoothing_width
  use multispec_namelist_module, only: Dbar, Dtherm, &
                                       c_init_1, c_init_2, temp_type, is_ideal_mixture

  implicit none

  integer         , intent(in   ) :: lo(3), hi(3), nspecies
  integer         , intent(in   ) :: rlo(3),rhi(3)
  integer         , intent(in   ) :: ulo(3),uhi(3), vlo(3),vhi(3), wlo(3),whi(3)
  double precision, intent(in   ) :: prob_lo(3)
  double precision, intent(in   ) :: prob_hi(3)
  double precision, intent(in   ) :: dx(3)
  double precision, intent(inout) :: rho(rlo(1):rhi(1),rlo(2):rhi(2),rlo(3):rhi(3),nspecies)
  double precision, intent(inout) ::   u(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
  double precision, intent(inout) ::   v(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
  double precision, intent(inout) ::   w(wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))

  ! local concentration
  double precision :: c(rlo(1):rhi(1),rlo(2):rhi(2),rlo(3):rhi(3),nspecies)

  ! local variables
  integer          :: i,j,k,n
  double precision :: half,x,y,z,rad,L(3),sumtot,c_loc,y1,z1,z2,r,r1,r2,m_e

  double precision :: random, coeff

  double precision :: rho_total, sinx, sinz

  half = 0.5d0

  L(1:3) = prob_hi(1:3)-prob_lo(1:3) ! Domain length

  select case (abs(prob_type))

  case (1) 

     !=============================================================
     ! bubble with radius = 1/4 of domain in x
     ! c=c_init_1(:) inside, c=c_init_2(:) outside
     ! can be discontinous or smooth depending on smoothing_width
     !=============================================================

     u = 0.d0
     v = 0.d0
     w = 0.d0

     rad = L(1)/4.d0

     !$omp parallel do private(i,j,k,x,y,z,r)
     do k=lo(3),hi(3)
        z = prob_lo(3) + (dble(k)+half)*dx(3) - half*(prob_lo(3)+prob_hi(3))
        do j=lo(2),hi(2)
           y = prob_lo(2) + (dble(j)+half)*dx(2) - half*(prob_lo(2)+prob_hi(2))
           do i=lo(1),hi(1)
              x = prob_lo(1) + (dble(i)+half)*dx(1) - half*(prob_lo(1)+prob_hi(1))

              r = sqrt(x**2 + y**2 + z**2)

              if (smoothing_width .eq. 0) then

                 ! discontinuous interface
                 if (r .lt. rad) then
                    c(i,j,k,1:nspecies) = c_init_1(1:nspecies)
                 else
                    c(i,j,k,1:nspecies) = c_init_2(1:nspecies)
                 end if

              else

                 ! smooth interface
                 c(i,j,k,1:nspecies-1) = c_init_1(1:nspecies-1) + &
                      (c_init_2(1:nspecies-1) - c_init_1(1:nspecies-1))* &
                      0.5d0*(1.d0 + tanh((r-rad)/(smoothing_width*dx(1))))

              end if

           end do
        end do
     end do
     !$omp end parallel do

  case (2) 

     !=========================================================
     ! constant concentration gradient along y
     ! c=c_init_1(:) on bottom, c=c_init_2(:) on top
     !=========================================================

     u = 0.d0
     v = 0.d0
     w = 0.d0

     !$omp parallel do private(i,j,k,x,y,z)
     do k=lo(3),hi(3)
        z = prob_lo(3) + (dble(k)+half)*dx(3) 
        do j=lo(2),hi(2)
           y = prob_lo(2) + (dble(j)+half)*dx(2)
           do i=lo(1),hi(1)
              x = prob_lo(1) + (dble(i)+half)*dx(1)

              c(i,j,k,1:nspecies) = c_init_1(1:nspecies) + &
                   (c_init_2(1:nspecies) - c_init_1(1:nspecies))*(y-prob_lo(2))/L(2)

           end do
        end do
     end do
     !$omp end parallel do

  case (12)
     
     !=============================================================
     ! Gaussian bubble centered in domain
     ! c=c_init_1(:) inside; c=c_init_2(:) outside
     ! lo- and hi-y walls move with prescribed velocity,
     ! see inhomogeneous_bc_val.f90
     !=============================================================
     
     u = 0.d0
     v = 0.d0
     w = 0.d0

     do k=lo(3),hi(3)
        z = prob_lo(3) + dx(3) * (dble(k)+0.5d0) - 0.5d0*(prob_lo(3)+prob_hi(3))
        do j=lo(2),hi(2)
           y = prob_lo(2) + dx(2) * (dble(j)+0.5d0) - 0.5d0*(prob_lo(2)+prob_hi(2))
           do i=lo(1),hi(1)
              x = prob_lo(1) + dx(1) * (dble(i)+0.5d0) - 0.5d0*(prob_lo(1)+prob_hi(1))

              r = sqrt (x**2 + y**2 + z**2)

              ! set c using Gaussian bump
              c(i,j,k,1:nspecies-1) = c_init_1(1:nspecies-1)*exp(-75.d0*r**2)

           enddo
        enddo
     enddo

  case (15)

     !=========================================================
     ! case +/-15: mostly for testing electrodiffusion
     ! Discontinuous band in central 1/2 (case 15)
     ! c=c_init_1(:) inside; c=c_init_2(:) outside
     !=========================================================

     u = 0.d0
     v = 0.d0
     w = 0.d0
     
     ! first quarter of domain
     z1 = (3*prob_lo(3) + prob_hi(3)) / 4.d0
        
     ! last quarter of domain
     z2 = (prob_lo(3) + 3*prob_hi(3)) / 4.d0

     if (smoothing_width .gt. 0.d0) then

        ! smoothed version
        do k=lo(3),hi(3)
             z = prob_lo(3) + dx(3)*(dble(k)+0.5d0) - z1
             do j=lo(2),hi(2)
                do i=lo(1),hi(1)
                   do n=1,nspecies

                      ! tanh smoothing in z
                      coeff=0.5d0*(tanh(z/(smoothing_width*dx(2)))+1.d0)*0.5d0*(tanh((-z+z2-z1)/(smoothing_width*dx(3)))+1.d0)
 
                      ! smooth between c_init_1(:) and c_init_2(:)
                      c_loc = c_init_2(n) + (c_init_1(n)-c_init_2(n))*coeff
                      c(i,j,k,n) = c_loc

                   end do
                end do
             end do
          end do

       end if
     
  case default

     call amrex_error("Desired prob_type not supported in 3D")

  end select

  do k=lo(3),hi(3)
     do j=lo(2),hi(2)
        do i=lo(1),hi(1)

           ! set final c_i such that sumtot(c_i) = 1 to within roundoff
           sumtot = 0.d0
           do n=1,nspecies-1
              sumtot = sumtot + c(i,j,k,n)
           end do
           c(i,j,k,nspecies) = 1.d0 - sumtot

           ! calculate rho_total from eos
           if (algorithm_type .eq. 6) then
              rho_total = rho0
           else
              sumtot = 0.d0
              do n=1,nspecies
                 ! sumtot represents rhoinv
                 sumtot = sumtot + c(i,j,k,n)/rhobar(n)
              end do
              rho_total = 1.d0/sumtot
           end if

           ! calculate rho_i
           rho(i,j,k,1:nspecies) = rho_total*c(i,j,k,1:nspecies)

        end do
     end do
  end do

end subroutine init_rho_and_umac

#endif
