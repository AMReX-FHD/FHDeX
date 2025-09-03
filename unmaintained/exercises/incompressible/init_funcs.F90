subroutine init_vel(lo, hi, vel, vello, velhi, dx, prob_lo, prob_hi, di, &
                    reallo, realhi) bind(C, name="init_vel")

  use amrex_fort_module, only : amrex_real
  use common_namelist_module, only: prob_type

  implicit none

  integer         , intent(in   ) :: lo(3), hi(3), vello(3), velhi(3), di
  real(amrex_real), intent(inout) :: vel(vello(1):velhi(1),vello(2):velhi(2),vello(3):velhi(3))
  real(amrex_real), intent(in   ) :: reallo(3), realhi(3)
  real(amrex_real), intent(in   ) :: prob_lo(3)
  real(amrex_real), intent(in   ) :: prob_hi(3)
  real(amrex_real), intent(in   ) :: dx(3)

  integer          :: i,j,k
  double precision :: pos(3),center(3),partdom,itVec(3),relpos(3),rad,rad2,zshft
  double precision :: L_hlf, k1, k1_inv, k2, k2_inv, r_a, r_b
  double precision :: pi, freq, amp, width1, perturb, slope, fun_ptrb
  double precision :: velx_bwd_temp, velx_fwd_temp

#if (AMREX_SPACEDIM == 2)
  zshft = 0.0d0
#elif (AMREX_SPACEDIM == 3)
  zshft = 0.5d0
#endif

  center = (realhi - reallo)/2d0

  !! IC parameters
  L_hlf = (realhi(1) - reallo(1))/2d0
  ! k1 & k2 determine steepness of velocity profile:
  k1 = 1d-2*L_hlf
  ! k1 = 1d-6*L_hlf
  k2 = k1
  k1_inv = 1/k1
  k2_inv = 1/k2

  ! Vortex:
  ! [r_a r_b] defines radial bounds of velocity bump:
  r_a = 0.35d0*L_hlf
  r_b = L_hlf - r_a

  ! Kelvin-Helmholtz:
  pi = acos(-1.d0)
  freq = 02.d0*pi/L_hlf
  amp = 2.0d-3*L_hlf
  ! amp = 2.0d-1*L_hlf
  width1 = L_hlf/2.0d0

  if (di .EQ. 0) then

     SELECT CASE (prob_type)
     CASE (0)
        vel = 0.d0
     CASE (1)
        !! Vortex:
        do k = lo(3), hi(3)
           do j = lo(2), hi(2)
              do i = lo(1), hi(1) + 1
                 itVec(1) = dble(i)*dx(1)
                 itVec(2) = (dble(j)+0.5d0)*dx(2)
                 itVec(3) = (dble(k)+zshft)*dx(3)

                 pos = reallo + itVec
                 relpos = pos - center
                 rad2 = DOT_PRODUCT(relpos(1:2),relpos(1:2))
                 rad = SQRT(rad2)

                 ! Multiply velocity magnitude by sin(theta)
                 vel(i,j,k) = 0.25d0*(1d0+tanh(k1_inv*(rad-r_a)))*(1d0+tanh(k2_inv*(r_b-rad))) &
                      *(relpos(2)/rad)

              end do
           end do
        end do
     CASE (2)
        !! KH, sine:
        do k = lo(3), hi(3)
           do j = lo(2), hi(2)
              do i = lo(1), hi(1) + 1
                 itVec(1) = dble(i)*dx(1)
                 itVec(2) = (dble(j)+0.5d0)*dx(2)
                 itVec(3) = (dble(k)+zshft)*dx(3)

                 pos = reallo + itVec
                 relpos = pos - center
                 rad2 = DOT_PRODUCT(relpos(1:2),relpos(1:2))
                 rad = SQRT(rad2)

                 perturb = amp*sin(freq*relpos(1))
                 fun_ptrb = 0.25d0*(1d0+tanh(k1_inv*(relpos(2) - (-width1/2.d0+perturb)))) &
                      *(1d0+tanh(k2_inv*((width1/2.d0+perturb) - relpos(2))))
                 vel(i,j,k) = fun_ptrb
              end do
           end do
        end do
     CASE (3)
        ! KH, smooth:
        do k = lo(3), hi(3)
           do j = lo(2), hi(2)
              do i = lo(1), hi(1) + 1
                 itVec(1) = dble(i)*dx(1)
                 itVec(2) = (dble(j)+0.5d0)*dx(2)
                 itVec(3) = (dble(k)+zshft)*dx(3)

                 pos = reallo + itVec
                 relpos = pos - center
                 rad2 = DOT_PRODUCT(relpos(1:2),relpos(1:2))
                 rad = SQRT(rad2)

                 fun_ptrb = 0.25d0*(1d0+tanh(k1_inv*(relpos(2) - (-width1/2.d0)))) &
                      *(1d0+tanh(k2_inv*((width1/2.d0) - relpos(2))))
                 vel(i,j,k) = fun_ptrb

                 ! vel(i,j,k) = 0.d0
              end do
           end do
        end do
     CASE DEFAULT
        print*, "Error: Invalid prob_type"
     END SELECT

  endif

  if (di .EQ. 1) then

     SELECT CASE (prob_type)
     CASE (0)
        vel = 0.d0
     CASE (1)
        !! Vortex:

        do k = lo(3), hi(3)
           do j = lo(2), hi(2) + 1
              do i = lo(1), hi(1)

                 itVec(1) = (dble(i)+0.5d0)*dx(1)
                 itVec(2) = dble(j)*dx(2)
                 itVec(3) = (dble(k)+zshft)*dx(3)

                 pos = reallo + itVec
                 relpos = pos - center
                 rad2 = DOT_PRODUCT(relpos(1:2),relpos(1:2))
                 rad = SQRT(rad2)

                 ! Multiply velocity magnitude by -cos(theta)
                 vel(i,j,k) = 0.25d0*(1d0+tanh(k1_inv*(rad-r_a)))*(1d0+tanh(k2_inv*(r_b-rad))) &
                      *(-relpos(1)/rad)

              end do
           end do
        end do

     CASE (2)
        !! KH, sine:
        do k = lo(3), hi(3)
           do j = lo(2), hi(2) + 1
              do i = lo(1), hi(1)

                 itVec(1) = (dble(i)+0.5d0)*dx(1)
                 itVec(2) = dble(j)*dx(2)
                 itVec(3) = (dble(k)+zshft)*dx(3)

                 pos = reallo + itVec
                 relpos = pos - center
                 rad2 = DOT_PRODUCT(relpos(1:2),relpos(1:2))
                 rad = SQRT(rad2)

                 perturb = amp*sin(freq*relpos(1))
                 slope = amp*freq*cos(freq*relpos(1))
                 fun_ptrb = 0.25d0*(1d0+tanh(k1_inv*(relpos(2) - (-width1/2.d0+perturb)))) &
                      *(1d0+tanh(k2_inv*((width1/2.d0+perturb) - relpos(2))))
                 vel(i,j,k) = slope*fun_ptrb

              end do
           end do
        end do
     CASE (3)
        !! KH, smooth:
        do k = lo(3), hi(3)
           do j = lo(2), hi(2) + 1
              do i = lo(1), hi(1)

                 vel(i,j,k) = 0.d0

              end do
           end do
        end do
     CASE DEFAULT
        print*, "Error: Invalid prob_type"
     END SELECT

  endif

  if (di .EQ. 2) then

     SELECT CASE (prob_type)
     CASE (0)
        vel = 0.d0
     CASE (1)
        do k = lo(3), hi(3) + 1
           do j = lo(2), hi(2)
              do i = lo(1), hi(1)

                 vel(i,j,k) = 0.0d0

              end do
           end do
        end do
     CASE (2)
        do k = lo(3), hi(3) + 1
           do j = lo(2), hi(2)
              do i = lo(1), hi(1)

                 vel(i,j,k) = 0.0d0

              end do
           end do
        end do
     CASE (3)
        do k = lo(3), hi(3) + 1
           do j = lo(2), hi(2)
              do i = lo(1), hi(1)

                 vel(i,j,k) = 0.0d0

              end do
           end do
        end do
     CASE DEFAULT
        print*, "Error: Invalid prob_type"
     END SELECT

  endif


end subroutine init_vel


  ! IMPORTANT: In the diffusion only code (init_rho), c_init specifies initial values for DENSITY
  ! In the low-Mach code (init_rho_and_umac), c_init specifies initial MASS FRACTIONS
  ! (should sum to unity!... but we overwrite the final concentration so sum(c_i)=1 before computing rho)
  ! The density follows from the EOS in the LM case so it cannot be specified
  ! Same applies to boundary conditions

  ! prob_types codes for init_lowmach:

  !=============================================================
  ! case 1:
  ! bubble with radius = 1/4 of domain in x
  ! c=c_init_1(:) inside, c=c_init_2(:) outside
  ! can be discontinous or smooth depending on smoothing_width

  !=========================================================
  ! case 2:
  ! constant concentration gradient along y
  ! c=c_init_1(:) on bottom, c=c_init_2(:) on top


#if(AMREX_SPACEDIM == 2)

subroutine init_rho_and_umac(lo,hi, &
                             rho,rlo,rhi,nspecies, &
                             u,ulo,uhi, v,vlo,vhi, &
                             dx, prob_lo, prob_hi, &
                             reallo, realhi) bind(C, name="init_rho_and_umac")

  use amrex_error_module
  use amrex_fort_module, only : amrex_real
  use common_namelist_module, only: prob_type, algorithm_type, T_init, &
                                    rho0, rhobar, initial_variance_mass, &
                                    smoothing_width, n_cells
  use multispec_namelist_module, only: Dbar, Dtherm, &
                                       c_init, temp_type, is_ideal_mixture

  implicit none

  integer         , intent(in   ) :: lo(2), hi(2), nspecies
  integer         , intent(in   ) :: rlo(2),rhi(2)
  integer         , intent(in   ) :: ulo(2),uhi(2), vlo(2),vhi(2)
  double precision, intent(in   ) :: reallo(2), realhi(2)
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
  double precision :: half,x,y,rad,L(2),sumtot,c_loc,y1,r,r1,r2,m_e,l1,l2

  double precision :: random

  double precision :: rho_total, sinx, sinz

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

     ! print *, "Hack: smoothing width = ", smoothing_width
     ! print *, "Hack: nspecies = ", nspecies
     ! print *, "Hack: dx = ", dx(1), " dy = ", dx(2)
     ! print *, "Hack: c_init 1 = ", c_init_1(1:nspecies)
     ! print *, "Hack: c_init 2 = ", c_init_2(1:nspecies)

     !$omp parallel do private(i,j,k,x,y,z,r)
     do j=lo(2),hi(2)
        y = prob_lo(2) + (dble(j)+half)*dx(2) - half*(prob_lo(2)+prob_hi(2))
        do i=lo(1),hi(1)
           x = prob_lo(1) + (dble(i)+half)*dx(1) - half*(prob_lo(1)+prob_hi(1))

           r = sqrt(x**2 + y**2)

           ! print *, "Hack: x = ", x, " y = ", y

           if (smoothing_width .eq. 0) then

              ! discontinuous interface
              if (r .lt. rad) then
                 c(i,j,1:nspecies) = c_init_1(1:nspecies)
              else
                 c(i,j,1:nspecies) = c_init_2(1:nspecies)
              end if

              ! print *, "Hack: c = ", c(i,j,1:nspecies)
              ! print *, "Hack: r = ", r, "Hack: r_ad = ", rad

           else

              ! smooth interface
              c(i,j,1:nspecies-1) = c_init_1(1:nspecies-1) + &
                   (c_init_2(1:nspecies-1) - c_init_1(1:nspecies-1))* &
                   0.5d0*(1.d0 + tanh((r-rad)/(smoothing_width*dx(1))))

              ! print *, "Hack: c = ", c(i,j,1:nspecies-1)

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

     l1 = 0.35*L(2)
     l2 = l1 + 0.3*L(2)

     ! print *, "Hack: smoothing width = ", smoothing_width
     ! print *, "Hack: nspecies = ", nspecies
     ! print *, "Hack: dx = ", dx(1), " dy = ", dx(2)
     ! print *, "Hack: c_init 1 = ", c_init_1(1:nspecies)
     ! print *, "Hack: c_init 2 = ", c_init_2(1:nspecies)

     !$omp parallel do private(i,j,k,x,y,z,r)
     do j=lo(2),hi(2)
        y = prob_lo(2) + (dble(j)+half)*dx(2)

        !print *, "y: ", y
        do i=lo(1),hi(1)
           x = prob_lo(1) + (dble(i)+half)*dx(1)


            !print *, "Hack: x = ", x, " y = ", y

           if (smoothing_width .eq. 0) then

              ! discontinuous interface
              if (y .lt. l1) then
                 c(i,j,1:nspecies) = c_init_1(1:nspecies)
              elseif(y .lt. l2) then
                 c(i,j,1:nspecies) = c_init_2(1:nspecies)
              else
                 c(i,j,1:nspecies) = c_init_1(1:nspecies)
              endif

              ! print *, "Hack: c = ", c(i,j,1:nspecies)
              ! print *, "Hack: r = ", r, "Hack: r_ad = ", rad

           else

              ! smooth interface
              c(i,j,1:nspecies-1) = c_init_1(1:nspecies-1) + &
                   (c_init_2(1:nspecies-1) - c_init_1(1:nspecies-1))* &
                   (1/(1+Exp(-smoothing_width*(y-l1))) - 1/(1+Exp(-smoothing_width*(y-l2))))

!              if((j .lt. 16) .or. (j .gt. n_cells(2)-17) ) then
!                c(i,j,1:nspecies-1) = c_init_1(1:nspecies-1)
!              endif
              ! print *, "Hack: c = ", c(i,j,1:nspecies-1)

           end if

        end do
        !print *, "c: ", c(i,j,1:nspecies)
     end do
     !$omp end parallel do

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
                             dx, prob_lo, prob_hi, &
                             reallo, realhi) bind(C, name="init_rho_and_umac")

  use amrex_error_module
  use amrex_fort_module, only : amrex_real
  use common_namelist_module, only: prob_type, algorithm_type, T_init, &
                                    rho0, rhobar, initial_variance_mass, &
                                    smoothing_width
  use multispec_namelist_module, only: Dbar, Dtherm, &
                                       c_init, temp_type, is_ideal_mixture

  implicit none

  integer         , intent(in   ) :: lo(3), hi(3), nspecies
  integer         , intent(in   ) :: rlo(3),rhi(3)
  integer         , intent(in   ) :: ulo(3),uhi(3), vlo(3),vhi(3), wlo(3),whi(3)
  double precision, intent(in   ) :: reallo(3), realhi(3)
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
  double precision :: half,x,y,z,rad,L(3),sumtot,c_loc,y1,r,r1,r2,m_e

  double precision :: random

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
