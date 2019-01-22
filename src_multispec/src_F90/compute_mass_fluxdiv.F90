module compute_mass_fluxdiv_module

  use multifab_module
  use define_bc_module
  use bc_module
  use div_and_grad_module
  use diffusive_mass_fluxdiv_module
  use stochastic_mass_fluxdiv_module
  use compute_mixture_properties_module
  use external_force_module
  use ml_layout_module
  use mass_flux_utilities_module
  use convert_stag_module
  use electrodiffusive_mass_fluxdiv_module
  use probin_common_module, only: variance_coef_mass, nspecies, molmass, &
                                  shift_cc_to_boundary, k_B
  use probin_charged_module, only: use_charged_fluid, charge_per_mass

  implicit none

  private

  public :: compute_mass_fluxdiv

contains

  ! compute diffusive and stochastic mass fluxes
  ! includes barodiffusion and thermodiffusion
  subroutine compute_mass_fluxdiv(mla,rho,rhotot,gradp_baro,Temp, &
                                  diff_mass_fluxdiv,stoch_mass_fluxdiv, &
                                  diff_mass_flux,stoch_mass_flux, &
                                  dt,stage_time,dx,weights,the_bc_tower, &
                                  charge,grad_Epot,Epot,permittivity, &
                                  zero_initial_Epot_in)
    ! Donev: Add a logical flag for whether to do electroneutral or not
       
    type(ml_layout), intent(in   )   :: mla
    type(multifab) , intent(inout)   :: rho(:)
    type(multifab) , intent(inout)   :: rhotot(:)
    type(multifab) , intent(in   )   :: gradp_baro(:,:)
    type(multifab) , intent(in   )   :: Temp(:)
    type(multifab) , intent(inout)   :: diff_mass_fluxdiv(:)
    type(multifab) , intent(inout)   :: stoch_mass_fluxdiv(:)
    type(multifab) , intent(inout)   :: diff_mass_flux(:,:)
    type(multifab) , intent(inout)   :: stoch_mass_flux(:,:)
    real(kind=dp_t), intent(in   )   :: dt
    real(kind=dp_t), intent(in   )   :: stage_time 
    real(kind=dp_t), intent(in   )   :: dx(:,:)
    real(kind=dp_t), intent(in   )   :: weights(:) 
    type(bc_tower) , intent(in   )   :: the_bc_tower
    type(multifab) , intent(inout), optional :: charge(:)
    type(multifab) , intent(inout), optional :: grad_Epot(:,:)
    type(multifab) , intent(inout), optional :: Epot(:)
    type(multifab) , intent(in   ), optional :: permittivity(:)
    logical        , intent(in   ), optional :: zero_initial_Epot_in

    ! local variables
    type(multifab) :: rhoWchi(mla%nlevel)        ! rho*W*chi*Gama
    type(multifab) :: molarconc(mla%nlevel)      ! molar concentration
    type(multifab) :: molmtot(mla%nlevel)        ! total molar mass
    type(multifab) :: Hessian(mla%nlevel)        ! Hessian-matrix
    type(multifab) :: Gama(mla%nlevel)           ! Gama-matrix
    type(multifab) :: D_bar(mla%nlevel)          ! D_bar-matrix
    type(multifab) :: D_therm(mla%nlevel)        ! DT-matrix
    type(multifab) :: zeta_by_Temp(mla%nlevel)   ! for Thermo-diffusion 
    type(multifab) :: sqrtLonsager_fc(mla%nlevel,mla%dim) ! cholesky factored Lonsager on face

    integer         :: n,i,dm,nlevs

    logical :: zero_initial_Epot

    type(bl_prof_timer), save :: bpt

    call build(bpt,"compute_mass_fluxdiv")

    nlevs = mla%nlevel  ! number of levels 
    dm    = mla%dim     ! dimensionality

    zero_initial_Epot = .true.
    if (present(zero_initial_Epot_in)) then
       zero_initial_Epot = zero_initial_Epot_in
    end if
      
    ! build cell-centered multifabs for nspecies and ghost cells contained in rho.
    do n=1,nlevs
       call multifab_build(rhoWchi(n),      mla%la(n), nspecies**2, rho(n)%ng)
       call multifab_build(molarconc(n),    mla%la(n), nspecies,    rho(n)%ng)
       call multifab_build(molmtot(n),      mla%la(n), 1,           rho(n)%ng)
       call multifab_build(Hessian(n),      mla%la(n), nspecies**2, rho(n)%ng)
       call multifab_build(Gama(n),         mla%la(n), nspecies**2, rho(n)%ng)
       call multifab_build(D_bar(n),        mla%la(n), nspecies**2, rho(n)%ng)
       call multifab_build(D_therm(n),      mla%la(n), nspecies,    rho(n)%ng)
       call multifab_build(zeta_by_Temp(n), mla%la(n), nspecies,    rho(n)%ng)
       do i=1,dm
          call multifab_build_edge(sqrtLonsager_fc(n,i), mla%la(n), nspecies**2, 0, i)
       end do
    end do
 
    call compute_rhotot(mla,rho,rhotot,ghost_cells_in=.true.)
 
    ! compute molmtot, molarconc (primitive variables) for 
    ! each-cell from rho(conserved) 
    call compute_molconc_molmtot(mla,rho,rhotot,molarconc,molmtot)
      
    ! populate D_bar and Hessian matrix 
    call compute_mixture_properties(mla,rho,rhotot,D_bar,D_therm,Hessian)

    ! compute Gama from Hessian
    call compute_Gama(mla,molarconc,Hessian,Gama)
   
    ! compute rho*W*chi and zeta/Temp
    call compute_rhoWchi(mla,rho,rhotot,molarconc,rhoWchi,D_bar)
    call compute_zeta_by_Temp(mla,molarconc,D_bar,D_therm,Temp,zeta_by_Temp)

    ! compute diffusive mass fluxes, "-F = rho*W*chi*Gamma*grad(x) - ..."
    call diffusive_mass_fluxdiv(mla,rho,rhotot,molarconc,rhoWchi,Gama, &
                                diff_mass_fluxdiv,Temp,zeta_by_Temp,gradp_baro, &
                                diff_mass_flux,dx,the_bc_tower)

    ! compute external forcing for manufactured solution and add to diff_mass_fluxdiv
    ! we should move this to occur before the call to compute_mass_fluxdiv and into
    ! the advance_timestep routines
    call external_source(mla,rho,diff_mass_fluxdiv,dx,stage_time)

    ! compute stochastic fluxdiv 
    if (variance_coef_mass .ne. 0.d0) then

       ! compute face-centered cholesky-factored Lonsager^(1/2)
       call compute_sqrtLonsager_fc(mla,rho,rhotot,sqrtLonsager_fc,dx)

       call stochastic_mass_fluxdiv(mla,rho,rhotot, &
                                    sqrtLonsager_fc,stoch_mass_fluxdiv,stoch_mass_flux,&
                                    dx,dt,weights,the_bc_tower%bc_tower_array)

    end if

    if (use_charged_fluid) then
       ! we pass in these multifabs for explicit electrodiffusion
       ! implicit electrodiffusion is handled elsewhere
       if (present(charge) .and. &
           present(grad_Epot) .and. &
           present(Epot) .and. &
           present(permittivity)) then
          call electrodiffusive_mass_fluxdiv(mla,rho,Temp,rhoWchi, &
                                             diff_mass_flux,diff_mass_fluxdiv, &
                                             stoch_mass_flux, &
                                             dx,the_bc_tower, &
                                             charge,grad_Epot,Epot, &
                                             permittivity,dt,zero_initial_Epot)
       end if
    end if

    if (any(shift_cc_to_boundary(:,:) .eq. 1)) then

       ! need d+ = gamma*(grad w+)_n + a grad(phi)_n
       ! need sum (a_i grad(phi)_n)
       ! need each element of chi
       ! in each boundary cell solve the linear system
       ! overwrite total fluxes on boundary faces
       ! update diff_mass_fluxdiv near boundary faces
       call mixed_boundary_flux(mla,rho,rhotot,grad_Epot, &
                                diff_mass_flux,diff_mass_fluxdiv,Temp,dx)

    end if

    ! free the multifab allocated memory
    do n=1,nlevs
       call multifab_destroy(rhoWchi(n))
       call multifab_destroy(molarconc(n))
       call multifab_destroy(molmtot(n))
       call multifab_destroy(Hessian(n))
       call multifab_destroy(Gama(n))
       call multifab_destroy(D_bar(n))
       call multifab_destroy(D_therm(n))
       call multifab_destroy(zeta_by_Temp(n))
       do i=1,dm
          call multifab_destroy(sqrtLonsager_fc(n,i))
       end do
    end do

    call destroy(bpt)

  end subroutine compute_mass_fluxdiv

  ! need d+ = gamma*(grad w+)_n + a grad(phi)_n
  ! need sum (a_i grad(phi)_n)
  ! need each element of chi
  ! in each boundary cell solve the linear system
  ! overwrite total fluxes on boundary faces
  ! update diff_mass_fluxdiv near boundary faces
  subroutine mixed_boundary_flux(mla,rho,rhotot,grad_Epot, &
                                 diff_mass_flux,diff_mass_fluxdiv,Temp,dx)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: rho(:)
    type(multifab) , intent(inout) :: rhotot(:)
    type(multifab) , intent(in   ) :: grad_Epot(:,:)
    type(multifab) , intent(inout) :: diff_mass_flux(:,:)
    type(multifab) , intent(inout) :: diff_mass_fluxdiv(:)
    type(multifab) , intent(in   ) :: Temp(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)

    ! local variables
    integer :: i,dm,n,nlevs
    integer :: ng_1,ng_2,ng_3,ng_4,ng_5,ng_6
    integer :: lo(mla%dim),hi(mla%dim)

    ! pointers into multifabs
    real(kind=dp_t), pointer :: dp1(:,:,:,:)
    real(kind=dp_t), pointer :: dp2(:,:,:,:)
    real(kind=dp_t), pointer :: dp3x(:,:,:,:)
    real(kind=dp_t), pointer :: dp3y(:,:,:,:)
    real(kind=dp_t), pointer :: dp3z(:,:,:,:)
    real(kind=dp_t), pointer :: dp4x(:,:,:,:)
    real(kind=dp_t), pointer :: dp4y(:,:,:,:)
    real(kind=dp_t), pointer :: dp4z(:,:,:,:)
    real(kind=dp_t), pointer :: dp5(:,:,:,:)
    real(kind=dp_t), pointer :: dp6(:,:,:,:)

    dm = mla%dim
    nlevs = mla%nlevel

    ng_1 = rho(1)%ng
    ng_2 = rhotot(1)%ng
    ng_3 = grad_Epot(1,1)%ng
    ng_4 = diff_mass_flux(1,1)%ng
    ng_5 = diff_mass_fluxdiv(1)%ng
    ng_6 = Temp(1)%ng

    do n=1,nlevs
       do i=1,nfabs(rho(n))

          dp1  => dataptr(rho(n),i)
          dp2  => dataptr(rhotot(n),i)
          dp3x => dataptr(grad_Epot(n,1),i)
          dp3y => dataptr(grad_Epot(n,2),i)
          dp4x => dataptr(diff_mass_flux(n,1),i)
          dp4y => dataptr(diff_mass_flux(n,2),i)
          dp5  => dataptr(diff_mass_fluxdiv(n),i)
          dp6  => dataptr(Temp(n),i)
          lo = lwb(get_box(rho(n),i))
          hi = upb(get_box(rho(n),i))
          select case (dm)
          case (2)
             call mixed_boundary_flux_2d(dp1(:,:,1,:),ng_1, &
                                         dp2(:,:,1,1),ng_2, &
                                         dp3x(:,:,1,1),dp3y(:,:,1,1),ng_3, &
                                         dp4x(:,:,1,:),dp4y(:,:,1,:),ng_4, &
                                         dp5(:,:,1,:),ng_5, &
                                         dp6(:,:,1,1),ng_6,lo,hi,dx(n,:))
          case (3)

          end select
       end do
    end do


  end subroutine mixed_boundary_flux

  subroutine mixed_boundary_flux_2d(rho,ng_1,rhotot,ng_2,grad_Epotx,grad_Epoty,ng_3, &
                                    fluxx,fluxy,ng_4,fluxdiv,ng_5,Temp,ng_6,lo,hi,dx)

    integer         :: lo(:),hi(:),ng_1,ng_2,ng_3,ng_4,ng_5,ng_6
    real(kind=dp_t) ::        rho(lo(1)-ng_1:,lo(2)-ng_1:,:)
    real(kind=dp_t) ::     rhotot(lo(1)-ng_2:,lo(2)-ng_2:)
    real(kind=dp_t) :: grad_Epotx(lo(1)-ng_3:,lo(2)-ng_3:)
    real(kind=dp_t) :: grad_Epoty(lo(1)-ng_3:,lo(2)-ng_3:)
    real(kind=dp_t) ::      fluxx(lo(1)-ng_4:,lo(2)-ng_4:,:)
    real(kind=dp_t) ::      fluxy(lo(1)-ng_4:,lo(2)-ng_4:,:)
    real(kind=dp_t) ::    fluxdiv(lo(1)-ng_5:,lo(2)-ng_5:,:)
    real(kind=dp_t) ::       Temp(lo(1)-ng_6:,lo(2)-ng_6:)
    real(kind=dp_t) :: dx(2)

    ! local
    integer :: i,j
    integer :: comp
    
    real(kind=dp_t) :: molmtot
    real(kind=dp_t) :: molmtot_b ! boundary value
    real(kind=dp_t) :: molarconc(nspecies)
    real(kind=dp_t) :: molarconc_b(nspecies) ! boundary value
    real(kind=dp_t) :: chi(nspecies,nspecies)
    real(kind=dp_t) :: rhoWchi(nspecies,nspecies)
    real(kind=dp_t) :: D_bar(nspecies,nspecies)
    real(kind=dp_t) :: charge_coef(nspecies)
    real(kind=dp_t) :: n

    real(kind=dp_t) :: sumai, gradx(nspecies), W(nspecies)
    real(kind=dp_t) :: A(2,2), b(2), detinv, dpls, dmin, dzero

    ! hack - need to generalize this later to use shift_cc_to_boundary

    if (lo(2) .eq. 0) then
       
       j = 0
       do i=lo(1),hi(1)

          ! D_bar in valid cell
          call compute_D_bar_local(rho(i,j,:),rhotot(i,j),D_bar)

          ! molarconc and molm in valid cell
          call compute_molconc_molmtot_local(nspecies,molmass,rho(i,j,:),rhotot(i,j), &
                                             molarconc,molmtot)

          ! molarconc and molm in ghost cell
          call compute_molconc_molmtot_local(nspecies,molmass,rho(i,j-1,:),rhotot(i,j-1), &
                                             molarconc_b,molmtot_b)

      
          call compute_chi(nspecies,molmass,rho(i,j,:),rhotot(i,j),molarconc,chi,D_bar)

          call compute_rhoWchi_from_chi_local(rho(i,j,:),chi,rhoWchi)

          n = 0.d0
          do comp=1,nspecies
             n = n + rho(i,j,comp)/molmass(comp)
          end do
            
          sumai = 0.d0
          do comp=1,nspecies
             charge_coef(comp) = rho(i,j,comp)*charge_per_mass(comp)/(n*k_B*Temp(i,j))
             sumai = sumai + charge_coef(comp)*grad_Epoty(i,j)
             
             gradx(comp) = (molarconc(comp) - molarconc_b(comp))/(0.5d0*dx(2))
             W(comp) = rho(i,j,comp) / rhotot(i,j)
          end do

          dpls = gradx(1) + charge_coef(1)*grad_Epoty(i,j)

          ! form 2x2 system
          A(1,1) = chi(2,2)
          A(1,2) = chi(2,3)
          A(2,1) = 1
          A(2,2) = 1
          b(1) = -chi(2,1)*dpls
          b(2) = sumai - dpls

          ! solve 2x2 system
          detinv = 1.d0 / (A(1,1)*A(2,2) - A(1,2)*A(2,1))
          dmin = detinv * (A(2,2)*b(1) - A(1,2)*b(2))
          dzero = detinv * (-A(2,1)*b(1) + A(1,1)*b(2))

          ! update fluxes on face
          fluxy(i,j,1) = W(1)*(chi(1,1)*dpls + chi(1,2)*dmin + chi(1,3)*dzero)
          fluxy(i,j,2) = 0.d0
          fluxy(i,j,3) = -fluxy(i,j,1)

          ! update flux divergence on interior
          do comp=1,nspecies
             fluxdiv(i,j,comp) =   (fluxx(i+1,j,comp) - fluxx(i,j,comp)) / dx(1) &
                                 + (fluxy(i,j+1,comp) - fluxy(i,j,comp)) / dx(2)
          end do

       end do

    end if



  end subroutine mixed_boundary_flux_2d
  
end module compute_mass_fluxdiv_module
