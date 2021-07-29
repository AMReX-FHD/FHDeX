module compute_mixture_properties_module
  ! The purpose of the fluid model is to provide concentration-dependent transport coefficients
  ! which should change with different problem types and number of species
  ! depending on the exact physical system being simulated. 

  use amrex_error_module
  use common_namelist_module
  use multispec_namelist_module
  use matrix_utilities_module
 
  implicit none

  private

  public :: compute_D_bar_local
  
  ! The values for mixture_type distinguished at present have a constant H matrix for the thermodynamics and:
  ! 0 - no dependence of transport coefficients on composition
  ! 1 - binary mixture of water and glycerol
  ! 2 - dilute binary electrolyte solution (3 species) with solvent as last species
  !     including counter-ion cross-diffusion coefficient ~sqrt(w)  
  ! 3 - density dependent viscosity (used for bubble regression tests)
  
contains
  
  subroutine mixture_properties_mass(tlo, thi, &
                                     rho, rhlo, rhhi, &
                                     rhotot, rtlo, rthi, &
                                     D_bar, dblo, dbhi, &
                                     D_therm, dtlo, dthi, &
                                     Hessian, hslo, hshi) bind(C,name="mixture_properties_mass")

    integer         , intent(in   ) :: tlo(3),thi(3)
    integer         , intent(in   ) :: rhlo(3),rhhi(3), rtlo(3),rthi(3), dblo(3),dbhi(3), dtlo(3),dthi(3), hslo(3),hshi(3)
    double precision, intent(in   ) ::     rho(rhlo(1):rhhi(1),rhlo(2):rhhi(2),rhlo(3):rhhi(3),nspecies)     ! density
    double precision, intent(in   ) ::  rhotot(rtlo(1):rhhi(1),rtlo(2):rthi(2),rtlo(3):rthi(3))              ! total density
    double precision, intent(inout) ::   D_bar(dblo(1):dbhi(1),dblo(2):dbhi(2),dblo(3):dbhi(3),nspecies*nspecies)
    double precision, intent(inout) :: D_therm(dblo(1):dbhi(1),dblo(2):dbhi(2),dblo(3):dbhi(3),nspecies)
    double precision, intent(inout) :: Hessian(hslo(1):hshi(1),hslo(2):hshi(2),hslo(3):hshi(3),nspecies*nspecies)
    
    ! local varialbes
    integer :: i,j,k

    ! for specific box, now start loops over alloted cells 
    do k=tlo(3),thi(3)
       do j=tlo(2),thi(2)
          do i=tlo(1),thi(1)

             call mixture_properties_mass_local(rho(i,j,k,:),rhotot(i,j,k), &
                                                D_bar(i,j,k,:),D_therm(i,j,k,:),&
                                                Hessian(i,j,k,:))
          end do
       end do
    end do
   
  end subroutine mixture_properties_mass

  ! The default case should be to simply set D_bar, D_therm and H to constants, read from the input file
  subroutine mixture_properties_mass_local(rho,rhotot,D_bar,D_therm,Hessian)
   
    double precision, intent(in)   :: rho(nspecies)        
    double precision, intent(in)   :: rhotot
    double precision, intent(out)  :: D_bar(nspecies,nspecies) 
    double precision, intent(out)  :: D_therm(nspecies) 
    double precision, intent(out)  :: Hessian(nspecies,nspecies)
 
    ! local variables
    integer :: n,row,column
    double precision :: massfrac(nspecies)

    ! Local values of transport and thermodynamic coefficients (functions of composition!):
    double precision, dimension(nspecies*(nspecies-1)/2) :: H_offdiag_local ! off-diagonal components of symmetric matrices
    double precision, dimension(nspecies) :: H_diag_local ! Diagonal component
    
    massfrac = rho/rhotot;

    ! populate D_bar
    call compute_D_bar_local(rho,rhotot,D_bar)
    
    ! For now we only encode constant Hessian matrices since we do not have any thermodynamic 
    ! models coded up (Wilson, NTLR, UNIQUAC, etc.)
    H_diag_local(1:nspecies) = H_diag(1:nspecies)
    H_offdiag_local(1:nspecies*(nspecies-1)/2) = H_offdiag(1:nspecies*(nspecies-1)/2)
    
    ! Complete the process by filling the matrices using generic formulae -- this part should not change
    ! populate Hessian matrix 
    n=0; 
    do row=1, nspecies  
       do column=1, row-1
          n=n+1
          if(is_ideal_mixture .eq. 0) then
             Hessian(row, column) = H_offdiag_local(n)   ! positive semidefinite matrix read from input
             Hessian(column, row) = Hessian(row,column)  ! Hessian is symmetric
          end if
       end do
       
       ! populate diagonals 
       D_therm(row)    = Dtherm(row)    ! thermal diffcoeff's read from input    
       if(is_ideal_mixture .eq. 0) then 
          Hessian(row, row) = H_diag_local(row)     
       else 
          Hessian = 0.d0     ! set matrix to zero for ideal-mixture
       end if
    
    end do
    
  end subroutine mixture_properties_mass_local

  subroutine compute_D_bar_local(rho,rhotot,D_bar)

    double precision, intent(in)   :: rho(nspecies)        
    double precision, intent(in)   :: rhotot
    double precision, intent(out)  :: D_bar(nspecies,nspecies) 

    ! off-diagonal components of symmetric matrices
    double precision, dimension(nspecies*(nspecies-1)/2) :: D_bar_local
    
    integer :: n,row,column


    select case (abs(mixture_type))
    case (1)  ! water-glycerol
    
       ! we require nspecies=2
       ! Dbar(1) = chi0 in the binary notation
       if (nspecies .ne. 2) then
          call bl_error("mixture_properties_mass_local assumes nspecies=2 if mixture_type=3 (water-glycerol)")
       end if
       
       call chi_water_glycerol(D_bar_local(1), rho, rhotot)

    case (2) ! Electrolyte mixture
    
       if (nspecies .ne. 3) then
          call bl_error("mixture_properties_mass_local assumes nspecies=3 if mixture_type=2 (water-glycerol)")
       end if

       ! This is the leading-order correction for dilute solutions
       ! In particular, the cross-diffusion coefficient ~sqrt(concentration) as per renormalization theory
       ! The ordering of the values is D_12; D_13, D_23

       D_bar_local(1) = Dbar(1)*sqrt(rho(1)/rhotot) ! counter-ion cross coefficient D_12~sqrt(w)
       D_bar_local(2) = Dbar(2) ! D_13 = self diffusion of first ion
       D_bar_local(3) = Dbar(3) ! D_23 = self diffusion of second ion


    case default

       D_bar_local(1:nspecies*(nspecies-1)/2) = Dbar(1:nspecies*(nspecies-1)/2) ! Keep it constant


    end select

    
    ! Complete the process by filling the matrices using generic formulae -- this part should not change
    ! populate D_bar and Hessian matrix 
    n=0; 
    do row=1, nspecies  
       do column=1, row-1
          n=n+1
          D_bar(row, column) = D_bar_local(n)                ! SM-diffcoeff's read from input
          D_bar(column, row) = D_bar(row, column)     ! symmetric
       end do       
       ! populate diagonals 
       D_bar(row, row) = 0.d0           ! as self-diffusion is zero
    end do


  end subroutine compute_D_bar_local

  ! subroutine compute_eta_kappa(mla,eta,eta_ed,kappa,rho,rhotot,Temp,dx,the_bc_level)

  !   type(ml_layout), intent(in   ) :: mla
  !   type(multifab) , intent(inout) :: eta(:)
  !   type(multifab) , intent(inout) :: eta_ed(:,:) ! nodal (2d); edge-centered (3d)
  !   type(multifab) , intent(inout) :: kappa(:)
  !   type(multifab) , intent(in   ) :: rho(:)
  !   type(multifab) , intent(in   ) :: rhotot(:)
  !   type(multifab) , intent(in   ) :: Temp(:)
  !   double precision, intent(in   ) :: dx(:,:)
  !   type(bc_level) , intent(in   ) :: the_bc_level(:)

  !   integer :: nlevs,dm,i,n,ng_e,ng_r,ng_m,ng_t
  !   integer :: lo(mla%dim),hi(mla%dim)

  !   double precision, pointer :: ep(:,:,:,:)
  !   double precision, pointer :: rp(:,:,:,:)
  !   double precision, pointer :: mp(:,:,:,:)
  !   double precision, pointer :: tp(:,:,:,:)

  !   type(bl_prof_timer), save :: bpt

  !   call build(bpt,"compute_eta_kappa")

  !   ng_e  = eta(1)%ng
  !   ng_r  = rho(1)%ng
  !   ng_m = rhotot(1)%ng
  !   ng_t  = Temp(1)%ng

  !   nlevs = mla%nlevel
  !   dm = mla%dim

  !   do n=1,nlevs
  !      call multifab_setval(kappa(n), 1.d0, all=.true.)
  !   end do

  !   do n=1,nlevs
  !      do i=1,nfabs(eta(n))
  !         ep  => dataptr(eta(n), i)
  !         rp  => dataptr(rho(n), i)
  !         mp => dataptr(rhotot(n), i)
  !         tp  => dataptr(Temp(n), i)
  !         lo = lwb(get_box(eta(n), i))
  !         hi = upb(get_box(eta(n), i))
  !         select case (dm)
  !         case (2)
  !            call compute_eta_2d(ep(:,:,1,1),ng_e,rp(:,:,1,:),ng_r,mp(:,:,1,1),ng_m, &
  !                                tp(:,:,1,1),ng_t,lo,hi,dx(n,:))
  !         case (3)
  !            call compute_eta_3d(ep(:,:,:,1),ng_e,rp(:,:,:,:),ng_r,mp(:,:,:,1),ng_m, &
  !                                tp(:,:,:,1),ng_t,lo,hi,dx(n,:))
  !         end select
  !      end do
  !   end do

  !   if (dm .eq. 2) then
  !      call average_cc_to_node(nlevs,eta,eta_ed(:,1),1,tran_bc_comp,1,the_bc_level)
  !   else if (dm .eq. 3) then
  !      call average_cc_to_edge(nlevs,eta,eta_ed,1,tran_bc_comp,1,the_bc_level)
  !   end if

  !   call destroy(bpt)

  ! end subroutine compute_eta_kappa

  ! subroutine compute_eta_2d(eta,ng_e,rho,ng_r,rhotot,ng_m,Temp,ng_t,lo,hi,dx)

  !   ! compute eta in valid AND ghost regions
  !   ! the ghost cells for rho, Temp, etc., have already been filled properly

  !   integer        , intent(in   ) :: lo(:), hi(:), ng_e, ng_r, ng_m, ng_t
  !   double precision, intent(inout) ::    eta(lo(1)-ng_e:,lo(2)-ng_e:)
  !   double precision, intent(inout) ::    rho(lo(1)-ng_r:,lo(2)-ng_r:,:)
  !   double precision, intent(inout) :: rhotot(lo(1)-ng_m:,lo(2)-ng_m:)
  !   double precision, intent(inout) ::   Temp(lo(1)-ng_t:,lo(2)-ng_t:)
  !   double precision, intent(in   ) :: dx(:)

  !   ! local
  !   integer :: i,j

  !   select case (mixture_type)
  !   case (1)
       
  !      do j=lo(2)-ng_e,hi(2)+ng_e
  !         do i=lo(1)-ng_e,hi(1)+ng_e

  !            call eta_water_glycerol(eta(i,j),rho(i,j,:),rhotot(i,j),Temp(i,j))

  !         end do
  !      end do

  !   case (3)

  !      do j=lo(2)-ng_e,hi(2)+ng_e
  !         do i=lo(1)-ng_e,hi(1)+ng_e

  !            eta(i,j) = visc_coef*(-5.d0 + 6.d0*rhotot(i,j))
                
  !         end do
  !      end do

  !   case default

  !      eta = visc_coef

  !   end select

  ! end subroutine compute_eta_2d

  ! subroutine compute_eta_3d(eta,ng_e,rho,ng_r,rhotot,ng_m,Temp,ng_t,lo,hi,dx)

  !   ! compute eta in valid AND ghost regions
  !   ! the ghost cells for rho, Temp, etc., have already been filled properly

  !   integer        , intent(in   ) :: lo(:), hi(:), ng_e, ng_r, ng_m, ng_t
  !   double precision, intent(inout) ::    eta(lo(1)-ng_e:,lo(2)-ng_e:,lo(3)-ng_e:)
  !   double precision, intent(inout) ::    rho(lo(1)-ng_r:,lo(2)-ng_r:,lo(3)-ng_r:,:)
  !   double precision, intent(inout) :: rhotot(lo(1)-ng_m:,lo(2)-ng_m:,lo(3)-ng_m:)
  !   double precision, intent(inout) ::   Temp(lo(1)-ng_t:,lo(2)-ng_t:,lo(3)-ng_t:)
  !   double precision, intent(in   ) :: dx(:)

  !   ! local
  !   integer :: i,j,k

  !   select case (mixture_type)
  !   case (1)

  !      do k=lo(3)-ng_e,hi(3)+ng_e
  !         do j=lo(2)-ng_e,hi(2)+ng_e
  !            do i=lo(1)-ng_e,hi(1)+ng_e

  !               call eta_water_glycerol(eta(i,j,k),rho(i,j,k,:),rhotot(i,j,k),Temp(i,j,k))

  !            end do
  !         end do
  !      end do

  !      case (3)

  !         do k=lo(3)-ng_e,hi(3)+ng_e
  !            do j=lo(2)-ng_e,hi(2)+ng_e
  !               do i=lo(1)-ng_e,hi(1)+ng_e

  !                  eta(i,j,k) = visc_coef*(-5.d0 + 6.d0*rhotot(i,j,k))
                
  !               end do
  !            end do
  !         end do

  !   case default

  !      eta = visc_coef

  !   end select

  ! end subroutine compute_eta_3d

  !=================================================
  ! Water-glycerol mixtures near room temperature
  !=================================================
  
  subroutine chi_water_glycerol(chi,rho,rhotot)

    ! This only works for room temperature for now
    double precision, intent(inout) :: chi
    double precision, intent(in   ) :: rho(:)
    double precision, intent(in   ) :: rhotot

    ! local
    double precision :: c_loc

    ! mass fraction of glycerol
    c_loc = rho(1)/rhotot

    ! chi = chi0 * rational function
    chi = Dbar(1)*(1.024d0-1.001692692d0*c_loc)/(1.d0+0.6632641981d0*c_loc)
        
  end subroutine chi_water_glycerol

  ! subroutine eta_water_glycerol(eta,rho,rhotot,Temp)

  !   double precision, intent(inout) :: eta
  !   double precision, intent(in   ) :: rho(:)
  !   double precision, intent(in   ) :: rhotot
  !   double precision, intent(in   ) :: Temp
    
  !   ! local
  !   double precision :: c_loc, nu_g, nu_w, x, T

  !   type(bl_prof_timer), save :: bpt

  !   call build(bpt,"eta_water_glycerol")

  !   ! mass fraction of glycerol
  !   c_loc = rho(1)/rhotot

  !   ! convert temperature to Celsius
  !   T = Temp - 273.d0 

  !   ! viscosities of pure glycerol and water
  !   nu_g = exp(9.09309 - 0.11397*T + 0.00054*T**2)
  !   nu_w = exp(0.55908 - 0.03051*T + 0.00015*T**2)

  !   x = c_loc*(1.d0 + (1.d0-c_loc)*(-0.727770 - 0.04943*c_loc - 1.2038*c_loc**2))

  !   ! visc_coef is the scaling prefactor (for unit conversion or non-unity scaling)
  !   eta = visc_coef*exp(-x*(-log(nu_g)+log(nu_w)))*nu_w*rhotot

  !   call destroy(bpt)

  ! end subroutine eta_water_glycerol

end module compute_mixture_properties_module
