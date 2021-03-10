module mass_flux_utilities_module 

  use amrex_error_module
  use common_namelist_module
  use multispec_namelist_module
  use matrix_utilities_module
  use compute_mixture_properties_module

  implicit none

  private

contains

  subroutine compute_molconc_molmtot(tlo, thi, &
                                     rho, rhlo, rhhi, &
                                     rhotot, rtlo, rthi, &
                                     molarconc, mclo, mchi, &
                                     molmtot, mtlo, mthi) bind(C,name="compute_molconc_molmtot")
 
    integer         , intent(in   ) :: tlo(3),thi(3)
    integer         , intent(in   ) :: rhlo(3),rhhi(3), rtlo(3),rthi(3), mclo(3),mchi(3), mtlo(3),mthi(3)
    double precision, intent(in   ) ::       rho(rhlo(1):rhhi(1),rhlo(2):rhhi(2),rhlo(3):rhhi(3),nspecies) ! density
    double precision, intent(in   ) ::    rhotot(rtlo(1):rthi(1),rtlo(2):rthi(2),rtlo(3):rthi(3))          ! total density in each cell 
    double precision, intent(inout) :: molarconc(mclo(1):mchi(1),mclo(2):mchi(2),mclo(3):mchi(3),nspecies) ! molar concentration
    double precision, intent(inout) ::   molmtot(mtlo(1):mthi(1),mtlo(2):mthi(2),mtlo(3):mthi(3))          ! total molar mass 
    
    ! local variables
    integer          :: i,j,k
    
    ! for specific box, now start loops over alloted cells    
    do k=tlo(3), thi(3)
       do j=tlo(2), thi(2)
          do i=tlo(1), thi(1)

             call compute_molconc_molmtot_local(nspecies,molmass,rho(i,j,k,:),rhotot(i,j,k),&
                                                molarconc(i,j,k,:),molmtot(i,j,k))

          end do
       end do
    end do
 
  end subroutine compute_molconc_molmtot
  
  subroutine compute_molconc_molmtot_local(nspecies_in,molmass_in,rho,rhotot,molarconc,molmtot)

    integer,          intent(in)   :: nspecies_in
    double precision, intent(in)   :: molmass_in(nspecies_in)
    double precision, intent(in)   :: rho(nspecies_in)           ! density- last dim for #species
    double precision, intent(in)   :: rhotot                     ! total density in each cell
    double precision, intent(out)  :: molarconc(nspecies_in)     ! molar concentration
    double precision, intent(out)  :: molmtot                    ! total molar mass
    
    ! local variables
    integer          :: n
    double precision, dimension(nspecies_in) :: W                ! mass fraction w_i = rho_i/rho
    double precision  :: Sum_woverm

    ! calculate mass fraction and total molar mass (1/m=Sum(w_i/m_i))
    Sum_woverm=0.d0
    do n=1, nspecies_in
       W(n) = rho(n)/rhotot
       Sum_woverm = Sum_woverm + W(n)/molmass_in(n)
    end do
    molmtot = 1.0d0/Sum_woverm 

    ! calculate molar concentrations in each cell (x_i=m*w_i/m_i) 
    do n=1, nspecies_in
       molarconc(n) = molmtot*W(n)/molmass_in(n)
    end do
    
  end subroutine compute_molconc_molmtot_local

  subroutine compute_Gamma(tlo,thi, &
                           molarconc, mlo, mhi, &
                           Hessian, hlo, hhi, & 
                           Gamma, glo, ghi) bind(C,name="compute_Gamma")

    integer, intent(in   )          :: tlo(3),thi(3)
    integer, intent(in   )          :: mlo(3),mhi(3), hlo(3),hhi(3), glo(3),ghi(3) 
    double precision, intent(in   ) :: molarconc(mlo(1):mhi(1),mlo(2):mhi(2),mlo(3):mhi(3),nspecies) ! molar concentration 
    double precision, intent(in   ) ::   Hessian(hlo(1):hhi(1),hlo(2):hhi(2),hlo(3):hhi(3),nspecies*nspecies)
    double precision, intent(inout) ::     Gamma(glo(1):ghi(1),glo(2):ghi(2),glo(3):ghi(3),nspecies*nspecies)

    ! local varialbes
    integer          :: i,j,k

    ! for specific box, now start loops over alloted cells 
    do k=tlo(3),thi(3)
       do j=tlo(2),thi(2)
          do i=tlo(1),thi(1)

             call compute_Gamma_local(molarconc(i,j,k,:),Hessian(i,j,k,:),Gamma(i,j,k,:),nspecies)
          end do
       end do
    end do
   
  end subroutine compute_Gamma

  subroutine compute_Gamma_local(molarconc,Hessian,Gamma,nspecies)
   
    integer         , intent(in   )   :: nspecies
    double precision, intent(in   )   :: molarconc(nspecies)
    double precision, intent(in   )   :: Hessian(nspecies,nspecies)
    double precision, intent(inout)   :: Gamma(nspecies,nspecies)
 
    ! local variables
    integer                                        :: row,column
    double precision, dimension(nspecies,nspecies) :: I, X_xxT

    ! local to define Gamma matrix
    double precision :: w1,w2 

    if (use_multiphase .eq. 1 .and. nspecies .eq. 2) then
       
       ! Gamma = I
       w1 = molarconc(1)
       w2 = molarconc(2)
       if(abs(w1+w2-1.d0).gt.1.d-14)then
          write(6,*)" mole fractions dont' add up in gamma computation"
       endif
       if(w1.lt.0.d0)then
          w1 = 0.d0
          w2 = 1.d0
       endif
       if(w2.lt.0.d0)then
          w2 = 0.d0
          w1 = 1.d0
       endif

       Gamma(1,2)=w1*dfloat(n_gex*n_gex)*alpha_gex*(w1**(n_gex-1))*(w2**(n_gex-1))
       Gamma(2,1)=w2*dfloat(n_gex*n_gex)*alpha_gex*(w2**(n_gex-1))*(w1**(n_gex-1))
       Gamma(1,1)=1.d0+w1*dfloat(n_gex*(n_gex-1))*alpha_gex*(w1**(n_gex-2))*(w2**n_gex)
       Gamma(2,2)=1.d0+w2*dfloat(n_gex*(n_gex-1))*alpha_gex*(w2**(n_gex-2))*(w1**n_gex)

    else

       ! Identity matrix
       I = 0.d0
       do row=1,nspecies
          I(row, row) = 1.d0        
       end do

       ! populate X_xxT
       if (is_ideal_mixture .eq. 1) then
          X_xxT = 0.d0
       else
          do row=1, nspecies  
             ! diagonal entries
             X_xxT(row,row) = molarconc(row) - molarconc(row)**2 
             do column=1, row-1
                ! off-diagnoal entries
                X_xxT(row,column)   = -molarconc(row)*molarconc(column)  ! form x*transpose(x) off diagonals 
                X_xxT(column, row)  = X_xxT(row, column)                 ! symmetric
             end do
          end do
       end if

       ! compute Gamma 
       Gamma = I + matmul(X_xxT, Hessian)
       
    end if
 
  end subroutine compute_Gamma_local


  subroutine compute_rhoWchi( tlo, thi, &
                              rho, rhlo, rhhi, &
                              rhotot, rtlo, rthi, &
                              molarconc, mclo, mchi, &
                              rhoWchi, rclo, rchi, &
                              D_bar, dblo, dbhi) bind(C,name="compute_rhoWchi")
 
    integer         , intent(in   ) :: tlo(3),thi(3)
    integer         , intent(in   ) :: rhlo(3),rhhi(3), rtlo(3),rthi(3), mclo(3),mchi(3), rclo(3),rchi(3), dblo(3),dbhi(3)
    double precision, intent(in   ) ::       rho(rhlo(1):rhhi(1),rhlo(2):rhhi(2),rhlo(3):rhhi(3),nspecies)          ! densities
    double precision, intent(in   ) ::    rhotot(rtlo(1):rthi(1),rtlo(2):rthi(2),rtlo(3):rthi(3))                   ! total density
    double precision, intent(in   ) :: molarconc(mclo(1):mchi(1),mclo(2):mchi(2),mclo(3):mchi(3),nspecies)          ! molar concentration
    double precision, intent(inout) ::   rhoWchi(rclo(1):rchi(1),rclo(2):rchi(2),rclo(3):rchi(3),nspecies*nspecies)
    double precision, intent(in   ) ::     D_bar(dblo(1):dbhi(1),dblo(2):dbhi(2),dblo(3):dbhi(3),nspecies*nspecies) ! MS diff-coefs
    
    ! local variables
    integer :: i,j,k

    ! for specific box, now start loops over alloted cells 
    do k=tlo(3),thi(3)
       do j=tlo(2),thi(2)
          do i=tlo(1),thi(1)
       
             call compute_rhoWchi_local(rho(i,j,k,:),rhotot(i,j,k),molarconc(i,j,k,:),&
                                    rhoWchi(i,j,k,:),D_bar(i,j,k,:),nspecies)

          end do
       end do
    end do
   
  end subroutine compute_rhoWchi


  subroutine compute_rhoWchi_local(rho,rhotot,molarconc,rhoWchi,D_bar,nspecies)

    integer,          intent(in   ) :: nspecies
    double precision, intent(in   ) :: rho(nspecies)
    double precision, intent(in   ) :: rhotot
    double precision, intent(in   ) :: molarconc(nspecies)
    double precision, intent(inout) :: rhoWchi(nspecies,nspecies)
    double precision, intent(in   ) :: D_bar(nspecies,nspecies)

    ! local variables
    integer         :: row,column,k
    double precision :: W(nspecies)

    double precision :: chi(nspecies,nspecies)
    double precision :: Deff, tmp

    integer         :: ntrace           ! number of trace species with w_k < fractional_tolerance
    integer         :: nspecies_sub     ! dim of subsystem = nspecies-ntrace
    double precision :: molmtot_sub
    double precision :: rhotot_sub
    double precision :: w1,w2


    ! this is a mapping used to eliminate elements in D_bar we don't need (for D_bar_sub)
    ! and for expanding sqrtLonsager_sub into sqrtLonsager
    ! it will contain the numbers 1, 2, ..., (nspecies-ntrace)
    ! with zeros in elements corresponding to trace elements
    ! (example) for a 5-species system having trace species 2 and 5:
    !  species       1 2 3 4 5
    !  dest(species) 1 0 2 3 0
    integer :: dest(nspecies)

    if (use_multiphase .eq. 1 .and. nspecies .eq. 2) then

        w1 = molarconc(1)
        w2 = molarconc(2)
        if(w1.lt.0)then
           w1=0.d0
           w2=1.d0
        endif
        if(w2.lt.0)then
           w2=0.d0
           w1=1.d0
        endif

        rhoWchi(1,1) = w2*rho0*D_bar(1,2)
        rhoWchi(1,2) = -w1*rho0*D_bar(1,2)
        rhoWchi(2,1) = -w2*rho0*D_bar(1,2)
        rhoWchi(2,2) = w1*rho0*D_bar(1,2)

     else

       ! compute the number of trace species
       ! build the mapping for expanding/contracting arrays
       ntrace = 0
       do row=1, nspecies
          W(row) = rho(row)/rhotot
          if (W(row) .lt. fraction_tolerance) then
             ntrace = ntrace + 1
             dest(row) = 0
          else
             dest(row) = row - ntrace
          end if
       end do

       if (ntrace .eq. nspecies-1) then

          ! this is all trace species except for 1 (essentially pure solvent);
          ! set rhoWchi to zero
          rhoWchi(:,:) = 0.d0

       else if (ntrace .eq. 0) then

          ! there are no trace species
          ! hence, chi = chi_sub

          ! compute chi 
          call compute_chi(nspecies,molmass,rho,rhotot,molarconc,chi,D_bar)

          ! compute rho*W*chi
          do row=1, nspecies
             do column=1, nspecies
                rhoWchi(row,column) = rho(row)*chi(row,column)
             end do
          end do

       else

          ! if there are trace species, we consider a subsystem 
          ! consisting of non-trace species

          nspecies_sub = nspecies - ntrace
          call compute_chi_sub()

       end if

    end if
    
  contains
  
    subroutine compute_chi_sub() ! We make this a subroutine to put local arrays on the stack and not heap
       double precision :: molmass_sub(nspecies_sub)
       double precision :: rho_sub(nspecies_sub)
       double precision :: W_sub(nspecies_sub)
       double precision :: molarconc_sub(nspecies_sub)

       double precision :: D_bar_sub(nspecies_sub,nspecies_sub)
       double precision :: chi_sub(nspecies_sub,nspecies_sub)

       ! create a vector of non-trace densities and molmass for the subsystem
       do row=1, nspecies
          if (dest(row) .ne. 0) then
             molmass_sub(dest(row)) = molmass(row)
             rho_sub(dest(row)) = rho(row)
          end if
       end do

       ! renormalize total density and mass fractions
       rhotot_sub = sum(rho_sub)
       do row=1, nspecies_sub
          W_sub(row) = rho_sub(row)/rhotot_sub
       end do

       ! construct D_bar_sub by mapping the full D_bar into D_bar_sub
       ! you could read in only the lower diagonals, 
       ! reflect, and set the diagnals to zero if you want

       do row=1, nspecies
          if (dest(row) .eq. 0) then
             cycle
          end if
          do column=1, nspecies
             if (dest(column) .ne. 0) then
                D_bar_sub(dest(row),dest(column)) = D_bar(row,column)
             end if
          end do
       end do

       ! compute molarconc_sub and molmtot_sub
       call compute_molconc_molmtot_local(nspecies_sub,molmass_sub,rho_sub,rhotot_sub,molarconc_sub,molmtot_sub)

       ! compute chi_sub
       call compute_chi(nspecies_sub,molmass_sub,rho_sub,rhotot_sub,molarconc_sub,chi_sub,D_bar_sub)

       ! compute full rho*W*chi
       rhoWchi(:,:) = 0.d0
       do column=1, nspecies
          if (dest(column) .eq. 0) then         ! column of a trace species
             ! compute Deff
             Deff = 0.d0
             do k=1, nspecies
                if (dest(k) .ne. 0) then
                   Deff = Deff + molarconc_sub(dest(k))/D_bar(k,column)
                end if
             end do
             Deff = 1.d0/Deff

             ! assign rhoWchi
             do row=1, nspecies
                if (row .eq. column) then
                   rhoWchi(row,column) = rhotot_sub*Deff*molmass(row)/molmtot_sub
                else if (dest(row) .eq. 0) then
                   rhoWchi(row,column) = 0.d0
                else
                   tmp = 0.d0
                   do k=1, nspecies
                      if (dest(k) .ne. 0) then
                         tmp = tmp + chi_sub(dest(row),dest(k))*molarconc_sub(dest(k))/D_bar(k,column)
                      end if
                   end do
                   rhoWchi(row,column) = Deff*rho_sub(dest(row))*(tmp-molmass(column)/molmtot_sub)
                end if
             end do
          else                                  ! column of a non-trace species
             ! assign rhoWchi
             do row=1, nspecies
                if (dest(row) .eq. 0) then
                   rhoWchi(row,column) = 0.d0
                else
                   rhoWchi(row,column) = rho_sub(dest(row))*chi_sub(dest(row),dest(column))
                end if
             end do 
          end if
       end do
    
    end subroutine compute_chi_sub

  end subroutine compute_rhoWchi_local

  subroutine compute_chi(nspecies_in,molmass_in,rho,rhotot,molarconc,chi,D_bar)
   
    integer,          intent(in   ) :: nspecies_in
    double precision, intent(in   ) :: molmass_in(nspecies_in)
    double precision, intent(in   ) :: rho(nspecies_in)
    double precision, intent(in   ) :: rhotot
    double precision, intent(in   ) :: molarconc(nspecies_in)
    double precision, intent(inout) :: chi(nspecies_in,nspecies_in)
    double precision, intent(in   ) :: D_bar(nspecies_in,nspecies_in)

    ! local variables
    integer                         :: row,column
    double precision                :: Sum_knoti   
    double precision                :: tmp
    double precision                :: eepsilon

    ! vectors and matrices to be used by LAPACK 
    double precision, dimension(nspecies_in,nspecies_in) :: Lambda
    double precision, dimension(nspecies_in)             :: W

    eepsilon = 1.d-16
    
    ! if nspecies_in = 2, use analytic formulas
    ! note: nspecies_in = 1 (that is, ntrace = nspecies-1) is treated separately (chi=0)
    !       before this routine is called
    if (nspecies_in .eq. 2) then
       W(1) = rho(1)/rhotot
       W(2) = rho(2)/rhotot

       if (use_multiphase .eq. 1) then
          W(1) = max(min(W(1),1.d0),eepsilon)
          W(2) = max(min(W(2),1.d0),eepsilon)
       end if

       tmp = molmass_in(1)*W(2)+molmass_in(2)*W(1)
       tmp = D_bar(1,2)*tmp*tmp/molmass_in(1)/molmass_in(2)

       chi(1,1) = tmp*W(2)/W(1)
       chi(1,2) = -tmp
       chi(2,1) = -tmp
       chi(2,2) = tmp*W(1)/W(2)

       return
    end if

    ! compute chi either selecting inverse/pseudoinverse or iterative methods 
    if (use_lapack .eq. 1) then
       call amrex_error('compute_chi: use_lapack not supported')
    end if

    call Dbar2chi_iterative(nspecies_in,chi_iterations,D_bar(:,:),molarconc(:),molmass_in(:),chi(:,:))

  end subroutine compute_chi

  subroutine compute_zeta_by_Temp(tlo,thi, &
                                  molarconc,mclo,mchi, & 
                                  D_bar,dblo,dbhi, & 
                                  Temp,tplo,tphi, & 
                                  zeta_by_Temp,ztlo,zthi, &
                                  D_therm,dtlo,dthi) bind(C,name="compute_zeta_by_Temp")

    integer,          intent(in   ) :: tlo(3),thi(3)
    integer,          intent(in   ) :: mclo(3),mchi(3), dblo(3),dbhi(3), tplo(3),tphi(3), ztlo(3),zthi(3), dtlo(3),dthi(3)
    double precision, intent(in   ) ::    molarconc(mclo(1):mchi(1),mclo(2):mchi(2),mclo(3):mchi(3),nspecies) ! molar concentration 
    double precision, intent(in   ) ::        D_bar(dblo(1):dbhi(1),dblo(2):dbhi(2),dblo(3):dbhi(3),nspecies) ! MS diff-coefs 
    double precision, intent(in   ) ::         Temp(tplo(1):tphi(1),tplo(2):tphi(2),tplo(3):tphi(3))          ! Temperature 
    double precision, intent(inout) :: zeta_by_Temp(ztlo(1):zthi(1),ztlo(2):zthi(2),ztlo(3):zthi(3),nspecies) ! zeta/T
    double precision, intent(in   ) ::      D_therm(dtlo(1):dthi(1),dtlo(2):dthi(2),dtlo(3):dthi(3),nspecies) ! thermo diff-coefs 
    
    ! local variables
    integer          :: i,j,k

    ! for specific box, now start loops over alloted cells 
    do k=tlo(3),thi(3)
       do j=tlo(2),thi(2)
          do i=tlo(1),thi(1)
       
             call compute_zeta_by_Temp_local(molarconc(i,j,k,:),&
                                             D_bar(i,j,k,:),Temp(i,j,k),zeta_by_Temp(i,j,k,:),&
                                             D_therm(i,j,k,:),nspecies)

          end do
       end do
    end do
   
  end subroutine compute_zeta_by_Temp

  subroutine compute_zeta_by_Temp_local(molarconc,D_bar,Temp,zeta_by_Temp,D_therm,nspecies)
    
    integer,          intent(in   ) :: nspecies
    double precision, intent(in   ) :: molarconc(nspecies) 
    double precision, intent(in   ) :: D_bar(nspecies,nspecies) 
    double precision, intent(in   ) :: Temp
    double precision, intent(inout) :: zeta_by_Temp(nspecies)
    double precision, intent(in   ) :: D_therm(nspecies)

    ! local variables
    integer                          :: row,column
    double precision                 :: Sum_knoti   

    ! vectors and matrices to be used by LAPACK 
    double precision, dimension(nspecies,nspecies) :: Lambda

    ! compute zeta_by_Temp for thermodiffusion
    if(is_nonisothermal .eq. 1) then

       ! compute Lambda_ij matrix; molarconc is 
       ! expressed in terms of molmtot,mi,rhotot etc. 
       do row=1, nspecies  
          do column=1, row-1
             Lambda(row, column) = -molarconc(row)*molarconc(column)/D_bar(row,column)
             Lambda(column, row) = Lambda(row, column) 
          end do
       end do

       do row=1, nspecies
          Sum_knoti = 0.d0
          do column=1, nspecies
             if(column.ne.row) then
                Sum_knoti = Sum_knoti + Lambda(row,column)*(D_therm(row)-D_therm(column))
             end if
             zeta_by_Temp(row) = Sum_knoti/Temp
          end do
       end do

    end if

  end subroutine compute_zeta_by_Temp_local


  subroutine compute_sqrtLonsager_fc(lo,hi, &
                                     rho, rhlo, rhhi, &
                                     rhotot, rtlo, rthi, &
                                     sqrtLonsager_x, sxlo, sxhi, &
                                     sqrtLonsager_y, sylo, syhi, &
#if (AMREX_SPACEDIM == 3)
                                     sqrtLonsager_z, szlo, szhi, &
#endif
                                     dx) bind(C,name="compute_sqrtLonsager_fc")

    integer         , intent(in   ) :: lo(3), hi(3)
    integer         , intent(in   ) :: rhlo(3),rhhi(3), rtlo(3),rthi(3), sxlo(3),sxhi(3), sylo(3),syhi(3)
    double precision, intent(in   ) ::            rho(rhlo(1):rhhi(1),rhlo(2):rhhi(2),rhlo(3):rhhi(3),nspecies)
    double precision, intent(in   ) ::         rhotot(rtlo(1):rthi(1),rtlo(2):rthi(2),rtlo(3):rthi(3))
    double precision, intent(inout) :: sqrtLonsager_x(sxlo(1):sxhi(1),sxlo(2):sxhi(2),sxlo(3):sxhi(3),nspecies*nspecies)
    double precision, intent(inout) :: sqrtLonsager_y(sylo(1):syhi(1),sylo(2):syhi(2),sylo(3):syhi(3),nspecies*nspecies)
#if (AMREX_SPACEDIM == 3)
    integer         , intent(in   ) :: szlo(3), szhi(3)
    double precision, intent(inout) :: sqrtLonsager_z(szlo(1):szhi(1),szlo(2):szhi(2),szlo(3):szhi(3),nspecies*nspecies)
#endif
    double precision, intent(in   ) :: dx(3)

    ! local variables
    integer         :: i,j,k
    double precision :: rhoav(nspecies)
    
    ! x-faces
    do k=lo(3),hi(3)
    do j=lo(2),hi(2)
    do i=lo(1),hi(1)+1
          call compute_nonnegative_rho_av(rho(i-1,j,k,:), rho(i,j,k,:), rhoav, dx)
          call compute_sqrtLonsager_local(rhoav,sum(rhoav),sqrtLonsager_x(i,j,k,:))
    end do
    end do
    end do

    ! y-faces
    do k=lo(3),hi(3)
    do j=lo(2),hi(2)+1
    do i=lo(1),hi(1)
          call compute_nonnegative_rho_av(rho(i,j-1,k,:), rho(i,j,k,:), rhoav, dx)
          call compute_sqrtLonsager_local(rhoav,sum(rhoav),sqrtLonsager_y(i,j,k,:))
    end do
    end do
    end do

#if (AMREX_SPACEDIM == 3)
    ! z-faces
    do k=lo(3),hi(3)+1
    do j=lo(2),hi(2)
    do i=lo(1),hi(1)
          call compute_nonnegative_rho_av(rho(i,j,k-1,:), rho(i,j,k,:), rhoav, dx)
          call compute_sqrtLonsager_local(rhoav,sum(rhoav),sqrtLonsager_z(i,j,k,:))
    end do
    end do
    end do
#endif
   
  end subroutine compute_sqrtLonsager_fc


  subroutine compute_nonnegative_rho_av(rho1, rho2, rhoav, dx)
    double precision, intent(in   ) :: rho1(nspecies), rho2(nspecies) ! Densities in two neighboring cells
    double precision, intent(  out) :: rhoav(nspecies)                ! Face-centered average  
    double precision, intent(in   ) :: dx(:)

    double precision :: dv, value1, value2, tmp1, tmp2
    integer :: comp

    ! cell volume
#if (AMREX_SPACEDIM == 2)
    dv = product(dx(1:2))*cell_depth
#elif (AMREX_SPACEDIM == 3)
    dv = product(dx(1:3))
#endif

    ! special version for rtil
    if (use_multiphase .eq. 1 .and. nspecies .eq. 2)then
       value1 = rho1(1)
       value2 = rho2(1)
       if(value1.le.0.or.value2.le.0)then
          rhoav(1) = 0.d0
       else
          rhoav(1) =min( 0.5d0*(value1+value2), rho0)
       endif
       rhoav(2) = rho0 - rhoav(1)

    else

       do comp=1,nspecies
          value1 = rho1(comp)/molmass(comp) ! Convert to number density
          value2 = rho2(comp)/molmass(comp)

          select case(avg_type)
          case(1) ! Arithmetic with a C0-smoothed Heaviside
             if ( (value1 .le. 0.d0) .or. (value2 .le. 0.d0) ) then
                rhoav(comp)=0.d0
             else
                tmp1=min(dv*value1,1.d0)
                tmp2=min(dv*value2,1.d0)
                rhoav(comp)=molmass(comp)*(value1+value2)/2.d0*tmp1*tmp2
             end if
          case(2) ! Geometric
             rhoav(comp)=molmass(comp)*sqrt(max(value1,0.d0)*max(value2,0.d0))
          case(3) ! Harmonic
             ! What we want here is the harmonic mean of max(value1,0) and max(value2,0)
             ! Where we define the result to be zero if either one is zero
             ! But numerically we want to avoid here division by zero
             if ( (value1 .le. 10.d0*tiny(1.d0)) .or. (value2 .le. 10.d0*tiny(1.d0)) ) then
                rhoav(comp)=0.d0
             else
                rhoav(comp)=molmass(comp)*2.d0/(1.d0/value1+1.d0/value2)
             end if
          case(10) ! Arithmetic with (discontinuous) Heaviside
             if ( (value1 .le. 0.d0) .or. (value2 .le. 0.d0) ) then
                rhoav(comp)=0.d0
             else
                rhoav(comp)=molmass(comp)*(value1+value2)/2.d0
             end if
          case(11) ! Arithmetic with C1-smoothed Heaviside
             if ( (value1 .le. 0.d0) .or. (value2 .le. 0.d0) ) then
                rhoav(comp)=0.d0
             else
                tmp1=dv*value1
                if (tmp1<1.d0) then
                   tmp1=(3.d0-2.d0*tmp1)*tmp1**2
                else
                   tmp1=1.d0
                end if
                tmp2=dv*value2
                if (tmp2<1.d0) then
                   tmp2=(3.d0-2.d0*tmp2)*tmp2**2
                else
                   tmp2=1.d0
                end if
                rhoav(comp)=molmass(comp)*(value1+value2)/2.d0*tmp1*tmp2
             endif
          case(12) ! Arithmetic with C2-smoothed Heaviside
             if ( (value1 .le. 0.d0) .or. (value2 .le. 0.d0) ) then
                rhoav(comp)=0.d0
             else
                tmp1=dv*value1
                if (tmp1<1.d0) then
                   tmp1=(10.d0-15.d0*tmp1+6.d0*tmp1**2)*tmp1**3
                else
                   tmp1=1.d0
                end if
                tmp2=dv*value2
                if (tmp2<1.d0) then
                   tmp2=(10.d0-15.d0*tmp2+6.d0*tmp2**2)*tmp2**3
                else
                   tmp2=1.d0
                end if
                rhoav(comp)=molmass(comp)*(value1+value2)/2.d0*tmp1*tmp2
             endif
          case default
             call bl_error("compute_nonnegative_rho_av: invalid avg_type")
          end select
       end do

    end if
       
  end subroutine compute_nonnegative_rho_av

  ! This routine must be called with non-negative densities, i.e., after calling compute_nonnegative_rho_av
  subroutine compute_sqrtLonsager_local(rho,rhotot,sqrtLonsager)
   
    double precision, intent(in)   :: rho(nspecies)            
    double precision, intent(in)   :: rhotot
    double precision, intent(out)  :: sqrtLonsager(nspecies,nspecies) 

    ! local variables
    integer         :: row,column,info
    double precision :: W(nspecies)
    double precision :: rcond

    double precision :: molarconc(nspecies)
    double precision :: molmtot
    double precision :: chi(nspecies,nspecies)
    double precision :: D_bar(nspecies,nspecies)

    integer :: ntrace                   ! number of trace species with w_k < fractional_tolerance
    integer :: nspecies_sub             ! dim of subsystem = nspecies-ntrace 
    double precision :: molmtot_sub
    double precision :: rhotot_sub

    ! this is a mapping used to eliminate elements in D_bar we don't need (for D_bar_sub)
    ! and for expanding sqrtLonsager_sub into sqrtLonsager
    ! it will contain the numbers 1, 2, ..., (nspecies-ntrace)
    ! with zeros in elements corresponding to trace elements
    ! (example) for a 5-species system having trace species 2 and 5:
    !  species       1 2 3 4 5
    !  dest(species) 1 0 2 3 0
    integer :: dest(nspecies)
  
    ! compute the number of trace species
    ! build the mapping for expanding/contracting arrays
    ntrace = 0
    do row=1, nspecies
       W(row) = rho(row)/rhotot
       if (W(row) .lt. fraction_tolerance) then
          ntrace = ntrace + 1
          dest(row) = 0
       else
          dest(row) = row - ntrace
       end if
    end do

    if (ntrace .eq. nspecies-1) then

       ! this is all trace species except for 1 (essentially pure solvent);
       ! set sqrtLonsager to zero
       sqrtLonsager(:,:) = 0.d0

    else if (ntrace .eq. 0) then

       ! there are no trace species

       ! compute molarconc and molmtot
       call compute_molconc_molmtot_local(nspecies,molmass,rho,rhotot,molarconc,molmtot)

       ! compute D_bar
       call compute_D_bar_local(rho,rhotot,D_bar)

       ! compute chi
       call compute_chi(nspecies,molmass,rho,rhotot,molarconc,chi,D_bar)

       ! compute Onsager matrix L (store in sqrtLonsager)
       do column=1, nspecies
          do row=1, nspecies
             sqrtLonsager(row, column) = molmtot*rhotot*W(row)*chi(row,column)*W(column)/k_B
          end do
       end do

       ! compute cell-centered Cholesky factor, sqrtLonsager
       call choldc(sqrtLonsager,nspecies)   

    else

       ! if there are trace species, we consider a subsystem 
       ! consisting of non-trace species

       nspecies_sub = nspecies - ntrace
       call compute_sqrtLonsager_sub()

    end if

  contains
  
    subroutine compute_sqrtLonsager_sub() ! We make this a subroutine to use stack instead of heap
    
       double precision :: molmass_sub(nspecies_sub)
       double precision :: rho_sub(nspecies_sub)
       double precision :: W_sub(nspecies_sub)
       double precision :: molarconc_sub(nspecies_sub)

       double precision :: D_bar_sub(nspecies_sub,nspecies_sub)
       double precision :: chi_sub(nspecies_sub,nspecies_sub)
       double precision :: sqrtLonsager_sub(nspecies_sub,nspecies_sub)

       ! create a vector of non-trace densities and molmass for the subsystem
       do row=1, nspecies
          if (dest(row) .ne. 0) then
             molmass_sub(dest(row)) = molmass(row)
             rho_sub(dest(row)) = rho(row)
          end if
       end do

       ! renormalize total density and mass fractions
       rhotot_sub = sum(rho_sub)
       do row=1, nspecies_sub
          W_sub(row) = rho_sub(row)/rhotot_sub
       end do

       ! first, compute the full D_bar
       ! then, construct D_bar_sub by mapping the full D_bar into D_bar_sub
       ! you could read in only the lower diagonals, 
       ! reflect, and set the diagnals to zero if you want

       call compute_D_bar_local(rho,rhotot,D_bar)

       do row=1, nspecies
          if (dest(row) .eq. 0) then
             cycle
          end if
          do column=1, nspecies
             if (dest(column) .ne. 0) then
                D_bar_sub(dest(row),dest(column)) = D_bar(row,column)
             end if
          end do
       end do
       
       ! compute molarconc_sub and molmtot_sub
       call compute_molconc_molmtot_local(nspecies_sub,molmass_sub,rho_sub,rhotot_sub,molarconc_sub,molmtot_sub)

       ! compute chi_sub
       call compute_chi(nspecies_sub,molmass_sub,rho_sub,rhotot_sub,molarconc_sub,chi_sub,D_bar_sub)

       ! compute Onsager matrix L_sub (store in sqrtLonsager_sub)
       do column=1, nspecies_sub
          do row=1, nspecies_sub
             sqrtLonsager_sub(row, column) = &
                  molmtot_sub*rhotot_sub*W_sub(row)*chi_sub(row,column)*W_sub(column)/k_B
          end do
       end do

       ! compute cell-centered Cholesky factor, sqrtLonsager_sub
       call choldc(sqrtLonsager_sub,nspecies_sub)   

       ! expand sqrtLonsager_sub into sqrtLonsager
       sqrtLonsager(:,:) = 0.d0
       do row=1, nspecies
          if (dest(row) .eq. 0) then
             cycle
          end if
          do column=1, nspecies
             if (dest(column) .ne. 0) then
                sqrtLonsager(row,column) = sqrtLonsager_sub(dest(row),dest(column))
             end if
          end do
       end do
       
    end subroutine compute_sqrtLonsager_sub

    subroutine chol_lapack(sqrtL,nspecies_in)
      integer, intent (in)            :: nspecies_in
      double precision, intent (inout) :: sqrtL(nspecies_in,nspecies_in)

      integer :: row, column
       
      call dpotrf_f95(sqrtL,'L', rcond, 'I', info)
      !stop "LAPACK95 dpotrf_f95 disabled"
    
      ! remove all upper-triangular entries and NXN entry that lapack doesn't set to zero 
      do row=1, nspecies_in
        do column=row+1, nspecies_in
          sqrtL(row,column) = 0.d0
        end do
      end do
      sqrtL(nspecies_in,nspecies_in) = 0.d0
    end subroutine chol_lapack 
       
  end subroutine compute_sqrtLonsager_local


!   subroutine compute_rhoWchi_from_chi(mla,rho,chi,rhoWchi)
 
!     type(ml_layout), intent(in   )  :: mla
!     type(multifab) , intent(in   )  :: rho(:)
!     type(multifab) , intent(in   )  :: chi(:) 
!     type(multifab) , intent(inout)  :: rhoWchi(:) 

!     ! local variables
!     integer :: lo(mla%dim), hi(mla%dim)
!     integer :: n,i,dm,nlevs,ng_1,ng_3,ng_4
 
!     ! pointer for rho(nspecies), molarconc(nspecies) 
!     double precision, pointer        :: dp1(:,:,:,:)  ! for rho    
!     double precision, pointer        :: dp3(:,:,:,:)  ! for chi
!     double precision, pointer        :: dp4(:,:,:,:)  ! for rhoWchi

!     type(mfiter) :: mfi
!     type(box) :: tilebox
!     integer :: tlo(mla%dim), thi(mla%dim)

!     type(bl_prof_timer), save :: bpt

!     call build(bpt,"compute_rhoWchi_from_chi")

!     dm = mla%dim        ! dimensionality
!     ng_1 = rho(1)%ng    ! number of ghost cells 
!     ng_3 = chi(1)%ng
!     ng_4 = rhoWchi(1)%ng
!     nlevs = mla%nlevel  ! number of levels 
 
!     !$omp parallel private(n,i,mfi,tilebox,tlo,thi) &
!     !$omp private(dp1,dp3,dp4,lo,hi)

!     ! loop over all boxes 
!     do n=1,nlevs
!        call mfiter_build(mfi, rho(n), tiling=.true.)

!        do while (more_tile(mfi))
!           i = get_fab_index(mfi)

!           tilebox = get_growntilebox(mfi,rhoWchi(n)%ng)
!           tlo = lwb(tilebox)
!           thi = upb(tilebox)

! !       do i=1,nfabs(rho(n))
!           dp1 => dataptr(rho(n), i)
!           dp3 => dataptr(chi(n), i)
!           dp4 => dataptr(rhoWchi(n), i)
!           lo  =  lwb(get_box(rho(n), i))
!           hi  =  upb(get_box(rho(n), i))
          
!           select case(dm)
!           case (2)
!              call compute_rhoWchi_from_chi_2d(dp1(:,:,1,:),dp3(:,:,1,:),dp4(:,:,1,:), &
!                                      ng_1,ng_3,ng_4,lo,hi,tlo,thi) 
!           case (3)
!              call compute_rhoWchi_from_chi_3d(dp1(:,:,:,:),dp3(:,:,:,:),dp4(:,:,:,:), &
!                                      ng_1,ng_3,ng_4,lo,hi,tlo,thi) 
!           end select
!        end do
!     end do
!     !$omp end parallel

!     call destroy(bpt)

!   end subroutine compute_rhoWchi_from_chi
  
!   subroutine compute_rhoWchi_from_chi_2d(rho,chi,rhoWchi,ng_1,ng_3,ng_4,glo,ghi,tlo,thi)
  
!     integer          :: glo(2), ghi(2), ng_1,ng_3,ng_4,tlo(2),thi(2)
!     double precision, intent(in   ) ::     rho(glo(1)-ng_1:,glo(2)-ng_1:,:) ! density; last dimension for species
!     double precision, intent(in   ) ::     chi(glo(1)-ng_3:,glo(2)-ng_3:,:) ! last dimension for nspecies^2
!     double precision, intent(in   ) :: rhoWchi(glo(1)-ng_4:,glo(2)-ng_4:,:) ! last dimension for nspecies^2

!     ! local variables
!     integer          :: i,j
  
!     ! for specific box, now start loops over alloted cells 
!     do j=tlo(2),thi(2)
!        do i=tlo(1),thi(1)
        
!           call compute_rhoWchi_from_chi_local(rho(i,j,:),chi(i,j,:),rhoWchi(i,j,:))

!        end do
!     end do

!   end subroutine compute_rhoWchi_from_chi_2d

!   subroutine compute_rhoWchi_from_chi_3d(rho,chi,rhoWchi,ng_1,ng_3,ng_4,glo,ghi,tlo,thi)

!     integer          :: glo(3), ghi(3), ng_1,ng_3,ng_4,tlo(3),thi(3)
!     double precision, intent(in   ) ::     rho(glo(1)-ng_1:,glo(2)-ng_1:,glo(3)-ng_1:,:) ! density; last dimension for species
!     double precision, intent(in   ) ::     chi(glo(1)-ng_3:,glo(2)-ng_3:,glo(3)-ng_3:,:) ! last dimension for nspecies^2
!     double precision, intent(in   ) :: rhoWchi(glo(1)-ng_4:,glo(2)-ng_4:,glo(3)-ng_4:,:) ! last dimension for nspecies^2
    
!     ! local variables
!     integer          :: i,j,k

!     ! for specific box, now start loops over alloted cells 
!     do k=tlo(3),thi(3)
!        do j=tlo(2),thi(2)
!           do i=tlo(1),thi(1)
       
!              call compute_rhoWchi_from_chi_local(rho(i,j,k,:),chi(i,j,k,:),rhoWchi(i,j,k,:))
              
!          end do
!       end do
!     end do
   
!   end subroutine compute_rhoWchi_from_chi_3d
  
!   subroutine compute_rhoWchi_from_chi_local(rho,chi,rhoWchi)
   
!     double precision, intent(in)   :: rho(nspecies)            
!     double precision, intent(in)   :: chi(nspecies,nspecies)   ! rank conversion done 
!     double precision, intent(out)  :: rhoWchi(nspecies,nspecies) 
 
!     ! local variables
!     integer                       :: row,column

!     ! populate rho*W*chi = rho_i*chi
!     do row=1, nspecies
!        do column=1, nspecies
!           rhoWchi(row,column) = rho(row)*chi(row,column)  
!        end do
!     end do

!   end subroutine compute_rhoWchi_from_chi_local

!   subroutine compute_baro_coef(mla,baro_coef,rho,rhotot,Temp)
 
!     type(ml_layout), intent(in   )  :: mla
!     type(multifab) , intent(inout)  :: baro_coef(:)
!     type(multifab) , intent(in   )  :: rho(:)
!     type(multifab) , intent(in   )  :: rhotot(:) 
!     type(multifab) , intent(in   )  :: Temp(:) 

!     ! local
!     integer :: lo(mla%dim), hi(mla%dim)
!     integer :: n,i,dm,nlevs,ng_1,ng_2,ng_3,ng_4

!     double precision, pointer        :: dp1(:,:,:,:)  ! for baro_coef
!     double precision, pointer        :: dp2(:,:,:,:)  ! for rho
!     double precision, pointer        :: dp3(:,:,:,:)  ! for rhotot
!     double precision, pointer        :: dp4(:,:,:,:)  ! for Temp

!     type(bl_prof_timer), save :: bpt

!     call build(bpt,"compute_baro_coef")

!     dm = mla%dim
!     nlevs = mla%nlevel

!     ng_1 = baro_coef(1)%ng
!     ng_2 = rho(1)%ng
!     ng_3 = rhotot(1)%ng
!     ng_4 = Temp(1)%ng

!     do n=1,nlevs
!        do i=1,nfabs(rho(n))
!           dp1 => dataptr(baro_coef(n),i)
!           dp2 => dataptr(rho(n),i)
!           dp3 => dataptr(rhotot(n),i)
!           dp4 => dataptr(Temp(n),i)
!           lo = lwb(get_box(rho(n),i))
!           hi = upb(get_box(rho(n),i))
!           select case(dm)
!           case (2)
!              call compute_baro_coef_2d(dp1(:,:,1,:),dp2(:,:,1,:),dp3(:,:,1,1),dp4(:,:,1,1), &
!                                         ng_1,ng_2,ng_3,ng_4,lo,hi) 
!           case (3)
!              call bl_error("compute_baro_coef_3d not written yet")
!           end select
!        end do
!     end do

!     call destroy(bpt)

!   end subroutine compute_baro_coef

!   subroutine compute_baro_coef_2d(baro_coef,rho,rhotot,Temp,ng_1,ng_2,ng_3,ng_4,lo,hi)
 
!     integer        , intent(in   ) :: lo(2), hi(2), ng_1, ng_2, ng_3, ng_4
!     double precision, intent(inout) :: baro_coef(lo(1)-ng_1:,lo(2)-ng_1:,:)
!     double precision, intent(in   ) ::       rho(lo(1)-ng_2:,lo(2)-ng_2:,:)
!     double precision, intent(in   ) ::    rhotot(lo(1)-ng_3:,lo(2)-ng_3:)
!     double precision, intent(in   ) ::      Temp(lo(1)-ng_4:,lo(2)-ng_4:)
    
!     ! local variables
!     integer :: i,j,comp
!     double precision :: n

!     do j=lo(2)-ng_1,hi(2)+ng_1
!        do i=lo(1)-ng_1,hi(1)+ng_1

!           n = 0.d0
!           do comp=1,nspecies
!              n = n + rho(i,j,comp)/molmass(comp)
!           end do

!           do comp=1,nspecies
!              baro_coef(i,j,comp) = rho(i,j,comp)/rhobar(comp) / (n*k_B*Temp(i,j))
!           end do

!        end do
!     end do

!   end subroutine compute_baro_coef_2d
  
end module mass_flux_utilities_module
