module trans_module

  use amrex_fort_module, only : amrex_real
  use common_namelist_module, only : ngc, nvars, nprimvars, diameter, max_species, molmass, k_b, nspecies

  implicit none

  private

  public :: trans_coeffs, get_trans_coeff

contains

  subroutine trans_coeffs(lo,hi, prim, eta, zeta, kappa) bind(C,name="trans_coeffs")

      integer         , intent(in   ) :: lo(3),hi(3)

      real(amrex_real), intent(in   ) :: prim(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3), nprimvars)

      real(amrex_real), intent(inout) :: eta(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3))
      real(amrex_real), intent(inout) :: zeta(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3))
      real(amrex_real), intent(inout) :: kappa(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3))

      ! real(amrex_real), intent(inout) :: chi(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3),nspecies)
      ! real(amrex_real), intent(inout) :: Dij(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3),nspecies,nspecies)

      integer :: i,j,k,l
      real(amrex_real) :: mgrams(MAX_SPECIES)
      real(amrex_real) :: vconst(nspecies), tconst(nspecies), R(nspecies), gamma1, gamma2, rootT, specaveta, specavkappa

      mgrams = molmass/(6.022140857e23)

      R = k_B/mgrams

      gamma1 = 1.2700
      gamma2 = 1.9223

      do l = 1, nspecies

        vconst(l) = sqrt(R(l))*(mgrams(l)/(diameter(l)**2))*(gamma1/4)*sqrt(1/3.14159265359)
        tconst(l) = R(l)*sqrt(R(l))*(mgrams(l)/(diameter(l)**2))*(5*gamma2/8)*sqrt(1/3.14159265359)

      enddo

      !Hard sphere for now. Check these defs
      !Per Sone (Sone Grad-Hilbert, Bird Chapman Enskog? Check this)

      ! print *, "Hack: dims = ", lo(1),hi(1), lo(2),hi(2), lo(3),hi(3)

      do k = lo(3),hi(3)
        do j = lo(2),hi(2)
          do i = lo(1),hi(1)

            specaveta = 0
            specavkappa = 0

            do l = 1, nspecies
              specaveta = specaveta + prim(i,j,k,6+l)*vconst(l)
              specavkappa = specavkappa + prim(i,j,k,6+l)*tconst(l)
            enddo

            rootT = sqrt(prim(i,j,k,5))

            zeta(i,j,k) = 0 !no bulk viscosity for now

            eta(i,j,k) = rootT*specaveta
            kappa(i,j,k) = rootT*specavkappa
            
            ! if((i.eq.0).and.(j.eq.0).and.(k.eq.0)) then
            !    print *, "transcoef: Hack = ", i, j, k, specaveta, specavkappa, rootT, prim(i,j,k,5)
            ! endif

            if(kappa(i,j,k) .ne. kappa(i,j,k)) then
               print *, "NAN! kappa ", i, j, k, prim(i,j,k,5)
               call exit()
            endif

            ! print*, "Hack (trans_coeffs): eta = ", eta(i,j,k)
            ! stop

          enddo
        enddo
      enddo

  end subroutine trans_coeffs

  subroutine get_trans_coeff(temp, fracvec, eta, zeta, kappa) bind(C,name="get_trans_coeff")


      double precision, intent(in   ) :: temp, fracvec(nspecies)
      double precision, intent(inout) :: eta, zeta, kappa

      integer :: l 
      real(amrex_real) :: mgrams(MAX_SPECIES)
      real(amrex_real) :: vconst(nspecies), tconst(nspecies), R(nspecies), gamma1, gamma2, rootT, specaveta, specavkappa

      mgrams = molmass/(6.022140857e23)
      R = k_B/mgrams

      gamma1 = 1.2700
      gamma2 = 1.9223

      do l = 1, nspecies

        vconst(l) = sqrt(R(l))*(mgrams(l)/(diameter(l)**2))*(gamma1/4)*sqrt(1/3.14159265359)
        tconst(l) = R(l)*sqrt(R(l))*(mgrams(l)/(diameter(l)**2))*(5*gamma2/8)*sqrt(1/3.14159265359)

      enddo

      specaveta = 0
      specavkappa = 0

      do l = 1, nspecies
        specaveta = specaveta + fracvec(l)*vconst(l)
        specavkappa = specavkappa + fracvec(l)*tconst(l)
      enddo

      rootT = sqrt(temp)

      zeta = 0 !no bulk viscosity for now

      eta = rootT*specaveta
      kappa = rootT*specavkappa

      if(kappa .ne. kappa) then
         print *, "NAN! kappa bc", temp, specavkappa
         call exit()
      endif


  end subroutine get_trans_coeff

end module trans_module






