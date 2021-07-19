module matrix_utilities_module

  use amrex_error_module
  use common_namelist_module
  use multispec_namelist_module

  implicit none

  private

  public :: Dbar2chi_iterative, choldc
  
contains

    ! nspecies_local is number of species
    ! num_iterations is the number of terms in the sum to use: 3-5 are reasonable values
    ! D_bar is matrix of Maxwell-Stefan binary diffusion coefficient
    ! chi is the multispecies diffusion matrix
    ! Xk is mole fractions --- MUST NOT BE ZERO
    subroutine Dbar2chi_iterative(nspecies_local,num_iterations,D_bar,Xk,molmass_local,chi) bind(C,name="Dbar2chi_iterative")
      integer, intent(in) :: nspecies_local
      integer, intent(in) :: num_iterations
      double precision, intent(in) :: D_bar(1:nspecies_local,1:nspecies_local)
      double precision, intent(in) :: Xk(1:nspecies_local), molmass_local(1:nspecies_local)
      double precision, intent(out) :: chi(1:nspecies_local,1:nspecies_local)
      
      ! Local variables
      double precision :: term1, term2, MWmix
      double precision :: Di(1:nspecies_local)
      double precision :: Deltamat(1:nspecies_local,1:nspecies_local), Zmat(1:nspecies_local,1:nspecies_local)
      double precision, dimension(1:nspecies_local,1:nspecies_local) :: Pmat, Jmat
      double precision, dimension(1:nspecies_local) :: Minv, Mmat
      double precision, dimension(1:nspecies_local,1:nspecies_local) :: PJ, matrix1, matrix2
      double precision :: scr
      double precision :: Ykp(1:nspecies_local), Xkp(1:nspecies_local)

      integer :: i, j, k, ii, jj
      
      ! mole fractions correction
      ! Turned this off since it should be done in the caller
      do ii = 1, nspecies_local
       Xkp(ii) = Xk(ii)
      end do

      ! molecular weight of mixture - EGLIB
      Mwmix = 0.0d0
      do ii = 1, nspecies_local
       MWmix = MWmix + Xkp(ii)*molmass_local(ii)
      end do

      ! mass fractions correction - EGLIB
      do ii = 1, nspecies_local
       Ykp(ii) = molmass_local(ii)/MWmix*Xkp(ii)
      end do

      ! Find Di matrix 
      do i = 1, nspecies_local
       term2 = 0.0d0
       do j = 1, nspecies_local
        if(j.ne.i) then
          term2 = term2 + Xkp(j)/D_bar(i,j)
        end if
       end do   
       Di(i) = (1.d0-Ykp(i))/term2 
      end do   

      ! Compute Mmat and Minv
      do i = 1, nspecies_local
       Mmat(i) = Xkp(i)/Di(i)
       Minv(i) = Di(i)/Xkp(i)
      end do
      

      ! Compute P matrix
      Pmat = 0.0d0
      do i = 1, nspecies_local
       do j = 1, nspecies_local
         Pmat(i,j) = - Ykp(j) 
         if(i.eq.j) then
          Pmat(i,j) =  Pmat(i,j) + 1.0d0  
         end if
       end do
      end do

      ! Compute Deltamat
      Deltamat = 0.0d0 
      do i = 1, nspecies_local
       do j = 1, nspecies_local
         if(i.eq.j) then
          term1 = 0.0d0
          do k = 1, nspecies_local
           if(k.ne.i) then
            term1 = term1 + Xkp(i)*Xkp(k)/D_bar(i,k)
           end if
          end do  
          Deltamat(i,i) = term1
         else
          Deltamat(i,j) = -Xkp(i)*Xkp(j)/D_bar(i,j) 
         end if  
          Zmat(i,j) = -Deltamat(i,j)
       end do
      end do  

      ! Compute Zmat
      do i = 1, nspecies_local
        Zmat(i,i) = Zmat(i,i) + Mmat(i)
      end do  

      ! Compute Jmat
      do i = 1, nspecies_local
       do j = 1, nspecies_local
         Jmat(i,j) = Minv(i)*Zmat(i,j)
        end do
       end do

      ! Compute PJ
      PJ = 0.0d0
      do i = 1, nspecies_local
       do j = 1, nspecies_local
        do k = 1, nspecies_local
         PJ(i,j) = PJ(i,j) + Pmat(i,k)*Jmat(k,j)
        end do
       end do
      end do

      ! Compute P M^-1 Pt; store it in matrix2
      do i = 1, nspecies_local
       do j = 1, nspecies_local
        scr = 0.d0
        do k = 1, nspecies_local
         scr = scr + Pmat(i,k)*Minv(k)*Pmat(j,k) 
            ! notice the change in indices for Pmat to represent Pmat^t
        end do
         matrix2(i,j) = scr
         chi(i,j) = scr
       end do
      end do


      do jj = 1,num_iterations
       do i = 1, nspecies_local
        do j = 1, nspecies_local
         scr = 0.d0
         do k = 1, nspecies_local
            scr = scr + PJ(i,k)*chi(k,j)
         end do
          matrix1(i,j) = scr+matrix2(i,j)
        end do
       end do 
       chi=matrix1
      end do

write(*,*) "HACK: chi row 1 " , chi(1,1) , " " , chi(1,2) , " " , chi(1,3) 
write(*,*) "HACK: chi row 2 " , chi(2,1) , " " , chi(2,2) , " " , chi(2,3) 
write(*,*) "HACK: chi row 3 " , chi(3,1) , " " , chi(3,2) , " " , chi(3,3) 

  end subroutine

   ! a is input matrix.  
   ! upon return the lower triangle and diagonal are overwritten by the cholesky factor
   subroutine choldc(a,np)
       integer :: np
       double precision, intent(inout) :: a(np,np)

       double precision :: p(np), dij(np,np)
       double precision :: dd(np,np)
       double precision :: yy(np), mwmix 

       integer :: i, j, k, ii, jj
       double precision :: sum1
       double precision :: small_number = 0.0d0 ! Some tolerance

       integer :: idiag,ising

       do i = 1, np

           ising = 0

        do j = i, np

           sum1 = a(i,j)

           do k = i-1, 1, -1

              sum1 = sum1 - a(i,k)*a(j,k)

           end do

           if(i.eq.j) then

             if(sum1.le.small_number) then

             p(i) = 0.d0

             ising = 1

             else

             p(i) = sqrt(sum1)

             end if

           else

             if(ising.eq.0)then

                a(j,i) = sum1/p(i)

             else

                a(j,i) = 0.d0

             end if

           end if

        end do

       end do


       do i = 1, np

          do j = i+1, np

           a(i,j) = 0.0d0 ! Zero upper triangle

          end do
          
          a(i,i) = p(i)

       end do

    end subroutine

end module matrix_utilities_module
