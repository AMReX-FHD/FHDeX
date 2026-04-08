subroutine advect_2d(time, lo, hi, &
     &            uin_p , uip_lo, uip_hi, &
     &            uin_f , uif_lo, uif_hi, &
     &            uout, uo_lo, uo_hi, &
     &            ptSp, ptsp_lo, ptsp_hi,&
     &            ptSf, ptsf_lo, ptsf_hi,&
     &            ifacep, ifp_lo, ifp_hi,&
     &            ifacef, iff_lo, iff_hi,&
     &            vx_p  , vxp_lo, vxp_hi, &
     &            vy_p  , vyp_lo, vyp_hi, &
     &            vx_f  , vxf_lo, vxf_hi, &
     &            vy_f  , vyf_lo, vyf_hi, &
     &            flxx, fx_lo, fx_hi, &
     &            flxy, fy_lo, fy_hi, &
     &            flxx1, fx1_lo, fx1_hi, &
     &            flxy1, fy1_lo, fy1_hi, &
     &            flxx2, fx2_lo, fx2_hi, &
     &            flxy2, fy2_lo, fy2_hi, &
     &            dx,dt,nu, correct) bind(C, name="advect_2d")

  use amrex_mempool_module, only : bl_allocate, bl_deallocate
  use compute_flux_module_2d, only : compute_flux_2d

  implicit none
  ! correct = 0 if we are making a predictor concentration and correct =1 if we are solving for a corrected concentration
  integer, intent(in) :: correct
  ! dx - grid size
  ! dt - time step
  ! parameters nu is the diffusion coeffiecent
  double precision, intent(in) :: dx(2), dt, time,nu
  ! work region

  integer, intent(in) :: lo(2), hi(2)
  integer, intent(in) :: uip_lo(2), uip_hi(2)
  integer, intent(in) :: uif_lo(2), uif_hi(2)
  integer, intent(in) :: uo_lo(2), uo_hi(2)
  integer, intent(in) :: ptsp_lo(2), ptsp_hi(2)
  integer, intent(in) :: ptsf_lo(2), ptsf_hi(2)
  integer, intent(in) :: ifp_lo(2), ifp_hi(2)
  integer, intent(in) :: iff_lo(2), iff_hi(2)
  integer, intent(in) :: vxf_lo(2), vxf_hi(2)
  integer, intent(in) :: vyf_lo(2), vyf_hi(2)
  integer, intent(in) :: vxp_lo(2), vxp_hi(2)
  integer, intent(in) :: vyp_lo(2), vyp_hi(2)
  integer, intent(in) :: fx_lo(2), fx_hi(2)
  integer, intent(in) :: fy_lo(2), fy_hi(2)
  integer, intent(in) :: fx1_lo(2), fx1_hi(2)
  integer, intent(in) :: fy1_lo(2), fy1_hi(2)
  integer, intent(in) :: fx2_lo(2), fx2_hi(2)
  integer, intent(in) :: fy2_lo(2), fy2_hi(2)
  ! ** IN: uin_p - previous concentration (c^n)
  !        uin_f - predicted concentration (c^*), (denoted f since it is only used in determining the flux)
  !        ptSp/ ptsf     - previous/ predicted point sources
  !        ifacep/ ifacef - previous/ predicted locating of interface
  !        vx_p/ vx_f     - previous/ predicted x component of velocity
  !        vy_p/ vy_f     - previous/ predicted y component of velocity
  !  NOTE: if correct=0 the previous/ predicted quantities are the same

  double precision, intent(in   ) :: uin_p (uip_lo(1):uip_hi(1),uip_lo(2):uip_hi(2))
  double precision, intent(in   ) :: uin_f (uif_lo(1):uif_hi(1),uif_lo(2):uif_hi(2))
  double precision, intent(in) :: ptSp(ptsp_lo(1):ptsp_hi(1),ptsp_lo(2):ptsp_hi(2))
  double precision, intent(in) :: ptSf(ptsf_lo(1):ptsf_hi(1),ptsf_lo(2):ptsf_hi(2))
  integer, intent(in) :: ifacep(ifp_lo(1):ifp_hi(1),ifp_lo(2):ifp_hi(2))
  integer, intent(in) :: ifacef(iff_lo(1):iff_hi(1),iff_lo(2):iff_hi(2))
  double precision, intent(in   ) :: vx_p  (vxp_lo(1):vxp_hi(1),vxp_lo(2):vxp_hi(2))
  double precision, intent(in   ) :: vy_p  (vyp_lo(1):vyp_hi(1),vyp_lo(2):vyp_hi(2))
  double precision, intent(in   ) :: vx_f  (vxf_lo(1):vxf_hi(1),vxf_lo(2):vxf_hi(2))
  double precision, intent(in   ) :: vy_f  (vyf_lo(1):vyf_hi(1),vyf_lo(2):vyf_hi(2))
  ! ** Out: uout - predicted concentration c^*( if correct = 0 ) or corrected concentration c^n+1 (if correct = 1)
  !         flxx1, flxy1  - the x and y componetents of the previous concentration fluxes
  !         flxx2, flxy2  - the x and y componetents of the predicted concentration fluxes
  !         flxx, flxy    - the x and y componetents of the average of the previous and predicted concentration fluxes
  double precision, intent(inout) :: uout(uo_lo(1):uo_hi(1),uo_lo(2):uo_hi(2))
  double precision, intent(  out) :: flxx(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2))
  double precision, intent(  out) :: flxy(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2))

  double precision, intent(  out) :: flxx1(fx1_lo(1):fx1_hi(1),fx1_lo(2):fx1_hi(2))
  double precision, intent(  out) :: flxy1(fy1_lo(1):fy1_hi(1),fy1_lo(2):fy1_hi(2))

  double precision, intent(  out) :: flxx2(fx2_lo(1):fx2_hi(1),fx2_lo(2):fx2_hi(2))
  double precision, intent(  out) :: flxy2(fy2_lo(1):fy2_hi(1),fy2_lo(2):fy2_hi(2))

  integer :: i, j
  integer :: glo(2), ghi(2)
  double precision :: dtdx(2), umax, vmax, conmax_in, conmax_out
  double precision :: flxx11, flxy11, flxx22, flxy22
  ! Some compiler may not support 'contiguous'.  Remove it in that case.
  double precision, dimension(:,:), pointer, contiguous :: &
       conx1, conx1_y, cony1, cony1_x, slope1
  double precision, dimension(:,:), pointer, contiguous :: &
       conx2, conx2_y, cony2, cony2_x, slope2

  dtdx = dt/dx

  glo = lo - 1
  ghi = hi + 1


  ! edge states
  call bl_allocate(conx1  ,glo(1), ghi(1), glo(2), ghi(2))
  call bl_allocate(conx1_y,glo(1), ghi(1), glo(2), ghi(2))
  call bl_allocate(cony1  ,glo(1), ghi(1), glo(2), ghi(2))
  call bl_allocate(cony1_x,glo(1), ghi(1), glo(2), ghi(2))
  ! slope
  call bl_allocate(slope1,glo(1), ghi(1), glo(2), ghi(2) )


  call bl_allocate(conx2  ,glo(1), ghi(1), glo(2), ghi(2))
  call bl_allocate(conx2_y,glo(1), ghi(1), glo(2), ghi(2))
  call bl_allocate(cony2  ,glo(1), ghi(1), glo(2), ghi(2))
  call bl_allocate(cony2_x,glo(1), ghi(1), glo(2), ghi(2))
  ! slope
  call bl_allocate(slope2,glo(1), ghi(1), glo(2), ghi(2))

  ! We like to allocate these **pointers** here and then pass them to a function
  ! to remove their pointerness for performance, because normally pointers could
  ! be aliasing.  We need to use pointers instead of allocatable arrays because
  ! we like to use AMReX's bl_allocate to allocate memeory instead of the intrinsic
  ! allocate.  Bl_allocate is much faster than allocate inside OMP.
  ! Note that one MUST CALL BL_DEALLOCATE.

  ! check if CFL condition is violated.

 umax = maxval(abs(vx_f))
  vmax = maxval(abs(vy_f))

  conmax_in=maxval(abs(uin_p))
  if ( umax*dt .ge. dx(1) .or. &
       vmax*dt .ge. dx(2) ) then
!!! TEMPORARILY UNTIL couple fluid velocity
     call bl_error("CFL violation. Use smaller adv.cfl.")
  end if

  ! call a function to compute flux
  call compute_flux_2d(lo, hi, dt, dx, &
                       uin_p, uip_lo, uip_hi, &
                       vx_p, vxp_lo, vxp_hi, &
                       vy_p, vyp_lo, vyp_hi, &
                       flxx1, fx1_lo, fx1_hi, &
                       flxy1, fy1_lo, fy1_hi, &
                       conx1, conx1_y, &
                       cony1, cony1_x, &
                       slope1, glo, ghi)

  call compute_flux_2d(lo, hi, dt, dx, &
                       uin_f, uif_lo, uif_hi, &
                       vx_f, vxf_lo, vxf_hi, &
                       vy_f, vyf_lo, vyf_hi, &
                       flxx2, fx2_lo, fx2_hi, &
                       flxy2, fy2_lo, fy2_hi, &
                       conx2, conx2_y,  &
                       cony2, cony2_x,  &
                       slope2, glo, ghi)

  ! Fluxes with out slopes, not sure why both are here
  do    j = lo(2), hi(2)
     do i = lo(1), hi(1)+1
        flxx1(i,j) = conx1(i,j) * vx_p(i,j)-nu*(uin_p(i,j)-uin_p(i-1,j))/dx(1)
        flxx2(i,j) = conx2(i,j) * vx_f(i,j)-nu*(uin_f(i,j)-uin_f(i-1,j))/dx(1)
     end do
  end do
  !
  do    j = lo(2), hi(2)+1
     do i = lo(1), hi(1)
        flxy1(i,j) = cony1(i,j) * vy_p(i,j)-nu*(uin_p(i,j)-uin_p(i,j-1))/dx(2)
        flxy2(i,j) = cony2(i,j) * vy_f(i,j)-nu*(uin_f(i,j)-uin_f(i,j-1))/dx(2)
     end do
  end do

  ! Do a conservative update
     do    j = lo(2), hi(2)
        do i = lo(1), hi(1)

           flxx11=flxx1(i,j)
           flxy11=flxy1(i,j)
           flxx22=flxx2(i,j)
           flxy22=flxy2(i,j)
           ! Enforce a no flux boundary condition at the interface of the particle, ie if we are evaluating a flux inside the sphere set the flux equal to zero

           if (ifacep(i,j) .eq. 2) then
           uout(i,j)= uin_p(i,j)
           else
           if (ifacep(i,j) .eq. 1) then
               if (ifacep(i+1,j) .eq. 2)then
               flxx1(i+1,j)=0
               else if (ifacep(i-1,j) .eq. 2) then
               flxx1(i,j)=0
               end if

               if (ifacep(i,j+1) .eq. 2)then
               flxy1(i,j+1)=0
               else if (ifacep(i,j-1) .eq. 2) then
               flxy1(i,j)=0
               end if
           endif
           if (ifacef(i,j) .eq. 1) then
               if (ifacef(i+1,j) .eq. 2)then
               flxx2(i+1,j)=0
               else if (ifacef(i-1,j) .eq. 2) then
               flxx2(i,j)=0
               end if

               if (ifacef(i,j+1) .eq. 2)then
               flxy2(i,j+1)=0
               else if (ifacef(i,j-1) .eq. 2) then
               flxy2(i,j)=0
               end if
           end if
           uout(i,j) = uin_p(i,j) + &
              0.5*(( (flxx1(i,j) - flxx1(i+1,j)) * dtdx(1) &
                + (flxy1(i,j) - flxy1(i,j+1)) * dtdx(2)) &
                + ( (flxx2(i,j) - flxx2(i+1,j)) * dtdx(1) &
                + (flxy2(i,j) - flxy2(i,j+1)) * dtdx(2))) &
                + dt*0.5*(ptSp(i,j)+ptSf(i,j))
             end if

            flxx(i,j)=0.5*(flxx1(i,j)+flxx2(i,j))
            flxy(i,j)=0.5*(flxy1(i,j)+flxy2(i,j))
        enddo
     enddo
  ! Old code for AMR
  ! Scale by face area in order to correctly reflx
     do    j = lo(2), hi(2)
        do i = lo(1), hi(1)+1
           flxx(i,j) = flxx(i,j) * (dt * dx(2))
        enddo
     enddo
     do    j = lo(2), hi(2)+1
        do i = lo(1), hi(1)
           flxy(i,j) = flxy(i,j) * (dt * dx(1))
        enddo
     enddo

  ! deallocate pointers
  call bl_deallocate(conx1  )
  call bl_deallocate(conx1_y)
  call bl_deallocate(cony1  )
  call bl_deallocate(cony1_x)

  call bl_deallocate(slope1)

  call bl_deallocate(conx2  )
  call bl_deallocate(conx2_y)
  call bl_deallocate(cony2  )
  call bl_deallocate(cony2_x)
  call bl_deallocate(slope2)

end subroutine advect_2d
