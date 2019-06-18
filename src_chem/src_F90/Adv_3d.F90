subroutine advect_3d(time, lo, hi, &
     &            uin , ui_lo, ui_hi, &
     &            uout, uo_lo, uo_hi, &
     &            ptS, pts_lo, pts_hi,&
     &            iface, if_lo, if_hi,&
     &            vx  , vx_lo, vx_hi, &
     &            vy  , vy_lo, vy_hi, &
     &            vz  , vz_lo, vz_hi, &
     &            flxx, fx_lo, fx_hi, &
     &            flxy, fy_lo, fy_hi, &
     &            flxz, fz_lo, fz_hi, &
     &            dx,dt,nu) bind(C, name="advect_3d")
  
  use amrex_mempool_module, only : bl_allocate, bl_deallocate
  use compute_flux_module_3d, only : compute_flux_3d

  implicit none

  integer, intent(in) :: lo(3), hi(3)
  double precision, intent(in) :: dx(3), dt, time,nu
  integer, intent(in) :: ui_lo(3), ui_hi(3)
  integer, intent(in) :: uo_lo(3), uo_hi(3)
  integer, intent(in) :: pts_lo(3), pts_hi(3)
  integer, intent(in) :: if_lo(3), if_hi(3)
  integer, intent(in) :: vx_lo(3), vx_hi(3)
  integer, intent(in) :: vy_lo(3), vy_hi(3)
  integer, intent(in) :: vz_lo(3), vz_hi(3)
  integer, intent(in) :: fx_lo(3), fx_hi(3)
  integer, intent(in) :: fy_lo(3), fy_hi(3)
  integer, intent(in) :: fz_lo(3), fz_hi(3)
  double precision, intent(in   ) :: uin (ui_lo(1):ui_hi(1),ui_lo(2):ui_hi(2),ui_lo(3):ui_hi(3))
  double precision, intent(inout) :: uout(uo_lo(1):uo_hi(1),uo_lo(2):uo_hi(2),uo_lo(3):uo_hi(3))
  double precision, intent(in) :: ptS(pts_lo(1):pts_hi(1),pts_lo(2):pts_hi(2),pts_lo(3):pts_hi(3))
  integer, intent(in) :: iface(if_lo(1):if_hi(1),if_lo(2):if_hi(2),if_lo(3):if_hi(3))
  double precision, intent(in   ) :: vx  (vx_lo(1):vx_hi(1),vx_lo(2):vx_hi(2),vx_lo(3):vx_hi(3))
  double precision, intent(in   ) :: vy  (vy_lo(1):vy_hi(1),vy_lo(2):vy_hi(2),vy_lo(3):vy_hi(3))
  double precision, intent(in   ) :: vz  (vz_lo(1):vz_hi(1),vz_lo(2):vz_hi(2),vz_lo(3):vz_hi(3))
  double precision, intent(  out) :: flxx(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2),fx_lo(3):fx_hi(3))
  double precision, intent(  out) :: flxy(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2),fy_lo(3):fy_hi(3))
  double precision, intent(  out) :: flxz(fz_lo(1):fz_hi(1),fz_lo(2):fz_hi(2),fz_lo(3):fz_hi(3))

  integer :: i, j, k
  integer :: glo(3), ghi(3)
  double precision :: dtdx(3), umax, vmax, wmax, conmax_in, conmax_out

  ! Some compiler may not support 'contiguous'.  Remove it in that case.
  double precision, dimension(:,:,:), pointer, contiguous :: &
       conx, conx_y, conx_z, cony, cony_x, cony_z, conz, conz_x, conz_y, slope

  dtdx = dt/dx

  glo = lo - 1
  ghi = hi + 1

  ! edge states
  call bl_allocate(conx  ,glo(1), ghi(1), glo(2), ghi(2), glo(3), ghi(3))
  call bl_allocate(conx_y,glo(1), ghi(1), glo(2), ghi(2), glo(3), ghi(3))
  call bl_allocate(conx_z,glo(1), ghi(1), glo(2), ghi(2), glo(3), ghi(3))
  call bl_allocate(cony  ,glo(1), ghi(1), glo(2), ghi(2), glo(3), ghi(3))
  call bl_allocate(cony_x,glo(1), ghi(1), glo(2), ghi(2), glo(3), ghi(3))
  call bl_allocate(cony_z,glo(1), ghi(1), glo(2), ghi(2), glo(3), ghi(3))
  call bl_allocate(conz  ,glo(1), ghi(1), glo(2), ghi(2), glo(3), ghi(3))
  call bl_allocate(conz_x,glo(1), ghi(1), glo(2), ghi(2), glo(3), ghi(3))
  call bl_allocate(conz_y,glo(1), ghi(1), glo(2), ghi(2), glo(3), ghi(3))
  ! slope
  call bl_allocate(slope,glo(1), ghi(1), glo(2), ghi(2), glo(3), ghi(3))  
  
  ! We like to allocate these **pointers** here and then pass them to a function
  ! to remove their pointerness for performance, because normally pointers could
  ! be aliasing.  We need to use pointers instead of allocatable arrays because
  ! we like to use AMReX's bl_allocate to allocate memeory instead of the intrinsic
  ! allocate.  Bl_allocate is much faster than allocate inside OMP.  
  ! Note that one MUST CALL BL_DEALLOCATE.

  ! check if CFL condition is violated.
  umax = maxval(abs(vx))
  vmax = maxval(abs(vy))
  wmax = maxval(abs(vz))

  conmax_in=maxval(abs(uin))
  if ( umax*dt .ge. dx(1) .or. &
       vmax*dt .ge. dx(2) .or. &
       wmax*dt .ge. dx(3) ) then
     print *, "umax = ", umax, ", vmax = ", vmax, ", wmax = ", wmax, ", dt = ", dt, ", dx = ", dx
!!! TEMPORARILY UNTIL couple fluid velocity
     call bl_error("CFL violation. Use smaller adv.cfl.")
  end if

  ! call a function to compute flux
  call compute_flux_3d(lo, hi, dt, dx, &
                       uin, ui_lo, ui_hi, &
                       vx, vx_lo, vx_hi, &
                       vy, vy_lo, vy_hi, &
                       vz, vz_lo, vz_hi, &
                       flxx, fx_lo, fx_hi, &
                       flxy, fy_lo, fy_hi, &
                       flxz, fz_lo, fz_hi, &
                       conx, conx_y, conx_z, &
                       cony, cony_x, cony_z, &
                       conz, conz_x, conz_y, &
                       slope, glo, ghi,nu)

  ! Do a conservative update
  do       k = lo(3), hi(3)
     do    j = lo(2), hi(2)
        do i = lo(1), hi(1)
           if (iface(i,j,k) .eq. 2) then
           uout(i,j,k)= uin(i,j,k)
           else
           if (iface(i,j,k) .eq. 1) then
               if (iface(i+1,j,k) .eq. 2)then
               flxx(i+1,j,k)=0
               else if (iface(i-1,j,k) .eq. 2) then
               flxx(i,j,k)=0
               end if

               if (iface(i,j+1,k) .eq. 2)then
               flxy(i,j+1,k)=0
               else if (iface(i,j-1,k) .eq. 2) then
               flxy(i,j,k)=0
               end if

               if (iface(i,j,k+1) .eq. 2)then
               flxz(i,j,k+1)=0
               else if (iface(i,j,k-1) .eq. 2) then
               flxz(i,j,k)=0
               end if
           end if 
           uout(i,j,k) = uin(i,j,k) + &
              ( (flxx(i,j,k) - flxx(i+1,j,k)) * dtdx(1) &
                + (flxy(i,j,k) - flxy(i,j+1,k)) * dtdx(2) &
                + (flxz(i,j,k) - flxz(i,j,k+1)) * dtdx(3) )+dt*ptS(i,j,k)
           endif
        enddo
     enddo
  enddo
  conmax_out=maxval(abs(uout))

  ! Scale by face area in order to correctly reflx
  do       k = lo(3), hi(3)
     do    j = lo(2), hi(2)
        do i = lo(1), hi(1)+1
           flxx(i,j,k) = flxx(i,j,k) * (dt * dx(2)*dx(3))
        enddo
     enddo
  enddo
  do       k = lo(3), hi(3)
     do    j = lo(2), hi(2)+1 
        do i = lo(1), hi(1)
           flxy(i,j,k) = flxy(i,j,k) * (dt * dx(1)*dx(3))
        enddo
     enddo
  enddo
  do       k = lo(3), hi(3)+1
     do    j = lo(2), hi(2)
        do i = lo(1), hi(1)
           flxz(i,j,k) = flxz(i,j,k) * (dt * dx(1)*dx(2))
        enddo
     enddo
  enddo

  call bl_deallocate(conx  )
  call bl_deallocate(conx_y)
  call bl_deallocate(conx_z)
  call bl_deallocate(cony  )
  call bl_deallocate(cony_x)
  call bl_deallocate(cony_z)
  call bl_deallocate(conz  )
  call bl_deallocate(conz_x)
  call bl_deallocate(conz_y)
  call bl_deallocate(slope)

end subroutine advect_3d
