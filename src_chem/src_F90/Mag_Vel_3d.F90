subroutine get_MagVel_3d( lo, hi, &
     &            u, u_lo, u_hi, &
     &            v, v_lo, v_hi,&
     &            w, w_lo, w_hi,&
     &            MagVel, m_lo, m_hi, & 
     &            zcen, z_lo, z_hi, & 
     &            iface, if_lo, if_hi,&
     &            dx, prob_lo) bind(C, name="get_MagVel_3d")
  


  implicit none
 
  ! problem specifications
  double precision, intent(in) :: dx(3), prob_lo(3)
  ! work region
  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) :: m_lo(3),  m_hi(3)

  integer, intent(in) :: if_lo(3), if_hi(3) 
  integer, intent(in) :: u_lo(3), u_hi(3)
  integer, intent(in) :: v_lo(3), v_hi(3)
  integer, intent(in) :: w_lo(3), w_hi(3)
  integer, intent(in) :: z_lo(3), z_hi(3)

  
  ! ** IN: u,v,w - x, y, z components of the velocity 
  !        iface - location of interface on grid
  double precision, intent(in   ) :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3))
  double precision, intent(in   ) :: v(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3))
  double precision, intent(in   ) :: w(w_lo(1):w_hi(1),w_lo(2):w_hi(2),w_lo(3):w_hi(3))
  
  
  integer, intent(in) :: iface(if_lo(1):if_hi(1),if_lo(2):if_hi(2),if_lo(3):if_hi(3))

  ! ** OUT: MagVel - magnitude of velocity
 double precision, intent(out) :: MagVel(m_lo(1):m_hi(1),m_lo(2):m_hi(2),m_lo(3):m_hi(3))
 double precision, intent(out) :: zcen(z_lo(1):z_hi(1),z_lo(2):z_hi(2),z_lo(3):z_hi(3))

  ! norm - normal of level set defined on cell centers
  ! uvel,vvel,wvel  - x, y, and z componets of velocity used to find the magnitude
  double precision  :: uvel, vvel, wvel
  integer :: i, j, k
do       k = lo(3), hi(3)
     do    j = lo(2), hi(2)
        do i = lo(1), hi(1)
           ! cell centered position in grid           
  
           if (iface(i,j,k)==2) then
           ! we don't care about the non-physical velocity inside the immersed boundary
           MagVel(i,j,k)=0.
           zcen(i,j,k)=0.
           else
           ! make sure velocities are cell centered
           uvel=(u(i-1,j,k)+u(i,j,k))/2
           vvel=(v(i,j-1,k)+v(i,j,k))/2
           wvel=(w(i,j,k-1)+w(i,j,k))/2
           MagVel(i,j,k)=sqrt(uvel**2+vvel**2+wvel**2)
           zcen(i,j,k)=wvel
           end if
        enddo
     enddo
  enddo

end subroutine get_MagVel_3d
