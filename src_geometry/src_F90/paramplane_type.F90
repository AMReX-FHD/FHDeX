module paramplane_type_module
  use amrex_fort_module, only: amrex_real, amrex_particle_real
  use iso_c_binding ,    only: c_int

  implicit none
  private
  
  public paramplane_t
  
  type, bind(C) :: paramplane_t
#if (BL_SPACEDIM == 3)
     real(amrex_particle_real) :: x0
     real(amrex_particle_real) :: y0
     real(amrex_particle_real) :: z0

     real(amrex_particle_real) :: ux
     real(amrex_particle_real) :: uy
     real(amrex_particle_real) :: uz

     real(amrex_particle_real) :: vx
     real(amrex_particle_real) :: vy
     real(amrex_particle_real) :: vz  


     real(amrex_particle_real) :: lnx
     real(amrex_particle_real) :: lny
     real(amrex_particle_real) :: lnz

     real(amrex_particle_real) :: rnx
     real(amrex_particle_real) :: rny
     real(amrex_particle_real) :: rnz

     real(amrex_particle_real) :: utop
     real(amrex_particle_real) :: vtop

     real(amrex_particle_real) :: costhetaleft
     real(amrex_particle_real) :: sinthetaleft
     real(amrex_particle_real) :: cosphileft
     real(amrex_particle_real) :: sinphileft

     real(amrex_particle_real) :: costhetaright
     real(amrex_particle_real) :: sinthetaright
     real(amrex_particle_real) :: cosphiright
     real(amrex_particle_real) :: sinphiright

     real(amrex_particle_real) :: velx
     real(amrex_particle_real) :: vely
     real(amrex_particle_real) :: velz

     real(amrex_particle_real) :: c0
     real(amrex_particle_real) :: resomg
     real(amrex_particle_real) :: agraph
     real(amrex_particle_real) :: bgraph
     real(amrex_particle_real) :: a0graph
     real(amrex_particle_real) :: b0graph
     real(amrex_particle_real) :: coltime

     real(amrex_particle_real) :: besslist(100);
     real(amrex_particle_real) :: dbesslist(100);

#endif

#if (BL_SPACEDIM == 2)
     real(amrex_particle_real) :: x0
     real(amrex_particle_real) :: y0

     real(amrex_particle_real) :: ux
     real(amrex_particle_real) :: uy

     real(amrex_particle_real) :: lnx
     real(amrex_particle_real) :: lny

     real(amrex_particle_real) :: rnx
     real(amrex_particle_real) :: rny

     real(amrex_particle_real) :: utop

     real(amrex_particle_real) :: costhetaleft
     real(amrex_particle_real) :: sinthetaleft
     real(amrex_particle_real) :: costhetaright
     real(amrex_particle_real) :: sinthetaright

#endif

     real(amrex_particle_real) :: porosityleft
     real(amrex_particle_real) :: specularityleft
     real(amrex_particle_real) :: temperatureleft
     real(amrex_particle_real) :: momentumleft

     real(amrex_particle_real) :: porosityright
     real(amrex_particle_real) :: specularityright
     real(amrex_particle_real) :: temperatureright
     real(amrex_particle_real) :: momentumright

     real(amrex_particle_real) :: periodicity

     real(amrex_particle_real) :: fxleft
     real(amrex_particle_real) :: fyleft
     real(amrex_particle_real) :: fzleft

     real(amrex_particle_real) :: fxright
     real(amrex_particle_real) :: fyright
     real(amrex_particle_real) :: fzright

     real(amrex_particle_real) :: fxleftav
     real(amrex_particle_real) :: fyleftav
     real(amrex_particle_real) :: fzleftav

     real(amrex_particle_real) :: fxrightav
     real(amrex_particle_real) :: fyrightav
     real(amrex_particle_real) :: fzrightav

     integer :: boundary

     real(amrex_particle_real) :: tension
     real(amrex_particle_real) :: stiffness
     real(amrex_particle_real) :: temperature

  end type paramplane_t

end module paramplane_type_module
