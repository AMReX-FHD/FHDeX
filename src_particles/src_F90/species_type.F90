module species_type_module
  use amrex_fort_module, only: amrex_real
  use iso_c_binding ,    only: c_int
  
  implicit none
  private
  
  public species_t
  
  type, bind(C) :: species_t

     integer(c_int)            :: type         
     integer(c_int)            :: total        
     integer(c_int)            :: ppb     

     real(amrex_real) :: m
     real(amrex_real) :: d
     real(amrex_real) :: T
     real(amrex_real) :: R
     real(amrex_real) :: q
     real(amrex_real) :: gamma1
     real(amrex_real) :: mu
     real(amrex_real) :: n0

     real(amrex_real) :: mfp
     real(amrex_real) :: P
     real(amrex_real) :: Neff
     real(amrex_real) :: cp
     real(amrex_real) :: propulsion
     real(amrex_real) :: total_diff
     real(amrex_real) :: dry_diff
     real(amrex_real) :: wet_diff
     real(amrex_real) :: sigma
     real(amrex_real) :: eepsilon

  end type species_t

contains
  
end module species_type_module
