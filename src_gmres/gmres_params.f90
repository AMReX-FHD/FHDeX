module gmres_params_module

  implicit none

  ! Begin the declarations of the ParmParse parameters

  integer, save :: precon_type

  ! End the declarations of the ParmParse parameters

contains

  subroutine read_gmres_params() bind(C, name="read_gmres_params")

    use amrex_parmparse_module, only: amrex_parmparse_build, &
                                      amrex_parmparse_destroy, &
                                      amrex_parmparse

    type (amrex_parmparse) :: pp

    ! default values

    ! read in from inputs file
    
    call amrex_parmparse_build(pp)

    call pp%query("precon_type",precon_type);

    call amrex_parmparse_destroy(pp)

  end subroutine read_gmres_params

end module gmres_params_module
