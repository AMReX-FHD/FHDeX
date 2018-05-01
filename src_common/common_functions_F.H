#ifndef _common_functions_F_H_
#define _common_functions_F_H_

#ifdef __cplusplus
extern "C" 
{
#endif

    // these are interfaces to fortran subroutines

    // read in fortran namelists into common_params_module
    void read_common_params (const char* probin_file, const int length);

    // copy contents of common_params_module to C++ common namespace
    void copy_common_params_to_c (amrex::Real* prob_lo, amrex::Real* prob_hi);

#ifdef __cplusplus
}
#endif

#endif
