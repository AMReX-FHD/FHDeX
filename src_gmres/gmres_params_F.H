#ifndef _gmres_params_F_H_
#define _gmres_params_F_H_

#ifdef __cplusplus
extern "C" 
{
#endif

    // these are interfaces to fortran subroutines

    // read in fortran namelist into common_params_module
    void read_gmres_params (const char* inputs_file, const int length);

    // copy contents of common_params_module to C++ common namespace
    void copy_gmres_params_to_c ();

#ifdef __cplusplus
}
#endif

#endif
