#ifndef _gmres_functions_F_H_
#define _gmres_functions_F_H_

#ifdef __cplusplus
extern "C" 
{
#endif

    // these are interfaces to fortran subroutines
    void read_gmres_params ();
    void finalize_gmres_params ();

#ifdef __cplusplus
}
#endif

#endif
