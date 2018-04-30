#ifndef _common_functions_F_H_
#define _common_functions_F_H_

#ifdef __cplusplus
extern "C" 
{
#endif

    // these are interfaces to fortran subroutines
    void read_common_params (const char* probin_file, const int length);

#ifdef __cplusplus
}
#endif

#endif
