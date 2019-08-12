#ifndef _debug_functions_F_H_
#define _debug_functions_F_H_

#ifdef __cplusplus
extern "C" 
{
#endif

    // these are interfaces to fortran subroutines

///////////////////////////
// in DEBUG.F90

void print_potential(const int* lo, const int* hi, 
		const amrex_real* pot, const int* potlo, const int* pothi, 
		const int *xstart, const int *ystart, const int *zstart); 

///////////////////////////



#ifdef __cplusplus
}
#endif

#endif
