#ifndef _compressible_test_functions_F_H_
#define _compressible_test_functions_F_H_

#ifdef __cplusplus
extern "C" 
{
#endif

  // these are interfaces to fortran subroutines

///////////////////////////
// in init_consvar.F90

  void init_consvar(const int* lo, const int* hi, 
  		    const amrex_real* cu, const int* culo, const int* cuhi,
		    const amrex_real* pu, const int* pulo, const int* puhi,
  		    const amrex_real* dx, 
  		    const amrex_real* reallo, const amrex_real* realhi);

///////////////////////////


#ifdef __cplusplus
}
#endif

#endif
