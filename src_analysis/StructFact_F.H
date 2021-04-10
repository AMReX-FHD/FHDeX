#ifndef _struct_fact_F_H_
#define _struct_fact_F_H_

#ifdef __cplusplus
extern "C" 
{
#endif


///////////////////////////
// in struct_fact.F90

  void fft_shift(const int* lo, const int* hi,
		 const amrex_real* dft, const int* dftlo, const int* dfthi,
                 const amrex_real* dft_temp, const int* dft_templo, const int* dft_temphi,
		 const int* zero_avg);


#ifdef __cplusplus
}
#endif

#endif
