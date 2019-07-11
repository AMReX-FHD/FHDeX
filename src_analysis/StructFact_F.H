#ifndef _struct_fact_F_H_
#define _struct_fact_F_H_

#ifdef __cplusplus
extern "C" 
{
#endif


///////////////////////////
// in struct_fact.F90

  void fft_shift(const int* lo, const int* hi,
		 const amrex_real* dft, const int* dftlo, const int* dfthi, const int* ncomp, 
		 const int* zero_avg);

  void sqrt_mf(const int* lo, const int* hi,
	       const amrex_real* mf, const int* mflo, const int* mfhi, const int* ncomp);

#ifdef __cplusplus
}
#endif

#endif
