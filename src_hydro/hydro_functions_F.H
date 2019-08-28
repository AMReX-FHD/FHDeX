#ifndef _hydro_functions_F_H_
#define _hydro_functions_F_H_

#ifdef __cplusplus
extern "C"
{
#endif

    // these are interfaces to fortran subroutines

///////////////////////////
// in hydro_bounds_mod.F90

void set_pressure_bc(const int * lo,     const int * hi,
                     const int * dom_lo, const int * dom_hi,
                     amrex::Real * data, const int * d_lo, const int * d_hi, const int * d_ncomp,
                     int ngc);

///////////////////////////

#ifdef __cplusplus
}
#endif

#endif
