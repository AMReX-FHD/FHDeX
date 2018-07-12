#ifndef _hydro_functions_F_H_
#define _hydro_functions_F_H_

#ifdef __cplusplus
extern "C" 
{
#endif

    // these are interfaces to fortran subroutines

///////////////////////////
// in mk_advective_m_fluxdiv.f90

  subroutine mk_advective_m_fluxdiv(lo,hi, &
				    umac, umaclo, umachi, &
				    vmac, vmaclo, vmachi, &
#if (AMREX_SPACEDIM == 3)
				    wmac, wmaclo, wmachi, &
#endif
				    mx, mxlo, mxhi, &
				    my, mylo, myhi, &
#if (AMREX_SPACEDIM == 3)
				    mz, mzlo, mzhi, &
#endif
				    m_updatex, mudxlo, mudxhi, &
				    m_updatey, mudylo, mudyhi, &
#if (AMREX_SPACEDIM == 3)
				    m_updatez, mudzlo, mudzhi, &
#endif
				    dx, increment);

///////////////////////////



#ifdef __cplusplus
}
#endif

#endif
