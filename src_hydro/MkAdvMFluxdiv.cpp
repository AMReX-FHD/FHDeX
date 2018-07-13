#include "hydro_functions.H"
#include "hydro_functions_F.H"

#include "common_functions.H"
#include "common_functions_F.H"
#include "common_namespace.H"

using namespace common;

void MkAdvMFluxdiv(const std::array<MultiFab, AMREX_SPACEDIM>& umac,
		   const std::array<MultiFab, AMREX_SPACEDIM>& m,
		   std::array<MultiFab, AMREX_SPACEDIM>& m_update,
		   const amrex::Real* dx,
		   const bool& increment)
{
    // Loop over boxes
    for (MFIter mfi(umac[0]); mfi.isValid(); ++mfi) {

        // Create cell-centered box from semi-nodal box
        const Box& validBox_cc = enclosedCells(mfi.validbox());

        mk_advective_m_fluxdiv(ARLIM_3D(validBox_cc.loVect()), ARLIM_3D(validBox_cc.hiVect()),
        		       BL_TO_FORTRAN_ANYD(umac[0][mfi]),
			       BL_TO_FORTRAN_ANYD(umac[1][mfi]),
#if (AMREX_SPACEDIM == 3)
			       BL_TO_FORTRAN_ANYD(umac[2][mfi]),
#endif
        		       BL_TO_FORTRAN_ANYD(m[0][mfi]),
			       BL_TO_FORTRAN_ANYD(m[1][mfi]),
#if (AMREX_SPACEDIM == 3)
			       BL_TO_FORTRAN_ANYD(m[2][mfi]),
#endif
        		       BL_TO_FORTRAN_ANYD(m_update[0][mfi]),
			       BL_TO_FORTRAN_ANYD(m_update[1][mfi]),
#if (AMREX_SPACEDIM == 3)
			       BL_TO_FORTRAN_ANYD(m_update[2][mfi]),
#endif
			       dx, &increment);
    }

}
