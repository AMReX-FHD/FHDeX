#include "hydro_functions.H"

void ConvertMToUmac(const std::array< MultiFab, AMREX_SPACEDIM >& s_fc,
                    std::array< MultiFab, AMREX_SPACEDIM >& umac,
                    std::array< MultiFab, AMREX_SPACEDIM >& m,
                    int m_to_umac) {

    BL_PROFILE_VAR("ConvertMToUmac()",ConvertMToUmac);

    if (m_to_umac == 1) {
        // compute umac = m / rho - NO GHOST CELLS
        for (int d=0; d<AMREX_SPACEDIM; d++) {
            MultiFab::Copy(   umac[d], m[d],    0, 0, 1, 0);
            MultiFab::Divide( umac[d], s_fc[d], 0, 0, 1, 0);
        }
    }
    else {
        // compute m = rho * umac - INCLUDING GHOST CELLS
        int ng = m[0].nGrow();
        for (int d=0; d<AMREX_SPACEDIM; d++) {
            MultiFab::Copy(     m[d], umac[d], 0, 0, 1, ng);
            MultiFab::Multiply( m[d], s_fc[d], 0, 0, 1, ng);
        }
    }
}
