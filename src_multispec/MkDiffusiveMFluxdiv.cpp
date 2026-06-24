#include "multispec_functions.H"
#include "gmres_functions.H"

void MkDiffusiveMFluxdiv(std::array<MultiFab, AMREX_SPACEDIM> & m_update,
                         const std::array<MultiFab, AMREX_SPACEDIM> & umac,
                         const MultiFab& eta,
                         const std::array<MultiFab, NUM_EDGE> & eta_ed,
                         const MultiFab& kappa,
                         const Geometry& geom,
                         const Real* dx,
                         const int& increment)
{
    BL_PROFILE_VAR("MkDiffusiveMFluxdiv()",MkDiffusiveMFluxdiv);

    BoxArray ba = eta.boxArray();
    DistributionMapping dmap = eta.DistributionMap();

    std::array< MultiFab, AMREX_SPACEDIM > alpha_fc;
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        alpha_fc[d].define(convert(ba,nodal_flag_dir[d]), dmap, 1, 0);
        // set alpha_fc to zero
        alpha_fc[d].setVal(0.);
    }

    if (increment == 1) {

        std::array< MultiFab, AMREX_SPACEDIM > Lphi_fc;
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            Lphi_fc[d].define(convert(ba,nodal_flag_dir[d]), dmap, 1, 0);
        }

        // compute -L(phi)
        // we could compute +L(phi) but then we'd have to multiply beta and kappa by -1
        StagApplyOp(geom,eta,kappa,eta_ed,umac,Lphi_fc,alpha_fc,dx,1.);

        // subtract -L(phi) from m_update
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            MultiFab::Subtract(m_update[d],Lphi_fc[d],0,0,1,0);
        }

    }
    else {

        // compute -L(phi)
        // we could compute +L(phi) but then we'd have to multiply beta and kappa by -1
        StagApplyOp(geom,eta,kappa,eta_ed,umac,m_update,alpha_fc,dx,1.);

        // multiply m_update by -1
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            m_update[d].mult(-1,0);
        }
    }
}
