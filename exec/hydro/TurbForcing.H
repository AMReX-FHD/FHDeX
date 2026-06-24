#ifndef _TurbForcing_H_
#define _TurbForcing_H_

#include "common_functions.H"

class TurbForcing {

    std::array< MultiFab, AMREX_SPACEDIM> sines;
    std::array< MultiFab, AMREX_SPACEDIM> cosines;

    GpuArray<Real,132> forcing_U;

    Real forcing_a;
    Real forcing_b;

public:

    TurbForcing();

    void define(BoxArray ba_in, DistributionMapping dmap_in,
                const Real& a_in, const Real& b_in);

    void Initialize(const Geometry& geom_in);

    void AddTurbForcing(std::array< MultiFab, AMREX_SPACEDIM >& gmres_rhs_u,
                        const Real& dt,
                        const int& update_U);

    Real getU(const int& i);

    void setU(const int& i, Real x);

};

#endif
