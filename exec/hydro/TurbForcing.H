#ifndef _TurbForcing_H_
#define _TurbForcing_H_

#include "common_functions.H"

class TurbForcing {

    GpuArray<int,22> kx{1, 0, 0, 1, 1, 0, 1, 2, 0, 0, 2, 2, 1, 0, 1, 0, 2, 1, 1, 2, 2, 0};
    GpuArray<int,22> ky{0, 1, 0, 1, 0, 1, 1, 0, 2, 0, 1, 0, 2, 2, 0, 1, 1, 2, 1, 2, 0, 2};
    GpuArray<int,22> kz{0, 0, 1, 0, 1, 1, 1, 0, 0, 2, 0, 1, 0, 1, 2, 2, 1, 1, 2, 0, 2, 2};

    std::array< MultiFab, AMREX_SPACEDIM> sines;
    std::array< MultiFab, AMREX_SPACEDIM> cosines;

    GpuArray<Real,132> forcing_U;

    Real forcing_a;
    Real forcing_b;

public:


    TurbForcing(BoxArray ba_in, DistributionMapping dmap_in,
                const Real& a_in, const Real& b_in);

    void Initialize(Geometry geom_in);

    void AddTurbForcing(std::array< MultiFab, AMREX_SPACEDIM >& gmres_rhs_u,
                        const Real& dt,
                        const int& update_U);

    Real getU(const int& i);

    void setU(const int& i, Real x);

};

#endif
