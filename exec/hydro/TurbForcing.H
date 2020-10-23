#ifndef _TurbForcing_H_
#define _TurbForcing_H_

#include "common_functions.H"

class TurbForcing {

    Vector<int> kx{1, 0, 0, 1, 1, 0, 1, 2, 0, 0, 2, 2, 1, 0, 1, 0, 2, 1, 1, 2, 2, 0};
    Vector<int> ky{0, 1, 0, 1, 0, 1, 1, 0, 2, 0, 1, 0, 2, 2, 0, 1, 1, 2, 1, 2, 0, 2};
    Vector<int> kz{0, 0, 1, 0, 1, 1, 1, 0, 0, 2, 0, 1, 0, 1, 2, 2, 1, 1, 2, 0, 2, 2};

    std::array< MultiFab, AMREX_SPACEDIM> sines;
    std::array< MultiFab, AMREX_SPACEDIM> cosines;

    GpuArray<Real,132> rngs;

public:TurbForcing(BoxArray ba_in, DistributionMapping dmap_in, Geometry geom_in);

    void AddTurbForcing(std::array< MultiFab, AMREX_SPACEDIM >& gmres_rhs_u,
                        const Real& dt,
                        const int& reset_RNG);

};

#endif
