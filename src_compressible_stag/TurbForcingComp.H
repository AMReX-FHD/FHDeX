#ifndef _TurbForcingComp_H_
#define _TurbForcingComp_H_

#include "common_functions.H"

class TurbForcingComp {

public:

    TurbForcingComp();

    void define(BoxArray ba_in, DistributionMapping dmap_in,
                const Real& a_in, const Real& b_in, const Real& c_in, const Real& d_in, const Real& alpha_in);

    void Initialize(const Geometry& geom_in);

    void CalcTurbForcingComp(std::array< MultiFab, AMREX_SPACEDIM >& vel_f,
                        const Real& dt,
                        const int& update);

    std::tuple<amrex::Real, amrex::Real> getU(const int& i);

    void setU(const int& i, Real fs, Real fc);

private:

    Gpu::HostVector<int> KX{1, 0, 0, 1, 1, 0, 1, 2, 0, 0, 2, 2, 1, 0, 1, 0, 2, 1, 1, 2, 2, 0};
    Gpu::HostVector<int> KY{0, 1, 0, 1, 0, 1, 1, 0, 2, 0, 1, 0, 2, 2, 0, 1, 1, 2, 1, 2, 0, 2};
    Gpu::HostVector<int> KZ{0, 0, 1, 0, 1, 1, 1, 0, 0, 2, 0, 1, 0, 1, 2, 2, 1, 1, 2, 0, 2, 2};

    Gpu::HostVector<Real> ForcingS;
    Gpu::HostVector<Real> ForcingC;

    Gpu::HostVector<Real> ForcingSold;
    Gpu::HostVector<Real> ForcingCold;

    std::array< MultiFab, AMREX_SPACEDIM> sines;
    std::array< MultiFab, AMREX_SPACEDIM> cosines;

    Real forcing_a; // 1st coeff for solenoidal forcing
    Real forcing_b; // 2nd coeff for solenoidal forcing
    Real forcing_c; // 1st coeff for dilational forcing
    Real forcing_d; // 2nd coeff for dilational forcing
    Real alpha;     // fraction of solenoidal forcing

};

#endif
