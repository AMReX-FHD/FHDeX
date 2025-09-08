#ifndef _StochMomFlux_H_
#define _StochMomFlux_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Vector.H>

#include "common_functions.H"

class StochMomFlux {

    const int ncomp_ed = 2;
    int n_rngs;

    Geometry geom;

    // Vectors containing Gaussian white noise field at each stage:
    // the vector is for "n_rngs", which is the number of random number stages
    // mflux cell-centered
    Vector< MultiFab > mflux_cc;
    // mflux on nodes in 2d, on edges in 3d
    Vector<std::array< MultiFab, NUM_EDGE > >  mflux_ed;

    // temporary storage of linear combinations of random number stages
    // multiplied by variance*sqrt(eta*temperature):
    MultiFab mflux_cc_weighted;
    std::array< MultiFab, NUM_EDGE >  mflux_ed_weighted;

public:

    // initialize n_rngs, geom
    // build MultiFabs to hold random numbers
    StochMomFlux(BoxArray,DistributionMapping,Geometry,int);

    // fill mflux_cc and mflux_ed with random numbers
    void fillMomStochastic();

    // create weighted sum of stage RNGs and store in mflux_cc_weighted and mflux_ed_weighted
    void weightMomflux(Vector<amrex::Real>);

    // scale random numbers that lie on physical boundaries appropriately
    void MomFluxBC();

    // Multiply mflux_weighted by sqrt(eta*temperature)
    void multbyVarSqrtEtaTemp(const amrex::MultiFab&, const std::array< MultiFab, NUM_EDGE >&,
                              const amrex::MultiFab&, const std::array< MultiFab, NUM_EDGE >&,
                              const amrex::Real&);

    // compute stochastic momentum flux divergence
    void StochMomFluxDiv(std::array< amrex::MultiFab, AMREX_SPACEDIM >&,
                         const int&, const amrex::MultiFab&,
                         const std::array< amrex::MultiFab, NUM_EDGE >&,
                         const amrex::MultiFab&,
                         const std::array< amrex::MultiFab, NUM_EDGE >&,
                         const amrex::Vector< amrex::Real >&, const amrex::Real&);

    void StochMomFluxDivWide(std::array< amrex::MultiFab, AMREX_SPACEDIM >&,
                         const int&, const amrex::MultiFab&,
                         const std::array< amrex::MultiFab, NUM_EDGE >&,
                         const amrex::MultiFab&,
                         const std::array< amrex::MultiFab, NUM_EDGE >&,
                         const amrex::Vector< amrex::Real >&, const amrex::Real&);

    void StochMomFluxDivWideSplit(std::array< amrex::MultiFab, AMREX_SPACEDIM >&,
                         const int&, const amrex::MultiFab&,
                         const std::array< amrex::MultiFab, NUM_EDGE >&,
                         const amrex::MultiFab&,
                         const std::array< amrex::MultiFab, NUM_EDGE >&,
                         const amrex::Vector< amrex::Real >&, const amrex::Real&);

                             // compute stochastic momentum flux divergence
    void StochMomFluxDivOrder3(std::array< amrex::MultiFab, AMREX_SPACEDIM >&,
                         const int&, const amrex::MultiFab&,
                         const std::array< amrex::MultiFab, NUM_EDGE >&,
                         const amrex::MultiFab&,
                         const std::array< amrex::MultiFab, NUM_EDGE >&,
                         const amrex::Vector< amrex::Real >&, const amrex::Real&);

    void StochMomFluxDivOrder3Split(std::array< amrex::MultiFab, AMREX_SPACEDIM >&,
                         const int&, const amrex::MultiFab&,
                         const std::array< amrex::MultiFab, NUM_EDGE >&,
                         const amrex::MultiFab&,
                         const std::array< amrex::MultiFab, NUM_EDGE >&,
                         const amrex::Vector< amrex::Real >&, const amrex::Real&);

    // utility to write out random number MultiFabs to plotfiles
    void writeMFs(std::array< MultiFab, AMREX_SPACEDIM >&);

};

#endif
