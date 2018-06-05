#ifndef _FhdParticleContainer_H_
#define _FhdParticleContainer_H_

#include <AMReX.H>
#include <AMReX_Particles.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ArrayLim.H>

#include "gmres_functions.H"

#include "common_functions.H"

namespace amrex {

class FhdParticleContainer
    : public ParticleContainer<15, 2+BL_SPACEDIM>
{
public:

    using FhdParIter = ParIter<15, 2+BL_SPACEDIM>;

    FhdParticleContainer(const Geometry            & geom, 
                                  const DistributionMapping & dmap,
                                  const BoxArray            & ba);

    void InitParticles();

    void updateParticles(const Real dt, const Real* dx, const std::array<MultiFab, AMREX_SPACEDIM>& umac,
                                           std::array<MultiFab, AMREX_SPACEDIM>& umacNodal,
                                           const std::array<MultiFab, AMREX_SPACEDIM>& RealFaceCoords,
                                           const MultiFab& betaCC, //Not necessary but may use later
                                           MultiFab& betaNodal, //Not necessary but may use later
                                           const MultiFab& rho, //Not necessary but may use later
                                           const std::array<MultiFab, AMREX_SPACEDIM>& source,
                                           const std::array<MultiFab, AMREX_SPACEDIM>& sourceTemp);

    void WriteParticlesAscii(int n);

};
    
}

#endif
