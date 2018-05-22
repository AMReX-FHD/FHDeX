#ifndef _FhdParticleContainer_H_
#define _FhdParticleContainer_H_

#include <AMReX.H>
#include <AMReX_Particles.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ArrayLim.H>

namespace amrex {

class FhdParticleContainer
    : public ParticleContainer<BL_SPACEDIM>
{
public:

    using MyParIter = ParIter<BL_SPACEDIM>;

    FhdParticleContainer(const Geometry            & geom, 
                                  const DistributionMapping & dmap,
                                  const BoxArray            & ba);

    void InitParticles();

    void moveParticles(const Real dt);

    void WriteParticlesAscii(int n);

};
    
}

#endif
