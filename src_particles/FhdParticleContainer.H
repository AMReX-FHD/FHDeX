#ifndef _FhdParticleContainer_H_
#define _FhdParticleContainer_H_

#include <AMReX.H>
#include <AMReX_Particles.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include "AMReX_NeighborParticles.H"
#include <AMReX_ArrayLim.H>

namespace amrex {

class FhdParticleContainer
    : public ParticleContainer<BL_SPACEDIM + 1, 0>
{
public:

    FhdParticleContainer(const Geometry            & geom, 
                                  const DistributionMapping & dmap,
                                  const BoxArray            & ba);

    ///
    /// Init one particle per cell with random velocities
    ///
    void InitParticles();

};
    
}

#endif // _NEIGHBORLISTPARTICLECONTAINER_H_
