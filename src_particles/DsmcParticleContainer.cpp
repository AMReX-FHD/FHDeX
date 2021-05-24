#include "DsmcParticleContainer.H"

#include "particle_functions_K.H"
#include "paramplane_functions_K.H"
#include <math.h>


DsmcParticleContainer::DsmcParticleContainer(const Geometry & geom,
                              const DistributionMapping & dmap,
                              const BoxArray            & ba,
                              int ncells)
    : NeighborParticleContainer<Dsmc_realData::count, Dsmc_intData::count> (geom, dmap, ba, ncells)
{
    BL_PROFILE_VAR("DsmcParticleContainer()",DsmcParticleContainer);
}


