#include "DsmcParticleContainer.H"
#include <math.h>

void FhdParticleContainer::externalForce(Real dt)
{
    int lev = 0;

    const Real* dx = Geom(lev).CellSize();
    const Real* plo = Geom(lev).ProbLo();
    const Real* phi = Geom(lev).ProbHi();

    Real meanDrag, fluctDrag;

    for (FhdParIter pti(* this, lev); pti.isValid(); ++pti) {

        const int grid_id = pti.index();
        const int tile_id = pti.LocalTileIndex();
        const Box& tile_box  = pti.tilebox();

        auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];
        auto& particles = particle_tile.GetArrayOfStructs();
        const long np = particles.numParticles();

        for (int i = 0; i < np; ++ i)
        {
            ParticleType & part = particles[i];

            part.rdata(FHD_realData::velx) += (grav[0]*dt);
            part.rdata(FHD_realData::vely) += (grav[1]*dt);
            part.rdata(FHD_realData::velz) += (grav[2]*dt);

        }
    }
}