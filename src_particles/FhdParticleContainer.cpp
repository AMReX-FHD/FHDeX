#include "FhdParticleContainer.H"

using namespace amrex;

FhdParticleContainer::FhdParticleContainer(const Geometry            & geom,
                              const DistributionMapping & dmap,
                              const BoxArray            & ba)
    : ParticleContainer<BL_SPACEDIM + 1, 0> 
      (geom, dmap, ba)
{}

void FhdParticleContainer::InitParticles()
{
    const int lev = 0;
    const Geometry& geom = Geom(lev); //Linking to geom given to constructor?
    const Real* dx  = geom.CellSize();
    
    std::mt19937 mt(0451);
    std::uniform_real_distribution<double> dist(0, 1.0);

    for (MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi) 
    {

        const Box& validBox = mfi.validbox();

        const int* lovect = validBox.loVect();

        ParticleType p;

        const int grid_id = mfi.index();
        auto& particle_grid = GetParticles(lev)[std::make_pair(grid_id,0)];

        for (IntVect iv = validBox.smallEnd(); iv <= validBox.bigEnd(); validBox.next(iv))
        {
            p.id() = ParticleType::NextID();
            p.cpu() = ParallelDescriptor::MyProc();

            p.pos(0) = (lovect[0] + iv[0] + 0.5)*dx[0];
            p.pos(1) = (lovect[1] + iv[1] + 0.5)*dx[0];
#if (BL_SPACEDIM == 3)
            p.pos(2) = (lovect[0] + iv[0] + 0.5)*dx[0];
#endif

            p.rdata(0) = 1;
            p.rdata(1) = 1;
            p.rdata(2) = 1;
#if (BL_SPACEDIM == 3)
            p.rdata(3) = 1;
#endif

            particle_grid.push_back(p);

        }

    }
}


