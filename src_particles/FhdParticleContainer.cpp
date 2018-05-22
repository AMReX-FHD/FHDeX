#include "FhdParticleContainer.H"
#include "particle_functions_F.H"

using namespace amrex;

FhdParticleContainer::FhdParticleContainer(const Geometry            & geom,
                              const DistributionMapping & dmap,
                              const BoxArray            & ba)
    : ParticleContainer<BL_SPACEDIM> (geom, dmap, ba)
{}

void FhdParticleContainer::InitParticles()
{
    const int lev = 0;
    const Geometry& geom = Geom(lev); //Linking to geom given to constructor?
    const Real* dx  = geom.CellSize();
    
    std::mt19937 mt(0451);
    std::uniform_real_distribution<double> unit(0, 1.0);
    std::uniform_real_distribution<double> dist(-1.0, 1.0);

    int boxes = 0;

    for (MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi) 
    {

        const Box& validBox = mfi.validbox();

        const int* lovect = validBox.loVect();
        const int* hivect = validBox.hiVect();

        boxes++;

        //Print() << "Box: " << boxes << "\n";
        //Print() << "LoCoord: (" << lovect[0]*dx[0] << ", " << lovect[1]*dx[1] << ", " << lovect[2]*dx[2] << ")." << "\n";

        ParticleType p;

        const int grid_id = mfi.index();
        auto& particle_grid = GetParticles(lev)[std::make_pair(grid_id,0)];


//Place one particle in the centre of each cell
/*        for (IntVect iv = validBox.smallEnd(); iv <= validBox.bigEnd(); validBox.next(iv))
        {
            p.id() = ParticleType::NextID();
            p.cpu() = ParallelDescriptor::MyProc();

           // Print() << "Iterator: (" << iv[0] << ", " << iv[1] << ", " << iv[2] << ")." << "\n";

            p.pos(0) = (iv[0] + 0.5)*dx[0];
            p.pos(1) = (iv[1] + 0.5)*dx[1];
#if (BL_SPACEDIM == 3)
            p.pos(2) = (iv[2] + 0.5)*dx[2];
#endif
            p.rdata(0) = dist(mt);
            p.rdata(1) = dist(mt);
            p.rdata(2) = dist(mt);
#if (BL_SPACEDIM == 3)
            p.rdata(3) = dist(mt);
#endif

            particle_grid.push_back(p);

        }*/

//Place 2 particles (per box) randomly in the domain
        for(int i = 0; i<2; i++)
        {
            p.id() = ParticleType::NextID();
            p.cpu() = ParallelDescriptor::MyProc();

            p.pos(0) = (lovect[0]+(hivect[0] - lovect[0])*unit(mt))*dx[0];
            p.pos(1) = (lovect[1]+(hivect[1] - lovect[1])*unit(mt))*dx[1];
#if (BL_SPACEDIM == 3)
            p.pos(2) = (lovect[2]+(hivect[2] - lovect[2])*unit(mt))*dx[2];
#endif
            p.rdata(0) = dist(mt);
            p.rdata(1) = dist(mt);
            p.rdata(2) = dist(mt);
#if (BL_SPACEDIM == 3)
            //p.rdata(3) = dist(mt);
#endif
            particle_grid.push_back(p);

        }



    }

    Redistribute(); //Not necessary?
}

void FhdParticleContainer::moveParticles(const Real dt)
{
    const int lev = 0;
    const RealBox& prob_domain = Geom(lev).ProbDomain();

    for (MyParIter pti(*this, lev); pti.isValid(); ++pti) 
    {
        AoS& parts = pti.GetArrayOfStructs();
        int Np = parts.size();

        Print() << "size: " << Np << "\n";

        move_particles(parts.data(), &Np, &dt, prob_domain.lo(), prob_domain.hi());
    }
}

void FhdParticleContainer::WriteParticlesAscii(int n)
{
    const std::string& pltfile = amrex::Concatenate("particles", n, 5);
    WriteAsciiFile(pltfile);
}

