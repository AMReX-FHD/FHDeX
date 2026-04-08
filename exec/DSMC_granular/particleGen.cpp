#include "INS_functions.H"
#include "common_functions.H"
#include "DsmcParticleContainer.H"
#include <sstream>
#include <string>
#include <fstream>

using namespace std;

void FhdParticleContainer::InitParticles(Real & dt)
{
    const int lev = 0;
    const Geometry& geom = Geom(lev);

    int pcount = 0;

    bool proc0_enter = true;

    std::array<Real, 3> vpart = {0., 0., 0.};
    Real spdmax = 0;
    Real umax = 0, vmax = 0, wmax = 0;
    Real spd;
    Real u[nspecies],v[nspecies],w[nspecies];

    for (MFIter mfi = MakeMFIter(lev, true); mfi.isValid(); ++mfi)
    {
        const Box& tile_box  = mfi.tilebox();
        const RealBox tile_realbox{tile_box, geom.CellSize(), geom.ProbLo()};
        const int grid_id = mfi.index();
        const int tile_id = mfi.LocalTileIndex();
        auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];
        auto& particles = particle_tile.GetArrayOfStructs();

        IntVect smallEnd = tile_box.smallEnd();
        IntVect bigEnd = tile_box.bigEnd();

        if(ParallelDescriptor::MyProc() == 0 && mfi.LocalTileIndex() == 0 && proc0_enter)
        {
            proc0_enter = false;
            if(particle_input>0)
            {
                std::ifstream particleFile("particles.dat");
                while(true)
                {
                    ParticleType p;
                    p.id()  = ParticleType::NextID();
                    p.cpu() = ParallelDescriptor::MyProc();
                    p.idata(FHD_intData::sorted) = -1;
                    particleFile >> p.pos(0);
                    particleFile >> p.pos(1);
                    particleFile >> p.pos(2);
                    particleFile >> vpart[0];
                    particleFile >> vpart[1];
                    particleFile >> vpart[2];
                    p.rdata(FHD_realData::velx) = vpart[0];
                    p.rdata(FHD_realData::vely) = vpart[1];
                    p.rdata(FHD_realData::velz) = vpart[2];
                    p.rdata(FHD_realData::timeFrac) = 1;

                    spd = sqrt(pow(vpart[0],2)+pow(vpart[1],2)+pow(vpart[2],2));
                    if(spd>spdmax){ spdmax=spd; }
                    // For calculating timstep from Courant number
                    if(vpart[0]>umax) { umax=vpart[0]; }
                    if(vpart[1]>vmax) { vmax=vpart[1]; }
                    if(vpart[2]>wmax) { wmax=vpart[2]; }

                    particleFile >> p.idata(FHD_intData::species);
                    int i_spec = p.idata(FHD_intData::species);
                    p.rdata(FHD_realData::R) = k_B/properties[i_spec].mass;
                    p.rdata(FHD_realData::boostx) = 0;
                    p.rdata(FHD_realData::boosty) = 0;
                    p.rdata(FHD_realData::boostz) = 0;
                    if( particleFile.eof() ) break;
                }
                particleFile.close();
            }
            else
            {
                // Initialize to bulk velocities
                for(int i_spec=0; i_spec < nspecies; i_spec++)
                {
                    u[i_spec] = 0.0; v[i_spec] = 0.0; w[i_spec] = 0.0;
                }

                for(int i_spec=0; i_spec < nspecies; i_spec++)
                {
                    // Standard deviation of velocity at temperature T_init
                    Real R     = k_B/properties[i_spec].mass;
                    Real stdev = sqrt(T_init[i_spec]*R);
                    for (int i_part=0; i_part<properties[i_spec].total;i_part++)
                    {
                        ParticleType p;
                        p.id()  = ParticleType::NextID();
                        p.cpu() = ParallelDescriptor::MyProc();
                        p.idata(FHD_intData::sorted) = -1;
                        p.idata(FHD_intData::species) = i_spec;
                        p.rdata(FHD_realData::R) = R;
                        p.rdata(FHD_realData::timeFrac) = 1;
                        p.pos(0) = prob_lo[0] + amrex::Random()*(prob_hi[0]-prob_lo[0]);
                        p.pos(1) = prob_lo[1] + amrex::Random()*(prob_hi[1]-prob_lo[1]);
                        p.pos(2) = prob_lo[2] + amrex::Random()*(prob_hi[2]-prob_lo[2]);
                        p.idata(FHD_intData::i) = -100;
                        p.idata(FHD_intData::j) = -100;
                        p.idata(FHD_intData::k) = -100;

                        vpart[0] = stdev*amrex::RandomNormal(0.,1.);
                        vpart[1] = stdev*amrex::RandomNormal(0.,1.);
                        vpart[2] = stdev*amrex::RandomNormal(0.,1.);
                        p.rdata(FHD_realData::velx) = vpart[0];
                        p.rdata(FHD_realData::vely) = vpart[1];
                        p.rdata(FHD_realData::velz) = vpart[2];
                        spd = sqrt(pow(vpart[0],2)+pow(vpart[1],2)+pow(vpart[2],2));

                        if(spd>spdmax){ spdmax=spd; }
                        // For calculating timstep from Courant number
                        if(vpart[0]>umax) { umax=vpart[0]; }
                        if(vpart[1]>vmax) { vmax=vpart[1]; }
                        if(vpart[2]>wmax) { wmax=vpart[2]; }

                        u[p.idata(FHD_intData::species)] += vpart[0];
                        v[p.idata(FHD_intData::species)] += vpart[1];
                        w[p.idata(FHD_intData::species)] += vpart[2];

                        p.rdata(FHD_realData::boostx) = 0;
                        p.rdata(FHD_realData::boosty) = 0;
                        p.rdata(FHD_realData::boostz) = 0;
                        particle_tile.push_back(p);
                        pcount++;
                    }
                }

                // Zero out the bulk velocities
                // Only do if no restart provided for particles
                // VCOM by species
                for(int i_spec=0; i_spec < nspecies; i_spec++)
                {
                    long npi = properties[i_spec].total;
                    u[i_spec] = u[i_spec]/npi;
                    v[i_spec] = v[i_spec]/npi;
                    w[i_spec] = w[i_spec]/npi;
                }

                const int np = particles.numParticles();
                for(int i = 0; i < np; i++)
                {
                    ParticleType & part = particles[i];
                    int i_spec = part.idata(FHD_intData::species);
                    part.rdata(FHD_realData::velx) = part.rdata(FHD_realData::velx)-u[i_spec];
                    part.rdata(FHD_realData::vely) = part.rdata(FHD_realData::vely)-v[i_spec];
                    part.rdata(FHD_realData::velz) = part.rdata(FHD_realData::velz)-w[i_spec];
                }
            }
        }
    }
    // Set guess of max relative velocity
    ParallelDescriptor::Bcast(&spdmax,1,ParallelDescriptor::IOProcessorNumber());
    mfvrmax.setVal(spdmax);

    if(ParallelDescriptor::MyProc() == 0) {
        dt  = umax*n_cells[0]/(prob_hi[0]-prob_lo[0]);
        dt += vmax*n_cells[1]/(prob_hi[1]-prob_lo[1]);
        dt += wmax*n_cells[2]/(prob_hi[2]-prob_lo[2]);
        dt  = 0.2/dt; // Courant number of 0.2
    }
    ParallelDescriptor::Bcast(&dt,1,ParallelDescriptor::IOProcessorNumber());
    amrex::Print() << "My dt " << dt << "\n";

    Redistribute();
    SortParticles();

    // Zero bulk velocities in each cell
    for (FhdParIter pti(* this, lev); pti.isValid(); ++pti)
    {
        const int grid_id = pti.index();
        const int tile_id = pti.LocalTileIndex();
        const Box& tile_box  = pti.tilebox();

        auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];
        auto& particles = particle_tile.GetArrayOfStructs();
        const long np = particles.numParticles();

        IntVect smallEnd = tile_box.smallEnd();
        IntVect bigEnd = tile_box.bigEnd();

        for (int i = smallEnd[0]; i <= bigEnd[0]; i++) {
        for (int j = smallEnd[1]; j <= bigEnd[1]; j++) {
        for (int k = smallEnd[2]; k <= bigEnd[2]; k++) {
            const IntVect& iv = {i,j,k};
            long imap = tile_box.index(iv);
            for (int niter=0; niter<3; niter++) {
            for (int ispec=0; ispec<nspecies; ispec++){
                Real ucom=0., vcom=0., wcom=0.;
                int np = m_cell_vectors[ispec][grid_id][imap].size();
                double lmass = properties[ispec].mass;
                for (int ip = 0; ip<np; ip++) {
                    int ipart = m_cell_vectors[ispec][grid_id][imap][ip];
                    ParticleType & part = particles[ipart];
                    ucom += part.rdata(FHD_realData::velx);
                    vcom += part.rdata(FHD_realData::vely);
                    wcom += part.rdata(FHD_realData::velz);
                }
                ucom /= (double)np;
                vcom /= (double)np;
                wcom /= (double)np;
                for (int ip = 0; ip<np; ip++)
                {
                    int ipart = m_cell_vectors[ispec][grid_id][imap][ip];
                    ParticleType & part = particles[ipart];
                    part.rdata(FHD_realData::velx) = part.rdata(FHD_realData::velx) - ucom;
                    part.rdata(FHD_realData::vely) = part.rdata(FHD_realData::vely) - vcom;
                    part.rdata(FHD_realData::velz) = part.rdata(FHD_realData::velz) - wcom;
                }
            }
            }
        }
        }
        }
    }
}

void FhdParticleContainer::ReInitParticles()
{
    const int lev = 0;
    const Geometry& geom = Geom(lev);

    int pcount = 0;
    bool proc0_enter = true;

    std::array<Real, 3> vpart = {0.,0.,0.};
    Real u,v,w;
    Real spd, spdmax = 0.;

    for (MFIter mfi = MakeMFIter(lev, true); mfi.isValid(); ++mfi)
    {
        const Box& tile_box = mfi.tilebox();
        const RealBox tile_realbox{tile_box, geom.CellSize(), geom.ProbLo()};
        const int grid_id = mfi.index();
        const int tile_id = mfi.LocalTileIndex();
        auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];
        auto& particles = particle_tile.GetArrayOfStructs();
        const long np = particles.numParticles();

        for (int i=0; i<np; ++i)
        {
            ParticleType & part = particles[i];
            part.id() = ParticleType::NextID();
            part.cpu() = ParallelDescriptor::MyProc();
            part.idata(FHD_intData::sorted) = -1;
            part.idata(FHD_intData::i) = -100;
            part.idata(FHD_intData::j) = -100;
            part.idata(FHD_intData::k) = -100;

            part.rdata(FHD_realData::timeFrac) = 1;
            Real u = part.rdata(FHD_realData::velx);
            Real v = part.rdata(FHD_realData::vely);
            Real w = part.rdata(FHD_realData::velz);
            spd = sqrt(pow(u,2.0)+pow(v,2.0)+pow(w,2.0));
            if(spd>spdmax)
            {
                spdmax=spd;
            }
        }
    }

    ParallelDescriptor::Bcast(&spdmax,1,ParallelDescriptor::IOProcessorNumber());
    mfvrmax.setVal(spdmax);

    Redistribute();
    SortParticles();
}