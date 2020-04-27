#ifndef _particle_functions_K_H_
#define _particle_functions_K_H_

#include <AMReX.H>
#include <FhdParticleContainer.H>
#include <common_namespace.H>

void repulsive_force (      FhdParticleContainer::ParticleType& part1,
                      const FhdParticleContainer::ParticleType& part2,
                      amrex::Real dx, amrex::Real dy, amrex::Real dz, amrex::Real dr2)
{
    using namespace amrex;

    const Real eps1 = part1.rdata(FHD_realData::eepsilon);
    const Real sig1 = part1.rdata(FHD_realData::sigma);
    const Real sig2 = part2.rdata(FHD_realData::sigma);
    const Real fac1 = std::pow((sig1 + sig2)/4.0, 12);
    const Real fac2 = std::pow((sig1 + sig2)/4.0, 6);
    const Real ff = eps1*4*(1./(dr2*dr2*dr2*dr2))*(-12*fac1/(dr2*dr2*dr2) + 6*fac2);
    
    part1.rdata(FHD_realData::forcex) -= dx*ff;
    part1.rdata(FHD_realData::forcey) -= dy*ff;
    part1.rdata(FHD_realData::forcez) -= dz*ff;
}

void compute_forces_nl (FhdParticleContainer::AoS& aos, int Np, int Nn,
                        amrex::NeighborList<FhdParticleContainer::ParticleType>& neighbor_list,
                        amrex::Real& rcount)
{
    using namespace amrex;

    auto nbor_data = neighbor_list.data();
    FhdParticleContainer::ParticleType* pstruct = aos().dataPtr();
    
    // need to access rmin, p_int_tog,

    AMREX_FOR_1D( Np, i, 
    {
        FhdParticleContainer::ParticleType& p1 = pstruct[i];

        for (const auto& p2 : nbor_data.getNeighbors(i))
        {       
            if ( (common::p_int_tog[p1.rdata(FHD_intData::species)] != 0) and 
                 (common::p_int_tog[p2.rdata(FHD_intData::species)] != 0) )
            {
                const Real dx = p1.pos(0) - p2.pos(0);
                const Real dy = p1.pos(1) - p2.pos(1);
                const Real dz = p1.pos(2) - p2.pos(2);
                
                Real r2 = dx*dx + dy*dy + dz*dz;
                r2 = amrex::max(r2, common::rmin*common::rmin);
                const Real r = sqrt(r2);
                
                if (r < (1.122*p1.rdata(FHD_realData::sigma)/2.0))
                {
                    rcount += 1;  // needs porting...
                    repulsive_force(p1, p2, dx, dy, dz, r2);
                }
            }
        }
    });
}

#endif
