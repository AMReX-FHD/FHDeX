#ifndef _particle_functions_K_H_
#define _particle_functions_K_H_

#include <AMReX.H>
#include <FhdParticleContainer.H>
#include <common_namespace.H>

/**
   Check if particle is near a boundary. Assumes boundary is parallel to either xy, zx, yz, plane
   
   Inputs : Particle part
            int dir: either +-1,+-2,+-3

   Returns: 
            int near_wall: 0=part's p3m_radius/2 does not overlap with wall, 1=yes  
 */
int near_wall_check (FhdParticleContainer::ParticleType& part, int dir)
{
    return 0;
}

AMREX_GPU_HOST_DEVICE
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
    
    // clearly we don't want to do this every time...
    Gpu::DeviceVector<int> int_tog(common::p_int_tog.size());
    for (std::size_t i = 0; i < common::p_int_tog.size(); ++i)
    {
        int_tog[i] = common::p_int_tog[i];
    }
    auto p_tog = int_tog.data();
    auto rmin = common::rmin;
    
    AMREX_FOR_1D( Np, i, 
    {
        FhdParticleContainer::ParticleType& p1 = pstruct[i];

        for (const auto& p2 : nbor_data.getNeighbors(i))
        {       
            if ( (p_tog[p1.idata(FHD_intData::species)] != 0) and 
                 (p_tog[p2.idata(FHD_intData::species)] != 0) )
            {
                const Real dx = p1.pos(0) - p2.pos(0);
                const Real dy = p1.pos(1) - p2.pos(1);
                const Real dz = p1.pos(2) - p2.pos(2);

                Real r2 = dx*dx + dy*dy + dz*dz;
                r2 = amrex::max(r2, rmin*rmin);
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
