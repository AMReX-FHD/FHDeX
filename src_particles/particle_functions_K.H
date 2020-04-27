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
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int near_wall_check (const FhdParticleContainer::ParticleType& part, int dir)
{
    using namespace amrex;
    using common::prob_hi;

    const int abs_dir = std::abs(dir);
    const Real p3m_rad = part.rdata(FHD_realData::p3m_radius);
    const int index = abs_dir - 1;

    int near_wall;
    
    if ((abs_dir == 1) or (abs_dir == 2) or (abs_dir == 3))
    {
        if (dir > 0)
        {
            if ((part.pos(index) + 0.5*p3m_rad) > prob_hi[index])
            {
                near_wall = 1;
            }
        }
        else
        {
            if ((part.pos(index) - 0.5*p3m_rad) < prob_lo[index])
            {
                near_wall = 1;
            }
        }
    }
    else
    {
        amrex::Abort("Incorrect function call!!!! \n");
    }
    
    return near_wall;
}

/*
  calculate location of an image charge reflected across some plane--assumes boundary
  is parallel to either xy,xz,yz plane
  inputs: particle p             : just need the particles position 
          int dir                : either +-1,+-2,+-3
          real(3) im_charge_pos  : output, 0=part's p3m_radius does not overlap with wall, 1=yes   

          abs(dir) > 0 means reflect part pos across plane of x(dir) = prob_hi(dir)
                   < 0 means reflect part pos across plane of x(dir) = prob_lo(dir)
*/
amrex::RealVect calc_im_charge_loc (const FhdParticleContainer::ParticleType& part, int dir)
{
    using namespace amrex;
    using common::prob_lo;
    using common::prob_hi;
    
    const int abs_dir = std::abs(dir);
    
    RealVect im_charge_pos;

    if (abs_dir == 1)
    {
        if (dir > 0)
        {
            // place image charge to the right
            im_charge_pos[1] = part.pos(1);
            im_charge_pos[2] = part.pos(2);
            im_charge_pos[0] = part.pos(0) + 2.0*(prob_hi[0]-part.pos(0));
        }
        else
        {
            // place image charge to the left
            im_charge_pos[1] = part.pos(1);
            im_charge_pos[2] = part.pos(2);
            im_charge_pos[0] = part.pos(0) - 2.0*(part.pos(0)-prob_lo[0]);
        }
    }
    else if (abs_dir == 2)
    {
        if (dir > 0)
        {
            // place image charge above
            im_charge_pos[0] = part.pos(0);
            im_charge_pos[2] = part.pos(2);
            im_charge_pos[1] = part.pos(1) + 2.0*(prob_hi[1]-part.pos(1));
        }
       else
       {
           // place image charge below
           im_charge_pos[0] = part.pos(0);
           im_charge_pos[2] = part.pos(2);
           im_charge_pos[1] = part.pos(1) - 2.0*(part.pos(1)-prob_lo[1]);
       }
    }
    else if (abs_dir == 3)
    {
        if (dir > 0)
        {
            // place image charge behind
            im_charge_pos[0] = part.pos(0);
            im_charge_pos[1] = part.pos(1);
            im_charge_pos[2] = part.pos(2) + 2.0*(prob_hi[2]-part.pos(2));
        }
        else
        {
            // place image charge in front
            im_charge_pos[0] = part.pos(0);
            im_charge_pos[1] = part.pos(1);
            im_charge_pos[2] = part.pos(2) - 2.0*(part.pos(2)-prob_lo[2]);
        }
    }
    else
    {
        amrex::Abort("Incorrect function call!!!! \n");
    }
    
    return im_charge_pos;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
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

    Gpu::DeviceScalar<Real> rcount_d(rcount);
    Real* prcount_d = rcount_d.dataPtr();
    
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
                    Gpu::Atomic::Add(prcount_d, 1.0);
                    repulsive_force(p1, p2, dx, dy, dz, r2);
                }
            }
        }
    });

    rcount = std::min(rcount, rcount_d.dataValue());
}

#endif
