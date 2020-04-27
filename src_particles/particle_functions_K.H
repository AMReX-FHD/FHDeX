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

    int near_wall = 0;
    
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
amrex::Real compute_p3m_force_mag (amrex::Real r, amrex::Real dx, amrex::Real dy, amrex::Real dz)
{
    using namespace amrex;
    using common::pkernel_es;

    // this is force per dx**2
    constexpr Real vals6[] =
        {0., 0.0130335, 0.0259347, 0.0384891, 0.0507132, 0.0624416, 0.0735093,
         0.0837511, 0.0931669, 0.101922, 0.109521, 0.116293, 0.121745,
         0.126535, 0.130335, 0.132978, 0.134134, 0.135621, 0.13529, 0.134629,
         0.133308, 0.130995, 0.129343, 0.126205, 0.122075, 0.118276, 0.115137,
         0.110677, 0.106547, 0.102748, 0.0979574, 0.0936624, 0.0898631,
         0.0860637, 0.0816036, 0.0782998, 0.0745005, 0.0707011, 0.0680581,
         0.0645891, 0.0616157, 0.0586423, 0.0561644, 0.0535214, 0.0513739,
         0.0487309, 0.0469138, 0.0450967, 0.0431145, 0.0414626, 0.0399759,
         0.038324, 0.0366721, 0.0355157, 0.034029, 0.0328727, 0.0317164,
         0.0307252, 0.0295689, 0.0285778, 0.0275866, 0.0267607, 0.0259347,
         0.0251088, 0.0242829, 0.0236221, 0.0229613, 0.0221354, 0.0214746,
         0.0209791};

    // these are in fractions of cell size dx
    constexpr Real points6[] =
        {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2, 1.3,
         1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2., 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7,
         2.8, 2.9, 3., 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4., 4.1,
         4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9, 5., 5.1, 5.2, 5.3, 5.4, 5.5,
         5.6, 5.7, 5.8, 5.9, 6., 6.1, 6.2, 6.3, 6.4, 6.5, 6.6, 6.7, 6.8, 6.9};

    // this is force per dx**2
    constexpr Real vals4[] =
        {0., 0.0191993, 0.0384104, 0.0561438, 0.0745141, 0.0914394, 0.106869,
         0.121719, 0.134211, 0.145262, 0.154958, 0.162921, 0.168841, 0.173264,
         0.176417, 0.17854, 0.177695, 0.175969, 0.174664, 0.17083, 0.165755,
         0.162663, 0.156037, 0.148626, 0.144648, 0.1381, 0.131473, 0.12615,
         0.120057, 0.112581, 0.106808, 0.100863, 0.0947845, 0.0905315, 0.0856055,
         0.0818209, 0.0768459, 0.0728331, 0.0694836, 0.0650267, 0.0627171, 0.0591556,
         0.0561682, 0.0539003, 0.0511791, 0.0492405, 0.0470903, 0.0453286, 0.043644, 0.0417381};

    // these are in fractions of cell size dx
    constexpr Real points4[] =
        {0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3,
         1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7,
         2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0, 4.1,
         4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9};

    int r_cell;
    Real r_norm, r_cell_frac, m;
    Real mag = 0.0;

    // Divison/multiplication by 0.1 below is bc the spacing in the table above is 0.1
    r_norm = r/dx;                   // separation dist in units of dx=dy=dz
    r_cell = std::floor(r_norm/0.1); // scaling by 10 allows r_cell to index the points/val arrays

    r_cell_frac = r_norm/0.1-r_cell;
    r_cell_frac = r_cell_frac*0.1;   // for use in point-slope formula below

    if (pkernel_es == 6)
    {
        // do linear interpolation of force between vals(i+1) and val(i) 
        m = (vals6[r_cell+1]-vals6[r_cell])/(points6[r_cell+1]-points6[r_cell]);
        
        mag  = m*r_cell_frac + vals6[r_cell];
    }
    else if (pkernel_es == 4)
    {
        // do linear interpolation of force between vals(i+1) and val(i) 
        m = (vals4[r_cell+1]-vals4[r_cell])/(points4[r_cell+1]-points4[r_cell]);
        
        mag  = m*r_cell_frac + vals4[r_cell];
    }
    else
    {
        amrex::Abort("P3M implemented only for pkernel 4 and 6! \n");
    }
    
    return mag;
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

template <typename T>
amrex::Gpu::DeviceVector<T> getDeviceVector (const amrex::Vector<T>& vec)
{
    amrex::Gpu::DeviceVector<T> dvec(vec.size());
    for (std::size_t i = 0; i < vec.size(); ++i)
    {
        dvec[i] = vec[i];
    }
    return dvec;
}

void compute_p3m_sr_correction_nl (FhdParticleContainer::AoS& aos, int Np, int Nn,
                                   amrex::NeighborList<FhdParticleContainer::ParticleType>& neighbor_list,
                                   const amrex::Real* dxp, amrex::Real& rcount)
{
    using namespace amrex;
    
    GpuArray<Real, 3> dx = {dxp[0], dxp[1], dxp[2]};
    
    Gpu::DeviceScalar<Real> rcount_d(rcount);
    Real* prcount_d = rcount_d.dataPtr();
    
    auto nbor_data = neighbor_list.data();
    FhdParticleContainer::ParticleType* pstruct = aos().dataPtr();

    Real ee = 1.0/(common::permittivity*4*3.142);
    Real dx2_inv = 1.0/(dx[0]*dx[0]); // assumes isotropic cells

    // clearly we don't want to do this every time...
    Gpu::DeviceVector<int> bc_es_lo_d = getDeviceVector(bc_es_lo);
    Gpu::DeviceVector<int> bc_es_hi_d = getDeviceVector(bc_es_hi);    

    auto p_bc_es_lo_d = bc_es_lo_d.data();
    auto p_bc_es_hi_d = bc_es_hi_d.data();
    
    for (int i = 0; i < Np; ++i) 
    {
        FhdParticleContainer::ParticleType& p1 = pstruct[i];

        Real r = 0.0;
        Real r2 = 0.0;
        RealVect dr(0., 0., 0.);

        Real correction_force_mag = 0.0;
        Real p3m_radius = p1.rdata(FHD_realData::p3m_radius);
        Real q = p1.rdata(FHD_realData::q);

        int near_wall_below = 0;
        if (p_bc_es_lo_d[1] != -1) near_wall_below = near_wall_check(p1, -2);

        int near_wall_above = 0;
        if (p_bc_es_hi_d[1] != -1) near_wall_above = near_wall_check(p1,  2);

        // image charge interactions for wall below
        if (near_wall_below)
        {
            RealVect im_charge_pos = calc_im_charge_loc(p1, -2);

            // compute sep vector btwn charge and its image
            dr[0] = p1.pos(0) - im_charge_pos[0];
            dr[1] = p1.pos(1) - im_charge_pos[1];
            dr[2] = p1.pos(2) - im_charge_pos[2];

            r2 = dr.dotProduct(dr);
            r = std::sqrt(r2);      // separation dist

            // perform coulomb and p3m interaction with image charge
            if ((p_bc_es_lo_d[1] == 1) and (r < p3m_radius))
            {
                // hom. dirichlet--image charge opposite that of particle
                
                // coloumb
                p1.rdata(FHD_realData::forcex) += ee*(dr[0]/r)*q*(-q)/r2;
                p1.rdata(FHD_realData::forcey) += ee*(dr[1]/r)*q*(-q)/r2;
                p1.rdata(FHD_realData::forcez) += ee*(dr[2]/r)*q*(-q)/r2;

                // p3m
                Real correction_force_mag = compute_p3m_force_mag(r, dx[0], dx[1], dx[2]);
                p1.rdata(FHD_realData::forcex) -= ee*(dr[0]/r)*q*(-q)*correction_force_mag*dx2_inv;
                p1.rdata(FHD_realData::forcey) -= ee*(dr[1]/r)*q*(-q)*correction_force_mag*dx2_inv;
                p1.rdata(FHD_realData::forcez) -= ee*(dr[2]/r)*q*(-q)*correction_force_mag*dx2_inv;
            }
            else if ((p_bc_es_lo_d[1] == 2) and (r < p3m_radius))
            {
                // hom. neumann  --image charge equal that of particle

                // coloumb
                p1.rdata(FHD_realData::forcex) += ee*(dr[0]/r)*q*q/r2;
                p1.rdata(FHD_realData::forcey) += ee*(dr[1]/r)*q*q/r2;
                p1.rdata(FHD_realData::forcez) += ee*(dr[2]/r)*q*q/r2;

                // p3m
                Real correction_force_mag = compute_p3m_force_mag(r, dx[0], dx[1], dx[2]);
                p1.rdata(FHD_realData::forcex) -= ee*(dr[0]/r)*q*q*correction_force_mag*dx2_inv;
                p1.rdata(FHD_realData::forcey) -= ee*(dr[1]/r)*q*q*correction_force_mag*dx2_inv;
                p1.rdata(FHD_realData::forcez) -= ee*(dr[2]/r)*q*q*correction_force_mag*dx2_inv;
            }
        }

        // image charge interactions for wall above
        if (near_wall_above)
        {
            RealVect im_charge_pos = calc_im_charge_loc(p1, 2);

            // compute sep vector btwn charge and its image
            dr[0] = p1.pos(0) - im_charge_pos[0];
            dr[1] = p1.pos(1) - im_charge_pos[1];
            dr[2] = p1.pos(2) - im_charge_pos[2];

            r2 = dr.dotProduct(dr);
            r = std::sqrt(r2);      // separation dist

            // perform coulomb and p3m interaction with image charge
            if ((p_bc_es_hi_d[1] == 1) and (r < p3m_radius))
            {
                // hom. dirichlet--image charge opposite that of particle
                
                // coloumb
                p1.rdata(FHD_realData::forcex) += ee*(dr[0]/r)*q*(-q)/r2;
                p1.rdata(FHD_realData::forcey) += ee*(dr[1]/r)*q*(-q)/r2;
                p1.rdata(FHD_realData::forcez) += ee*(dr[2]/r)*q*(-q)/r2;

                // p3m
                Real correction_force_mag = compute_p3m_force_mag(r, dx[0], dx[1], dx[2]);
                p1.rdata(FHD_realData::forcex) -= ee*(dr[0]/r)*q*(-q)*correction_force_mag*dx2_inv;
                p1.rdata(FHD_realData::forcey) -= ee*(dr[1]/r)*q*(-q)*correction_force_mag*dx2_inv;
                p1.rdata(FHD_realData::forcez) -= ee*(dr[2]/r)*q*(-q)*correction_force_mag*dx2_inv;
            }
            else if ((p_bc_es_hi_d[1] == 2) and (r < p3m_radius))
            {
                // hom. neumann  --image charge equal that of particle

                // coloumb
                p1.rdata(FHD_realData::forcex) += ee*(dr[0]/r)*q*q/r2;
                p1.rdata(FHD_realData::forcey) += ee*(dr[1]/r)*q*q/r2;
                p1.rdata(FHD_realData::forcez) += ee*(dr[2]/r)*q*q/r2;

                // p3m
                Real correction_force_mag = compute_p3m_force_mag(r, dx[0], dx[1], dx[2]);
                p1.rdata(FHD_realData::forcex) -= ee*(dr[0]/r)*q*q*correction_force_mag*dx2_inv;
                p1.rdata(FHD_realData::forcey) -= ee*(dr[1]/r)*q*q*correction_force_mag*dx2_inv;
                p1.rdata(FHD_realData::forcez) -= ee*(dr[2]/r)*q*q*correction_force_mag*dx2_inv;
            }
        }

        // loop through neighbor list
        for (const auto& p2 : nbor_data.getNeighbors(i))
        {
            Real q2 = p2.rdata(FHD_realData::q);

            dr[0] = p1.pos(0) - p2.pos(0);
            dr[1] = p1.pos(1) - p2.pos(1);
            dr[2] = p1.pos(2) - p2.pos(2);

            r2 = dr.dotProduct(dr);
            r = std::sqrt(r2);      // separation dist

            if (r < p3m_radius)
            {
                // coloumb
                p1.rdata(FHD_realData::forcex) += ee*(dr[0]/r)*q*q2/r2;
                p1.rdata(FHD_realData::forcey) += ee*(dr[1]/r)*q*q2/r2;
                p1.rdata(FHD_realData::forcez) += ee*(dr[2]/r)*q*q2/r2;

                //p3m
                Real correction_force_mag = compute_p3m_force_mag(r, dx[0], dx[1], dx[2]);
                p1.rdata(FHD_realData::forcex) -= ee*(dr[0]/r)*q*q2*correction_force_mag*dx2_inv;
                p1.rdata(FHD_realData::forcey) -= ee*(dr[1]/r)*q*q2*correction_force_mag*dx2_inv;
                p1.rdata(FHD_realData::forcez) -= ee*(dr[2]/r)*q*q2*correction_force_mag*dx2_inv;

                // If part near wall, check if NL part also near wall
                if (near_wall_below)
                {
                    int near_wall_below_nl_part = near_wall_check(p2, -2);
                    if (near_wall_below_nl_part)
                    {
                        RealVect im_charge_pos = calc_im_charge_loc(p2, -2);

                        // compute sep vector btwn charge and NL's image
                        dr[0] = p1.pos(0) - im_charge_pos[0];
                        dr[1] = p1.pos(1) - im_charge_pos[1];
                        dr[2] = p1.pos(2) - im_charge_pos[2];

                        r2 = dr.dotProduct(dr);
                        r = std::sqrt(r2);      // separation dist

                        // perform coulomb and p3m interaction with NL particle's image charge
                        if ((p_bc_es_lo_d[1] == 1) and (r < p3m_radius))
                        {
                            // hom. dirichlet--image charge opposite that of particle
                
                            // coloumb
                            p1.rdata(FHD_realData::forcex) += ee*(dr[0]/r)*q*(-q2)/r2;
                            p1.rdata(FHD_realData::forcey) += ee*(dr[1]/r)*q*(-q2)/r2;
                            p1.rdata(FHD_realData::forcez) += ee*(dr[2]/r)*q*(-q2)/r2;

                            // p3m
                            Real correction_force_mag = compute_p3m_force_mag(r, dx[0], dx[1], dx[2]);
                            p1.rdata(FHD_realData::forcex) -= ee*(dr[0]/r)*q*(-q2)*correction_force_mag*dx2_inv;
                            p1.rdata(FHD_realData::forcey) -= ee*(dr[1]/r)*q*(-q2)*correction_force_mag*dx2_inv;
                            p1.rdata(FHD_realData::forcez) -= ee*(dr[2]/r)*q*(-q2)*correction_force_mag*dx2_inv;
                        }
                        else if ((p_bc_es_lo_d[1] == 2) and (r < p3m_radius))
                        {
                            // hom. neumann  --image charge equal that of particle
                            
                            // coloumb
                            p1.rdata(FHD_realData::forcex) += ee*(dr[0]/r)*q*q2/r2;
                            p1.rdata(FHD_realData::forcey) += ee*(dr[1]/r)*q*q2/r2;
                            p1.rdata(FHD_realData::forcez) += ee*(dr[2]/r)*q*q2/r2;
                            
                            // p3m
                            Real correction_force_mag = compute_p3m_force_mag(r, dx[0], dx[1], dx[2]);
                            p1.rdata(FHD_realData::forcex) -= ee*(dr[0]/r)*q*q2*correction_force_mag*dx2_inv;
                            p1.rdata(FHD_realData::forcey) -= ee*(dr[1]/r)*q*q2*correction_force_mag*dx2_inv;
                            p1.rdata(FHD_realData::forcez) -= ee*(dr[2]/r)*q*q2*correction_force_mag*dx2_inv;
                        }
                    }

                    int near_wall_above_nl_part = near_wall_check(p2, 2);
                    if (near_wall_above_nl_part)
                    {
                        RealVect im_charge_pos = calc_im_charge_loc(p2, 2);

                        // compute sep vector btwn charge and NL's image
                        dr[0] = p1.pos(0) - im_charge_pos[0];
                        dr[1] = p1.pos(1) - im_charge_pos[1];
                        dr[2] = p1.pos(2) - im_charge_pos[2];

                        r2 = dr.dotProduct(dr);
                        r = std::sqrt(r2);      // separation dist
                    
                        // perform coulomb and p3m interaction with image charge
                        if ((p_bc_es_hi_d[1] == 1) and (r < p3m_radius))
                        {
                            // hom. dirichlet--image charge opposite that of particle
                            
                            // coloumb
                            p1.rdata(FHD_realData::forcex) += ee*(dr[0]/r)*q*(-q2)/r2;
                            p1.rdata(FHD_realData::forcey) += ee*(dr[1]/r)*q*(-q2)/r2;
                            p1.rdata(FHD_realData::forcez) += ee*(dr[2]/r)*q*(-q2)/r2;

                            // p3m
                            Real correction_force_mag = compute_p3m_force_mag(r, dx[0], dx[1], dx[2]);
                            p1.rdata(FHD_realData::forcex) -= ee*(dr[0]/r)*q*(-q2)*correction_force_mag*dx2_inv;
                            p1.rdata(FHD_realData::forcey) -= ee*(dr[1]/r)*q*(-q2)*correction_force_mag*dx2_inv;
                            p1.rdata(FHD_realData::forcez) -= ee*(dr[2]/r)*q*(-q2)*correction_force_mag*dx2_inv;
                        }
                        else if ((p_bc_es_hi_d[1] == 2) and (r < p3m_radius))
                        {
                            // hom. neumann  --image charge equal that of particle

                            // coloumb
                            p1.rdata(FHD_realData::forcex) += ee*(dr[0]/r)*q*q2/r2;
                            p1.rdata(FHD_realData::forcey) += ee*(dr[1]/r)*q*q2/r2;
                            p1.rdata(FHD_realData::forcez) += ee*(dr[2]/r)*q*q2/r2;

                            // p3m
                            Real correction_force_mag = compute_p3m_force_mag(r, dx[0], dx[1], dx[2]);
                            p1.rdata(FHD_realData::forcex) -= ee*(dr[0]/r)*q*2*correction_force_mag*dx2_inv;
                            p1.rdata(FHD_realData::forcey) -= ee*(dr[1]/r)*q*2*correction_force_mag*dx2_inv;
                            p1.rdata(FHD_realData::forcez) -= ee*(dr[2]/r)*q*2*correction_force_mag*dx2_inv;
                        }
                    }
                }                
            }
        }
    }

    rcount = std::min(rcount, rcount_d.dataValue());    
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
    Gpu::DeviceVector<int> p_int_tog_d = getDeviceVector(p_int_tog);
    auto pp_int_tog_d = p_int_tog_d.data();

    auto rmin = common::rmin;
    
    AMREX_FOR_1D( Np, i, 
    {
        FhdParticleContainer::ParticleType& p1 = pstruct[i];

        for (const auto& p2 : nbor_data.getNeighbors(i))
        {       
            if ( (pp_int_tog_d[p1.idata(FHD_intData::species)] != 0) and 
                 (pp_int_tog_d[p2.idata(FHD_intData::species)] != 0) )
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
