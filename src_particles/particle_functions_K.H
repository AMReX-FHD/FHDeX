#ifndef _particle_functions_K_H_
#define _particle_functions_K_H_

#include <AMReX.H>
#include <FhdParticleContainer.H>
#include <common_namespace.H>
#include <kernel_functions_K.H>
#include <math.h>

/**
   Check if particle is near a boundary. Assumes boundary is parallel to either xy, zx, yz, plane
   
   Inputs : Particle part
            int dir: either +-1,+-2,+-3

   Returns: 
            int near_wall: 0=part's p3m_radius/2 does not overlap with wall, 1=yes  
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int near_wall_check (const FhdParticleContainer::ParticleType& part, int dir, Real dist)
{
    using namespace amrex;
    using common::prob_hi;
    using common::prob_lo;
    
    const int abs_dir = std::abs(dir);
    //const Real p3m_rad = part.rdata(FHD_realData::p3m_radius);
    const int index = abs_dir - 1;

    int near_wall = 0;
    
    if ((abs_dir == 1) or (abs_dir == 2) or (abs_dir == 3))
    {
        if (dir > 0)
        {
            if ((part.pos(index) + dist) > prob_hi[index])
            {
                near_wall = 1;
            }
        }
        else
        {
            if ((part.pos(index) - dist) < prob_lo[index])
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
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::RealVect calc_im_charge_loc (const FhdParticleContainer::ParticleType& part, int dir) noexcept
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

    //Print() << "RCELLFRAC: " << r_cell_frac << ", " << r_norm << "\n";

    if (pkernel_es[0] == 6)
    {
        // do linear interpolation of force between vals(i+1) and val(i) 
        m = (vals6[r_cell+1]-vals6[r_cell])/(points6[r_cell+1]-points6[r_cell]);
        
        mag  = m*r_cell_frac + vals6[r_cell];
    }
    else if (pkernel_es[0] == 4)
    {
        // do linear interpolation of force between vals(i+1) and val(i) 
        m = (vals4[r_cell+1]-vals4[r_cell])/(points4[r_cell+1]-points4[r_cell]);
        
        mag  = m*r_cell_frac + vals4[r_cell];

        //Print()  << scientific << setprecision(15) << "MAGGPU: " << m << ", " << r_cell_frac << ", " << vals4[r_cell] << ", " << r_cell << "\n";
    }
    else
    {
        amrex::Abort("P3M implemented only for pkernel 4 and 6! \n");
    }
    
    return mag;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void force_function2 (      FhdParticleContainer::ParticleType& part1,
                      const FhdParticleContainer::ParticleType& part2,
                      const amrex::GpuArray<amrex::Real, 3>& domsize)
{
    using namespace amrex;
    using common::permittivity;
    using common::images;
    
    Real ee = (1.0/(permittivity*4*3.14159265));

    Real dx0[3];
    dx0[0] = part1.pos(0) - part2.pos(0);
    dx0[1] = part1.pos(1) - part2.pos(1);
    dx0[2] = part1.pos(2) - part2.pos(2);

    int imag = (images == 0) ? 1 : images;
    
    Real maxdist = 0.99*amrex::min(imag * domsize[0],
                                   imag * domsize[1],
                                   imag * domsize[2]);

    Real q1 = part1.rdata(FHD_realData::q);
    Real q2 = part2.rdata(FHD_realData::q);
    
    for (int ii = -images; ii <= images; ++ii) {
        for (int jj = -images; jj <= images; ++jj) {
            for (int kk = -images; kk <= images; ++kk) {

                Real dx[3];
                dx[0] = dx0[0] - ii*domsize[0];
                dx[1] = dx0[1] - jj*domsize[1];
                dx[2] = dx0[2] - kk*domsize[2];

                Real dr2 = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
                if (dr2 != 0) //might find a more careful way of doing this
                {
                    Real rtdr2 = std::sqrt(dr2);
                    if (rtdr2 < maxdist)
                    {
                        part1.rdata(FHD_realData::forcex) += ee*(dx[0]/rtdr2)*q1*q2/dr2;
                        part1.rdata(FHD_realData::forcey) += ee*(dx[1]/rtdr2)*q1*q2/dr2;
                        part1.rdata(FHD_realData::forcez) += ee*(dx[2]/rtdr2)*q1*q2/dr2;
                    }
                }
            }
        }
    }
}



AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void repulsive_force_lj (      FhdParticleContainer::ParticleType& part1,
                      amrex::Real dx, amrex::Real dy, amrex::Real dz, amrex::Real dr, amrex::Real sig, amrex::Real eps)
{
//    using namespace amrex;

    const Real dr2 = dr*dr;
    const Real fac1 = std::pow(sig, 12);
    const Real fac2 = std::pow(sig, 6);
    const Real ff = eps*4*(1./(dr2*dr2*dr2*dr2))*(-12*fac1/(dr2*dr2*dr2) + 6*fac2);
    
    part1.rdata(FHD_realData::forcex) -= dx*ff;
    part1.rdata(FHD_realData::forcey) -= dy*ff;
    part1.rdata(FHD_realData::forcez) -= dz*ff;
}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void repulsive_force_wall (      FhdParticleContainer::ParticleType& part1,
                      amrex::Real dx, amrex::Real dy, amrex::Real dz, amrex::Real dr, amrex::Real sig, amrex::Real eps)
{
    using namespace amrex;

    //note that we are adjusting eps for the different wall potential

    const Real dr2 = dr*dr;
    const Real fac1 = std::pow((sig), 9);
    const Real fac2 = std::pow((sig), 3);
    const Real ff = eps*1.5*sqrt(3)*(1./(dr2*dr2*dr))*(-9*fac1/(dr2*dr2*dr2) + 3*fac2);
    
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

void compute_p3m_sr_correction_nl_gpu (FhdParticleContainer::AoS& aos, int Np, int Nn,
                                   amrex::NeighborList<FhdParticleContainer::ParticleType>& neighbor_list,
                                   const amrex::Real* dxp, amrex::Real& rcount, amrex::Real& rcountI)
{
    using namespace amrex;
    
    GpuArray<Real, 3> dx = {dxp[0], dxp[1], dxp[2]};
    
    Gpu::DeviceScalar<Real> rcount_d(rcount);
    Real* prcount_d = rcount_d.dataPtr();

    Gpu::DeviceScalar<Real> rcount_di(rcountI);
    Real* prcount_di = rcount_di.dataPtr();
    
    auto nbor_data = neighbor_list.data();
    FhdParticleContainer::ParticleType* pstruct = aos().dataPtr();

    Real ee = 1.0/(common::permittivity*4*3.141592653589793238);
    Real dx2_inv = 1.0/(dx[0]*dx[0]); // assumes isotropic cells

    // clearly we don't want to do this every time...
    Gpu::DeviceVector<int> bc_es_lo_d = getDeviceVector(bc_es_lo);
    Gpu::DeviceVector<int> bc_es_hi_d = getDeviceVector(bc_es_hi);    
    //Gpu::DeviceVector<int> p_int_tog_d = getDeviceVector(p_int_tog);

    auto p_bc_es_lo_d = bc_es_lo_d.data();
    auto p_bc_es_hi_d = bc_es_hi_d.data();
    //auto p_p_int_tog_d = p_int_tog_d.data();

    AMREX_FOR_1D( Np, i, 
    {
        FhdParticleContainer::ParticleType& p1 = pstruct[i];
        if(p1.idata(FHD_intData::pinned) == 0)
        {
            int spec1 = p1.idata(FHD_intData::species)-1;
            int spec2 = p1.idata(FHD_intData::species)-1;
            int p_matrix_val = spec1*nspecies + spec2;

            Real r = 0.0;
            Real r2 = 0.0;
            RealVect dr(0., 0., 0.);

            Real correction_force_mag = 0.0;
            Real p3m_radius = p1.rdata(FHD_realData::p3m_radius);
            Real q = p1.rdata(FHD_realData::q);

            int near_wall_below = 0;
            if (p_bc_es_lo_d[1] != -1) near_wall_below = near_wall_check(p1, -2, 0.5*p3m_radius);

            int near_wall_above = 0;
            if (p_bc_es_hi_d[1] != -1) near_wall_above = near_wall_check(p1,  2, 0.5*p3m_radius);

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
                    if(p_int_tog[p_matrix_val] != 0)
                    {
                        p1.rdata(FHD_realData::forcex) += ee*(dr[0]/r)*q*(-q)/r2;
                        p1.rdata(FHD_realData::forcey) += ee*(dr[1]/r)*q*(-q)/r2;
                        p1.rdata(FHD_realData::forcez) += ee*(dr[2]/r)*q*(-q)/r2;
                    }
                    // p3m
                    Real correction_force_mag = compute_p3m_force_mag(r, dx[0], dx[1], dx[2]);
                    Gpu::Atomic::Add(prcount_di, 1.0);
                    p1.rdata(FHD_realData::forcex) -= ee*(dr[0]/r)*q*(-q)*correction_force_mag*dx2_inv;
                    p1.rdata(FHD_realData::forcey) -= ee*(dr[1]/r)*q*(-q)*correction_force_mag*dx2_inv;
                    p1.rdata(FHD_realData::forcez) -= ee*(dr[2]/r)*q*(-q)*correction_force_mag*dx2_inv;

                    //Print() << "Below!\n";

                    //Print() << p1.id() << ": " << setprecision(15) << correction_force_mag << "\n";
                }
                else if ((p_bc_es_lo_d[1] == 2) and (r < p3m_radius))
                {
                    // hom. neumann  --image charge equal that of particle

                    // coloumb
                    if(p_int_tog[p_matrix_val] != 0)
                    {
                        p1.rdata(FHD_realData::forcex) += ee*(dr[0]/r)*q*q/r2;
                        p1.rdata(FHD_realData::forcey) += ee*(dr[1]/r)*q*q/r2;
                        p1.rdata(FHD_realData::forcez) += ee*(dr[2]/r)*q*q/r2;
                    }

                    // p3m
                    Real correction_force_mag = compute_p3m_force_mag(r, dx[0], dx[1], dx[2]);
                    Gpu::Atomic::Add(prcount_di, 1.0);
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
                    if(p_int_tog[p_matrix_val] != 0)
                    {
                        p1.rdata(FHD_realData::forcex) += ee*(dr[0]/r)*q*(-q)/r2;
                        p1.rdata(FHD_realData::forcey) += ee*(dr[1]/r)*q*(-q)/r2;
                        p1.rdata(FHD_realData::forcez) += ee*(dr[2]/r)*q*(-q)/r2;
                    }
                    // p3m
                    Real correction_force_mag = compute_p3m_force_mag(r, dx[0], dx[1], dx[2]);
                    Gpu::Atomic::Add(prcount_di, 1.0);
                    p1.rdata(FHD_realData::forcex) -= ee*(dr[0]/r)*q*(-q)*correction_force_mag*dx2_inv;
                    p1.rdata(FHD_realData::forcey) -= ee*(dr[1]/r)*q*(-q)*correction_force_mag*dx2_inv;
                    p1.rdata(FHD_realData::forcez) -= ee*(dr[2]/r)*q*(-q)*correction_force_mag*dx2_inv;
                }
                else if ((p_bc_es_hi_d[1] == 2) and (r < p3m_radius))
                {
                    // hom. neumann  --image charge equal that of particle

                    // coloumb
                    if(p_int_tog[p_matrix_val] != 0)
                    {
                        p1.rdata(FHD_realData::forcex) += ee*(dr[0]/r)*q*q/r2;
                        p1.rdata(FHD_realData::forcey) += ee*(dr[1]/r)*q*q/r2;
                        p1.rdata(FHD_realData::forcez) += ee*(dr[2]/r)*q*q/r2;
                    }

                    // p3m
                    Real correction_force_mag = compute_p3m_force_mag(r, dx[0], dx[1], dx[2]);
                    Gpu::Atomic::Add(prcount_di, 1.0);
                    p1.rdata(FHD_realData::forcex) -= ee*(dr[0]/r)*q*q*correction_force_mag*dx2_inv;
                    p1.rdata(FHD_realData::forcey) -= ee*(dr[1]/r)*q*q*correction_force_mag*dx2_inv;
                    p1.rdata(FHD_realData::forcez) -= ee*(dr[2]/r)*q*q*correction_force_mag*dx2_inv;

                }

            }

            // loop through neighbor list
            for (const auto& p2 : nbor_data.getNeighbors(i))
            {
                Real q2 = p2.rdata(FHD_realData::q);

                spec2 = p2.idata(FHD_intData::species)-1;
                p_matrix_val = spec1*nspecies + spec2;

                dr[0] = p1.pos(0) - p2.pos(0);
                dr[1] = p1.pos(1) - p2.pos(1);
                dr[2] = p1.pos(2) - p2.pos(2);

                r2 = dr.dotProduct(dr);
                r = std::sqrt(r2);      // separation dist

                if (r < p3m_radius)
                {
                    // coloumb
                    if(p_int_tog[p_matrix_val] != 0)
                    {
                        p1.rdata(FHD_realData::forcex) += ee*(dr[0]/r)*q*q2/r2;
                        p1.rdata(FHD_realData::forcey) += ee*(dr[1]/r)*q*q2/r2;
                        p1.rdata(FHD_realData::forcez) += ee*(dr[2]/r)*q*q2/r2;
                    }

                    //p3m
                    Real correction_force_mag = compute_p3m_force_mag(r, dx[0], dx[1], dx[2]);
                    Gpu::Atomic::Add(prcount_d, 1.0);
                    p1.rdata(FHD_realData::forcex) -= ee*(dr[0]/r)*q*q2*correction_force_mag*dx2_inv;
                    p1.rdata(FHD_realData::forcey) -= ee*(dr[1]/r)*q*q2*correction_force_mag*dx2_inv;
                    p1.rdata(FHD_realData::forcez) -= ee*(dr[2]/r)*q*q2*correction_force_mag*dx2_inv;

                    //Print() << scientific << setprecision(15) << correction_force_mag << ", " << r << ", " << dx[0] << std::endl;

                    // If part near wall, check if NL part also near wall
                    if (near_wall_below)
                    {
                        int near_wall_below_nl_part = near_wall_check(p2, -2, 0.5*p3m_radius);
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

                                //Print() << "Pre force: " << setprecision(17) << p1.rdata(FHD_realData::forcex) << ", " << p1.rdata(FHD_realData::forcey) << ", " << p1.rdata(FHD_realData::forcez) << "\n";

                                // coloumb
                                if(p_int_tog[p_matrix_val] != 0)
                                {
                                    p1.rdata(FHD_realData::forcex) += ee*(dr[0]/r)*q*(-q2)/r2;
                                    p1.rdata(FHD_realData::forcey) += ee*(dr[1]/r)*q*(-q2)/r2;
                                    p1.rdata(FHD_realData::forcez) += ee*(dr[2]/r)*q*(-q2)/r2;
                                }
                                   // Print() << "Coulomb: " << setprecision(17) << ee*(dr[0]/r)*q*(-q2)/r2 << ", " << ee*(dr[1]/r)*q*(-q2)/r2 << ", " << ee*(dr[2]/r)*q*(-q2)/r2 << "\n";

                                // p3m
                                Real correction_force_mag = compute_p3m_force_mag(r, dx[0], dx[1], dx[2]);
                                Gpu::Atomic::Add(prcount_di, 1.0);
                                p1.rdata(FHD_realData::forcex) -= ee*(dr[0]/r)*q*(-q2)*correction_force_mag*dx2_inv;
                                p1.rdata(FHD_realData::forcey) -= ee*(dr[1]/r)*q*(-q2)*correction_force_mag*dx2_inv;
                                p1.rdata(FHD_realData::forcez) -= ee*(dr[2]/r)*q*(-q2)*correction_force_mag*dx2_inv;

                                //Print() << "Below: " << p1.id() << ": " << setprecision(17) << correction_force_mag << ", " << p1.rdata(FHD_realData::forcex) << ", " << p1.rdata(FHD_realData::forcey) << ", " << p1.rdata(FHD_realData::forcez) << "\n";
                            }
                            else if ((p_bc_es_lo_d[1] == 2) and (r < p3m_radius))
                            {
                                // hom. neumann  --image charge equal that of particle
                                
                                // coloumb
                                if(p_int_tog[p_matrix_val] != 0)
                                {
                                    p1.rdata(FHD_realData::forcex) += ee*(dr[0]/r)*q*q2/r2;
                                    p1.rdata(FHD_realData::forcey) += ee*(dr[1]/r)*q*q2/r2;
                                    p1.rdata(FHD_realData::forcez) += ee*(dr[2]/r)*q*q2/r2;
                                }
                                
                                // p3m
                                Real correction_force_mag = compute_p3m_force_mag(r, dx[0], dx[1], dx[2]);
                                Gpu::Atomic::Add(prcount_di, 1.0);
                                p1.rdata(FHD_realData::forcex) -= ee*(dr[0]/r)*q*q2*correction_force_mag*dx2_inv;
                                p1.rdata(FHD_realData::forcey) -= ee*(dr[1]/r)*q*q2*correction_force_mag*dx2_inv;
                                p1.rdata(FHD_realData::forcez) -= ee*(dr[2]/r)*q*q2*correction_force_mag*dx2_inv;
                            }
                        }
                    }

                   if (near_wall_above)
                   {
                        int near_wall_above_nl_part = near_wall_check(p2, 2, 0.5*p3m_radius);
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
                                if(p_int_tog[p_matrix_val] != 0)
                                {
                                    p1.rdata(FHD_realData::forcex) += ee*(dr[0]/r)*q*(-q2)/r2;
                                    p1.rdata(FHD_realData::forcey) += ee*(dr[1]/r)*q*(-q2)/r2;
                                    p1.rdata(FHD_realData::forcez) += ee*(dr[2]/r)*q*(-q2)/r2;
                                }
                                // p3m
                                Real correction_force_mag = compute_p3m_force_mag(r, dx[0], dx[1], dx[2]);
                                Gpu::Atomic::Add(prcount_di, 1.0);
                                p1.rdata(FHD_realData::forcex) -= ee*(dr[0]/r)*q*(-q2)*correction_force_mag*dx2_inv;
                                p1.rdata(FHD_realData::forcey) -= ee*(dr[1]/r)*q*(-q2)*correction_force_mag*dx2_inv;
                                p1.rdata(FHD_realData::forcez) -= ee*(dr[2]/r)*q*(-q2)*correction_force_mag*dx2_inv;

                               // Print() << "Above: " << p1.id() << ": " << setprecision(15) << correction_force_mag << "\n";
                            }
                            else if ((p_bc_es_hi_d[1] == 2) and (r < p3m_radius))
                            {
                                // hom. neumann  --image charge equal that of particle

                                // coloumb
                                if(p_int_tog[p_matrix_val] != 0)
                                {
                                    p1.rdata(FHD_realData::forcex) += ee*(dr[0]/r)*q*q2/r2;
                                    p1.rdata(FHD_realData::forcey) += ee*(dr[1]/r)*q*q2/r2;
                                    p1.rdata(FHD_realData::forcez) += ee*(dr[2]/r)*q*q2/r2;
                                }

                                // p3m
                                Real correction_force_mag = compute_p3m_force_mag(r, dx[0], dx[1], dx[2]);
                                Gpu::Atomic::Add(prcount_di, 1.0);
                                p1.rdata(FHD_realData::forcex) -= ee*(dr[0]/r)*q*q2*correction_force_mag*dx2_inv;
                                p1.rdata(FHD_realData::forcey) -= ee*(dr[1]/r)*q*q2*correction_force_mag*dx2_inv;
                                p1.rdata(FHD_realData::forcez) -= ee*(dr[2]/r)*q*q2*correction_force_mag*dx2_inv;


                            }
                        }
                    }                
                }


            }
        }
    });

    rcount = rcount_d.dataValue();
    rcountI = rcount_di.dataValue();
}

void compute_forces_nl_gpu (FhdParticleContainer::AoS& aos, int Np, int Nn,
                        amrex::NeighborList<FhdParticleContainer::ParticleType>& neighbor_list,
                        amrex::Real& rcount, amrex::Real& rdcount)
{
    using namespace amrex;

    Gpu::DeviceScalar<Real> rcount_d(rcount);
    Real* prcount_d = rcount_d.dataPtr();
    Gpu::DeviceScalar<Real> rdcount_d(rdcount);
    Real* prdcount_d = rdcount_d.dataPtr();
    
    auto nbor_data = neighbor_list.data();
    FhdParticleContainer::ParticleType* pstruct = aos().dataPtr();
    
    // clearly we don't want to do this every time...
    //Gpu::DeviceVector<int> p_int_tog_d = getDeviceVector(p_int_tog);
   // auto pp_int_tog_d = p_int_tog_d.data();

    auto rmin = common::rmin;

    // clearly we don't want to do this every time...
    Gpu::DeviceVector<int> bc_vel_lo_d = getDeviceVector(bc_vel_lo);
    Gpu::DeviceVector<int> bc_vel_hi_d = getDeviceVector(bc_vel_hi);    
    auto p_bc_vel_lo_d = bc_vel_lo_d.data();
    auto p_bc_vel_hi_d = bc_vel_hi_d.data();

//    Gpu::DeviceVector<Real> prob_lo_d = getDeviceVector(prob_lo);
//    Gpu::DeviceVector<Real> prob_hi_d = getDeviceVector(prob_hi);    
//    auto p_prob_lo_d = prob_lo_d.data();
//    auto p_prob_hi_d = prob_hi_d.data();


    
    AMREX_FOR_1D( Np, i, 
    {
        FhdParticleContainer::ParticleType& p1 = pstruct[i];

        if(p1.idata(FHD_intData::pinned) == 0)
        {
            int spec1 = p1.idata(FHD_intData::species)-1;

            for (const auto& p2 : nbor_data.getNeighbors(i))
            {       

                int spec2 = p2.idata(FHD_intData::species)-1;
                int p_matrix_val = spec1*nspecies + spec2;

                if ( p_int_tog[p_matrix_val] != 0)
                {
                    const Real dx = p1.pos(0) - p2.pos(0);
                    const Real dy = p1.pos(1) - p2.pos(1);
                    const Real dz = p1.pos(2) - p2.pos(2);

                    const Real r = sqrt(dx*dx + dy*dy + dz*dz);

                    const Real rMax = rmax[p_matrix_val]*sigma[p_matrix_val];       

                    if (r <= rMax)
                    {
                        const Real rMinMax = amrex::max(r, rmin[p_matrix_val]*sigma[p_matrix_val]);

                        repulsive_force_lj(p1, dx, dy, dz, rMinMax, sigma[p_matrix_val], eepsilon[p_matrix_val]);
                        Gpu::Atomic::Add(prcount_d, 1.0);
                    }
                }    
            }    
        }
    });

    if (sr_tog > 1)
    {
        GpuArray<Real, 3> plo = {prob_lo[0], prob_lo[1], prob_lo[2]};
        GpuArray<Real, 3> phi = {prob_hi[0], prob_hi[1], prob_hi[2]};

        if(bc_vel_lo[1] != -1)
        {
            AMREX_FOR_1D( Np, i, 
            {
                FhdParticleContainer::ParticleType& p1 = pstruct[i];
                int spec1 = p1.idata(FHD_intData::species)-1;

                if(p1.idata(FHD_intData::pinned) == 0)
                {

                    const Real dx = 0;
                    const Real dy = p1.pos(1) - plo[1];
                    const Real dz = 0;

                    const Real r = sqrt(dx*dx + dy*dy + dz*dz);

                    const Real rMax = rmax_wall[spec1]*sigma_wall[spec1];
                    
                    if (r <= rMax)
                    {
                        Gpu::Atomic::Add(prdcount_d, 1.0);

                        const Real rMinMax = amrex::max(r, rmin_wall[spec1]*sigma_wall[spec1]);
                        if(sr_tog == 2)
                        {
                            repulsive_force_wall(p1, dx, dy, dz, rMinMax, sigma_wall[spec1], eepsilon_wall[spec1]);
                        }else
                        {
                            repulsive_force_lj(p1, dx, dy, dz, rMinMax, sigma_wall[spec1], eepsilon_wall[spec1]);
                        }
                    }
                }
            });
        }
        

        if(bc_vel_hi[1] != -1)
        {
            AMREX_FOR_1D( Np, i, 
            {
                FhdParticleContainer::ParticleType& p1 = pstruct[i];
                int spec1 = p1.idata(FHD_intData::species)-1;

                if(p1.idata(FHD_intData::pinned) == 0)
                {

                    const Real dx = 0;
                    const Real dy = p1.pos(1) - phi[1];
                    const Real dz = 0;

                    const Real r = sqrt(dx*dx + dy*dy + dz*dz);

                    const Real rMax = rmax_wall[spec1]*sigma_wall[spec1];
                    
                    if (r <= rMax)
                    {
                        Gpu::Atomic::Add(prdcount_d, 1.0);

                        const Real rMinMax = amrex::max(r, rmin_wall[spec1]*sigma_wall[spec1]);
                        if(sr_tog == 2)
                        {
                            repulsive_force_wall(p1, dx, dy, dz, rMinMax, sigma_wall[spec1], eepsilon_wall[spec1]);
                        }else
                        {
                            repulsive_force_lj(p1, dx, dy, dz, rMinMax, sigma_wall[spec1], eepsilon_wall[spec1]);
                        }
                    }
                }
            });
        }


    }

    rcount = rcount_d.dataValue();
    rdcount = rdcount_d.dataValue();
}

template <typename F>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void get_weights (const FhdParticleContainer::ParticleType& p, F f,
                  const IntVect& nodal_flag, 
                  amrex::Real weights[3][2*F::ks], int indices[3][2*F::ks],
                  const amrex::GpuArray<amrex::Real, 3>& plo,
                  const amrex::GpuArray<amrex::Real, 3>& dx,
                  const amrex::GpuArray<amrex::Real, 3>& dxi)
{
    using namespace amrex;

    constexpr int ks = F::ks;

    int fi[3], fn[3];
    for (int idim = 0; idim < 3; ++idim) {
        Real fr = (p.pos(idim) - plo[idim]) * dxi[idim];
        fi[idim] = static_cast<int>(std::floor(fr));
        Real fd = fr - fi[idim];
        fn[idim] = (fd < 0.5) ? -1 : 0;
    }
        
    for (int idim = 0; idim < 3; ++idim) {
        for (int i = 0 ; i < 2*ks; ++i) {
            int index = fi[idim] + i -(ks-1) + (1 - nodal_flag[idim]) * fn[idim];
            Real rr = p.pos(idim) - (index + (1-nodal_flag[idim])*0.5)*dx[idim] + plo[idim];
            weights[idim][i] = f(rr*dxi[idim]);
            indices[idim][i] = index;
        }
    }
}

template <typename F>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void spread_op (amrex::Array4<amrex::Real> const& arr, F f, 
                amrex::Real w[3][2*F::ks], int indices[3][2*F::ks],
                amrex::Real value)
{
    constexpr int ks = F::ks;

    for (int k = 0 ; k < 2*ks; ++k) {
        for (int j = 0 ; j < 2*ks; ++j) {
            for (int i = 0 ; i < 2*ks; ++i) {
                amrex::Gpu::Atomic::Add( &arr(indices[0][i], indices[1][j], indices[2][k]),
                                         w[0][i]*w[1][j]*w[2][k]*value);
            }
        }
    }    
}

template <typename F>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real inter_op (amrex::Array4<amrex::Real> const& arr, F f,
               amrex::Real w[3][2*F::ks], int indices[3][2*F::ks])
{
    constexpr int ks = F::ks;

    Real output = 0;

    for (int k = 0 ; k < 2*ks; ++k) {
        for (int j = 0 ; j < 2*ks; ++j) {
            for (int i = 0 ; i < 2*ks; ++i) {
                int ii = indices[0][i];
                int jj = indices[1][j];
                int kk = indices[2][k];
                output += w[0][i]*w[1][j]*w[2][k]*arr(ii, jj, kk);
            }
        }
    }

    return output; 
}

//template <typename F>
//void emf(FhdParticleContainer::AoS& aos, FArrayBox& Ex, FArrayBox& Ey, FArrayBox& Ez, F f,
//                     const amrex::Real* plo_in, const amrex::Real* dx_in)
//{
//    using namespace amrex;
//    constexpr int ks = F::ks;
//    
//    GpuArray<Real, 3> plo = {plo_in[0], plo_in[1], plo_in[2]};
//    GpuArray<Real, 3> dx  = {dx_in[0],  dx_in[1],  dx_in[2] };    
//    GpuArray<Real, 3> dxi = {1.0/dx_in[0],  1.0/dx_in[1],  1.0/dx_in[2] };

//    Real volinv = 1.0/(dx_in[0]*dx_in[1]*dx_in[2]);
//    Real permittivity = common::permittivity;
//    
//    const auto Np = aos.numParticles();
//    const auto pstruct = aos().dataPtr();

//    auto arrEx = charge.array();
//    const IntVect nodal_flag = charge.box().type();
//    
//    amrex::ParallelFor(Np, [=] AMREX_GPU_HOST_DEVICE (int ip) noexcept
//    {
//        FhdParticleContainer::ParticleType& p = pstruct[ip];

//        Real qm = p.rdata(FHD_realData::q)/permittivity;
//        
//        Real w[3][2*ks];
//        int indices[3][2*ks];

//        get_weights(p, f, nodal_flag, w, indices, plo, dx, dxi);
//        inter_op(arr, f, w, indices, p, 1);
//    });
//}

template <typename F>
void collect_charge_gpu (FhdParticleContainer::AoS& aos, FArrayBox& charge, F f,
                     const amrex::Real* plo_in, const amrex::Real* dx_in)
{
    using namespace amrex;
    constexpr int ks = F::ks;
    
    GpuArray<Real, 3> plo = {plo_in[0], plo_in[1], plo_in[2]};
    GpuArray<Real, 3> dx  = {dx_in[0],  dx_in[1],  dx_in[2] };    
    GpuArray<Real, 3> dxi = {1.0/dx_in[0],  1.0/dx_in[1],  1.0/dx_in[2] };

    Real volinv = 1.0/(dx_in[0]*dx_in[1]*dx_in[2]);
    Real permittivity = common::permittivity;
    
    const auto Np = aos.numParticles();
    const auto pstruct = aos().dataPtr();

    auto arr = charge.array();
    const IntVect nodal_flag = charge.box().type();
    
    amrex::ParallelFor(Np, [=] AMREX_GPU_HOST_DEVICE (int ip) noexcept
    {
        FhdParticleContainer::ParticleType& p = pstruct[ip];

        Real qm = p.rdata(FHD_realData::q)/permittivity;
        
        Real w[3][2*ks];
        int indices[3][2*ks];

        get_weights(p, f, nodal_flag, w, indices, plo, dx, dxi);
        spread_op(arr, f, w, indices, qm*volinv);
    });
}

void collect_charge_gpu (FhdParticleContainer::AoS& aos, FArrayBox& charge,
                     const amrex::Real* plo, const amrex::Real* dx)
{
    if (common::pkernel_es[0] == 3)
    {
        collect_charge_gpu (aos, charge, Kernel3P(), plo, dx);
    }
    else if (common::pkernel_es[0] == 4)
    {
        collect_charge_gpu(aos, charge, Kernel4P(), plo, dx);
    }
    else if (common::pkernel_es[0] == 6)
    {
        collect_charge_gpu(aos, charge, Kernel6P(), plo, dx);
    }
}

template <typename F>
void emf_gpu (FhdParticleContainer::AoS& aos, F f,
                      FArrayBox& Ex, FArrayBox& Ey, FArrayBox& Ez,
                      const amrex::Real* plo_in, const amrex::Real* dx_in)
{
    GpuArray<Real, 3> plo = {plo_in[0], plo_in[1], plo_in[2]};
    GpuArray<Real, 3> dx  = {dx_in[0],  dx_in[1],  dx_in[2] };
    GpuArray<Real, 3> dxi = {1.0/dx_in[0],  1.0/dx_in[1],  1.0/dx_in[2] };    

    const auto Np = aos.numParticles();
    const auto pstruct = aos().dataPtr();

    auto Exarr = Ex.array();
    auto Eyarr = Ey.array();
    auto Ezarr = Ez.array();

    const IntVect Exflag = Ex.box().type();
    const IntVect Eyflag = Ey.box().type();
    const IntVect Ezflag = Ez.box().type();

    constexpr int ks = F::ks;
       
    amrex::ParallelFor(Np, [=] AMREX_GPU_HOST_DEVICE (int ip) noexcept
    {
        FhdParticleContainer::ParticleType& p = pstruct[ip];
        //Print() << p.idata(FHD_intData::species) << "\n";
                
        if(p.idata(FHD_intData::pinned) == 0)
        {
            Real w[3][2*ks];
            int indices[3][2*ks];
            //Print() << "Here!\n";
            //E field is cell centered so we can use the same weights/indicies for each component
            get_weights(p, f, Exflag, w, indices, plo, dx, dxi);
            p.rdata(FHD_realData::forcex) += p.rdata(FHD_realData::q)*inter_op(Exarr, f, w, indices);

            //get_weights(p, f, Eyflag, w, indices, plo, dx, dxi);
            p.rdata(FHD_realData::forcey) += p.rdata(FHD_realData::q)*inter_op(Eyarr, f, w, indices);

            //get_weights(p, f, Ezflag, w, indices, plo, dx, dxi);
            p.rdata(FHD_realData::forcez) += p.rdata(FHD_realData::q)*inter_op(Ezarr, f, w, indices);
        }
       
    });
    
}

void emf_gpu(FhdParticleContainer::AoS& aos,
                      FArrayBox& Ex, FArrayBox& Ey, FArrayBox& Ez,
                      const amrex::Real* plo_in, const amrex::Real* dx_in)
{
    if (common::pkernel_es[0] == 3)
    {
        emf_gpu(aos, Kernel4P(), Ex, Ey, Ez, plo_in, dx_in);
    }
    else if (common::pkernel_es[0] == 4)
    {
        emf_gpu(aos, Kernel4P(), Ex, Ey, Ez, plo_in, dx_in);
    }
    else if (common::pkernel_es[0] == 6)
    {
        emf_gpu(aos, Kernel6P(), Ex, Ey, Ez, plo_in, dx_in);
    }
}


template <typename F>
void spread_ions_fhd_gpu (FhdParticleContainer::AoS& aos, F f,
                      FArrayBox& sourcex, FArrayBox& sourcey, FArrayBox& sourcez,
                      const amrex::Real* plo_in,
                      const amrex::Real* dx_in)
{
    GpuArray<Real, 3> plo = {plo_in[0], plo_in[1], plo_in[2]};


    GpuArray<Real, 3> dx  = {dx_in[0],  dx_in[1],  dx_in[2] };    
    GpuArray<Real, 3> dxi = {1.0/dx_in[0],  1.0/dx_in[1],  1.0/dx_in[2] };

    Real volinv = 1.0/(dx_in[0]*dx_in[1]*dx_in[2]);


    const auto Np = aos.numParticles();
    const auto pstruct = aos().dataPtr();

    auto sourcexarr = sourcex.array();
    auto sourceyarr = sourcey.array();
    auto sourcezarr = sourcez.array();

    const IntVect sourcexflag = sourcex.box().type();
    const IntVect sourceyflag = sourcey.box().type();
    const IntVect sourcezflag = sourcez.box().type();


    constexpr int ks = F::ks;
       
    amrex::ParallelFor(Np, [=] AMREX_GPU_HOST_DEVICE (int ip) noexcept
    {
        FhdParticleContainer::ParticleType& p = pstruct[ip];

        Real forcex = p.rdata(FHD_realData::forcex)*volinv;
        Real forcey = p.rdata(FHD_realData::forcey)*volinv;
        Real forcez = p.rdata(FHD_realData::forcez)*volinv;
        
        Real w[3][2*ks];
        int indices[3][2*ks];
        
        get_weights(p, f, sourcexflag, w, indices, plo, dx, dxi);
        spread_op(sourcexarr, f, w, indices, forcex);

        get_weights(p, f, sourceyflag, w, indices, plo, dx, dxi);
        spread_op(sourceyarr, f, w, indices, forcey);

        get_weights(p, f, sourcezflag, w, indices, plo, dx, dxi);
        spread_op(sourcezarr, f, w, indices, forcez);
    });
    
}

void spread_ions_fhd_gpu(FhdParticleContainer::AoS& aos,
                      FArrayBox& sourcex, FArrayBox& sourcey, FArrayBox& sourcez,
                      const amrex::Real* plo_in,
                      const amrex::Real* dx_in)
{
    if (common::pkernel_es[0] == 3)
    {
        spread_ions_fhd_gpu(aos, Kernel4P(), sourcex, sourcey, sourcez, plo_in, dx_in);
    }
    else if (common::pkernel_es[0] == 4)
    {
        spread_ions_fhd_gpu(aos, Kernel4P(), sourcex, sourcey, sourcez, plo_in, dx_in);
    }
    else if (common::pkernel_es[0] == 6)
    {
        spread_ions_fhd_gpu(aos, Kernel6P(), sourcex, sourcey, sourcez, plo_in, dx_in);
    }
}

void mob_interp_gpu(Real z, Real a, Real* tmob, Real* nmob)
{
    if(wall_mob==1)
    {
        //Single plane approximation
        Real h = z/a;
        Real cn = 0.6452124576429801;
        Real ct = 0.635779246332329;
        *nmob = std::max(1.0 - 9.0/(8.0*(h+cn)) + 1.0/(2.0*pow(h+cn,3.0)) - 1.0/(8.0*pow(h+cn,5)) ,0.0);
        *tmob = std::max(1.0 - 9.0/(16.0*(h+ct)) + 2.0/(16.0*pow((h+ct),3.0)) -1.0/(16.0*pow((h+ct),5)),0.0);
    }else if(wall_mob==2)
    {
        Real h = z/a;
//    Measured using 256*256  !8.284925e-7, diff 1.326e-05
        *tmob = std::max(0.977222 - 10.3514/pow(2.34431 + h , 3) - 0.407416/pow(2.34431 + h , 1),0.0);
        *nmob = std::max(0.989214 + 152.697/pow(2.62177 + h,5) - 33.2256/pow(2.62177  + h,3) - 0.991641/(2.62177 + h),0.0);

    }
    else if(wall_mob==3)
    {
        Real h = (z-3.265E-7)/a;
//    OFFSET Measured using 256*256  !8.284925e-7, diff 1.326e-05
        *tmob = std::max(0.977222 - 10.3514/pow(2.34431 + h , 3) - 0.407416/pow(2.34431 + h , 1),0.0);
        *nmob = std::max(0.989214 + 152.697/pow(2.62177 + h,5) - 33.2256/pow(2.62177  + h,3) - 0.991641/(2.62177 + h),0.0);

    }
    else if(wall_mob==4)
    {
        //Single plane approximation shifted
        Real h = (z-3.265E-7)/a;
        Real cn = 0.6452124576429801;
        Real ct = 0.635779246332329;
        *nmob = std::max(1.0 - 9.0/(8.0*(h+cn)) + 1.0/(2.0*pow(h+cn,3.0)) - 1.0/(8.0*pow(h+cn,5)) ,0.0);
        *tmob = std::max(1.0 - 9.0/(16.0*(h+ct)) + 2.0/(16.0*pow((h+ct),3.0)) -1.0/(16.0*pow((h+ct),5)),0.0);

    }else if(wall_mob==5)
    {
        Real h = (z-3.265E-7)/a;
//  for freund problem
        *tmob = std::max(0.98557 + 13.9508/pow(1.82653  + h,5) - 6.83735/pow(1.82653   + h,3) - 1.00616/(1.82653  + h),0.0);
        *nmob = std::max(0.969255 + 163.299/pow(2.70459  + h,5) - 36.2912/pow(2.70459   + h,3) - 0.706417/(2.70459  + h),0.0);


    }else if(wall_mob==6)
    {
        Real h = (z-3.294E-7)/a;
//  for freund problem
        *tmob = std::max(0.98557 + 13.9508/pow(1.82653  + h,5) - 6.83735/pow(1.82653   + h,3) - 1.00616/(1.82653  + h),0.0);
        *nmob = std::max(0.969255 + 163.299/pow(2.70459  + h,5) - 36.2912/pow(2.70459   + h,3) - 0.706417/(2.70459  + h),0.0);

    }else if(wall_mob==7)
    {
        Real h = (z)/a;
//  for freund problem
        *tmob = std::max(0.98557 + 13.9508/pow(1.82653  + h,5) - 6.83735/pow(1.82653   + h,3) - 1.00616/(1.82653  + h),0.0);
        *nmob = std::max(0.969255 + 163.299/pow(2.70459  + h,5) - 36.2912/pow(2.70459   + h,3) - 0.706417/(2.70459  + h),0.0);

    }
}

void mob_interp_der_gpu(Real z, Real a, Real* tmobDer, Real* nmobDer)
{


    if(wall_mob==1)
    {
        //Single plane approximation
        Real cn = 0.6452124576429801;
        Real ct = 0.635779246332329;
        *nmobDer = 5.0/(8.*a*pow(cn + z/a,6)) - 3/(2.*a*pow(cn + z/a,4)) + 9/(8.*a*pow(cn + z/a,2));
        *tmobDer = 5/(16.*a*pow(ct + z/a,6)) - 3/(8.*a*pow(ct + z/a,4)) + 9/(16.*a*pow(ct + z/a,2));
    }else if(wall_mob==2)
    {
//  Measured using 256*256  !8.284925e-7, diff 1.326e-05
        Real h = z/a;

        *tmobDer = std::max(31.0541/(a*pow(2.34431 + h , 4)) + 0.407416/(a*pow(2.34431 + h , 2)),0.0);
        *nmobDer = std::max(-763.486/(a*pow(2.62177 + h,6)) + 99.6767/(a*pow(2.62177  + h,4)) + 0.991641/(a*pow(2.62177 + h,2)),0.0);
    }else if(wall_mob==3)
    {
//  Measured using 256*256  !8.284925e-7, diff 1.326e-05
        Real h = (z-3.265E-7)/a;

        *tmobDer = std::max(31.0541/(a*pow(2.34431 + h , 4)) + 0.407416/(a*pow(2.34431 + h , 2)),0.0);
        *nmobDer = std::max(-763.486/(a*pow(2.62177 + h,6)) + 99.6767/(a*pow(2.62177  + h,4)) + 0.991641/(a*pow(2.62177 + h,2)),0.0);
    }else if(wall_mob==4)
    {
        //Single plane approximation shifted
        Real cn = 0.6452124576429801;
        Real ct = 0.635779246332329;
        Real h = (z-3.265E-7)/a;

        *nmobDer = 5.0/(8.*a*pow(cn + h,6)) - 3/(2.*a*pow(cn + h,4)) + 9/(8.*a*pow(cn + h,2));
        *tmobDer = 5/(16.*a*pow(ct + h,6)) - 3/(8.*a*pow(ct + h,4)) + 9/(16.*a*pow(ct + h,2));
    }else if(wall_mob==5)
    {
//  for freund problem
        Real h = (z-3.265E-7)/a;

        *tmobDer = std::max(-69.754/(a*pow(1.82653 + h , 6)) + 20.5121/(a*pow(1.82653 + h , 4)) + 1.00616/(a*pow(1.82653 + h , 2)),0.0);
        *nmobDer = std::max(-816.495/(a*pow(2.70459 + h,6)) + 108.874/(a*pow(2.70459  + h,4)) + 0.706417/(a*pow(2.70459 + h,2)),0.0);
    }else if(wall_mob==6)
    {
//  for freund problem
        Real h = (z-3.294E-7)/a;

        *tmobDer = std::max(-69.754/(a*pow(1.82653 + h , 6)) + 20.5121/(a*pow(1.82653 + h , 4)) + 1.00616/(a*pow(1.82653 + h , 2)),0.0);
        *nmobDer = std::max(-816.495/(a*pow(2.70459 + h,6)) + 108.874/(a*pow(2.70459  + h,4)) + 0.706417/(a*pow(2.70459 + h,2)),0.0);
    }else if(wall_mob==7)
    {
//  for freund problem
        Real h = (z)/a;

        *tmobDer = std::max(-69.754/(a*pow(1.82653 + h , 6)) + 20.5121/(a*pow(1.82653 + h , 4)) + 1.00616/(a*pow(1.82653 + h , 2)),0.0);
        *nmobDer = std::max(-816.495/(a*pow(2.70459 + h,6)) + 108.874/(a*pow(2.70459  + h,4)) + 0.706417/(a*pow(2.70459 + h,2)),0.0);
    }


}

void get_mobility_diff_gpu(Real* nmob, Real* tmob, Real* nmobDer, Real* tmobDer, FhdParticleContainer::ParticleType& part, Real z)
{

    //using namespace common;

    Real awet = k_B*T_init[0]/(part.rdata(FHD_realData::wetDiff)*visc_coef*M_PI*6.0);
    Real atotal = k_B*T_init[0]/(part.rdata(FHD_realData::totalDiff)*visc_coef*M_PI*6.0);

    Real hwet = z/awet;
    Real htotal = z/atotal;

    Real tmobwet;
    Real nmobwet;
    Real tmobtotal;
    Real nmobtotal;

    Real tmobwetDer;
    Real nmobwetDer;
    Real tmobtotalDer;
    Real nmobtotalDer;

    //The mobility is dimensionless but the derivative is dimensional. Fix this at some point.

    mob_interp_gpu(z, awet, &tmobwet, &nmobwet);
    mob_interp_gpu(z, atotal, &tmobtotal, &nmobtotal);

    mob_interp_der_gpu(z, awet, &tmobwetDer, &nmobwetDer);
    mob_interp_der_gpu(z, atotal, &tmobtotalDer, &nmobtotalDer);

    *tmob = std::max((tmobtotal*part.rdata(FHD_realData::totalDiff) - tmobwet*part.rdata(FHD_realData::wetDiff))/part.rdata(FHD_realData::dryDiff),0.0);
    *nmob = std::max((nmobtotal*part.rdata(FHD_realData::totalDiff) - nmobwet*part.rdata(FHD_realData::wetDiff))/part.rdata(FHD_realData::dryDiff),0.0);

    *nmobDer = (nmobtotalDer*part.rdata(FHD_realData::totalDiff) - nmobwetDer*part.rdata(FHD_realData::wetDiff))/part.rdata(FHD_realData::dryDiff);
    *tmobDer = (tmobtotalDer*part.rdata(FHD_realData::totalDiff) - tmobwetDer*part.rdata(FHD_realData::wetDiff))/part.rdata(FHD_realData::dryDiff);
}


void get_explicit_mobility_gpu(Real* mob, Real* mobDer, FhdParticleContainer::ParticleType& part, const Real* plo, const Real* phi)
{                           

    Real nmob;
    Real tmob;
    Real nmobDer;
    Real tmobDer;
    Real z;

    mob[0] = 1;
    mob[1] = 1;
    mob[2] = 1;

    mobDer[0] = 0;
    mobDer[1] = 0;
    mobDer[2] = 0;


    if((bc_vel_lo[0] == 2) && (bc_vel_hi[0] == 2))
    {

       z = part.pos(0);

       if(z > (phi[0]-plo[0])/2.0)
       {
          z = phi[0] - z;
       }

       get_mobility_diff_gpu(&nmob, &tmob, &nmobDer, &tmobDer, part, z);

       mob[0] = nmob;
       mob[1] = tmob;               
       mob[2] = tmob;

       mobDer[0] = nmobDer;
       mobDer[1] = tmobDer;               
       mobDer[2] = tmobDer;

       z = part.pos(0);
       if(z > (phi[0]-plo[0])/2.0)
       {
           mobDer[0] = -mobDer[0];
           mobDer[1] = -mobDer[1];
           mobDer[2] = -mobDer[2];
       }
    }

    if((bc_vel_lo[1] == 2) && (bc_vel_hi[1] == 2))
    {
       z = part.pos(1);

       if(z > (phi[1]-plo[1])/2.0)
       {
          z = phi[1] - z;
       }

       get_mobility_diff_gpu(&nmob, &tmob, &nmobDer, &tmobDer, part, z);

       mob[0] = tmob;
       mob[1] = nmob;               
       mob[2] = tmob;

       mobDer[0] = tmobDer;
       mobDer[1] = nmobDer;               
       mobDer[2] = tmobDer;

       z = part.pos(1);
       if(z > (phi[1]-plo[1])/2.0)
       {
           mobDer[0] = -mobDer[0];
           mobDer[1] = -mobDer[1];
           mobDer[2] = -mobDer[2];
       }

        

    }

    if((bc_vel_lo[2] == 2) && (bc_vel_hi[2] == 2))
    {

       z = z = part.pos(2);;

       if(z > (phi[2]-plo[2])/2.0)
       {
          z = phi[2] - z;
       }

       get_mobility_diff_gpu(&nmob, &tmob, &nmobDer, &tmobDer, part, z);

       mob[0] = tmob;
       mob[1] = tmob;               
       mob[2] = nmob;

       mobDer[0] = tmobDer;
       mobDer[1] = tmobDer;               
       mobDer[2] = nmobDer;

       z = part.pos(2);
       if(z > (phi[2]-plo[2])/2.0)
       {
           mobDer[0] = -mobDer[0];
           mobDer[1] = -mobDer[1];
           mobDer[2] = -mobDer[2];
       }
    }

    //Print() << "z, nmob, pos: " << z << ", " << nmob << ", " << part.pos(2) << "\n";
    //cin.get();
}

void dry_gpu(Real dt, FhdParticleContainer::ParticleType& part, Real* dry_terms, Real* mb, Real* mobDir)
{
    Real normalrand[3];
    Real std[3];
    Real bfac[3];

    normalrand[0] = amrex::RandomNormal(0.,1.);
    normalrand[1] = amrex::RandomNormal(0.,1.);
    normalrand[2] = amrex::RandomNormal(0.,1.);

//    !std = sqrt(part%dry_diff*k_B*2d0*t_init(1))
    std[0] = sqrt(2.0*mb[0]*part.rdata(FHD_realData::dryDiff));
    std[1] = sqrt(2.0*mb[1]*part.rdata(FHD_realData::dryDiff));
    std[2] = sqrt(2.0*mb[2]*part.rdata(FHD_realData::dryDiff));

    //DRL: dry diffusion coef: part%dry_diff, temperature: t_init(1)

    bfac[0] = variance_coef_mom*std[0]*normalrand[0]/sqrt(dt);
    bfac[1] = variance_coef_mom*std[1]*normalrand[1]/sqrt(dt);
    bfac[2] = variance_coef_mom*std[2]*normalrand[2]/sqrt(dt);

    dry_terms[0] = mb[0]*part.rdata(FHD_realData::dryDiff)*part.rdata(FHD_realData::forcex)/(k_B*T_init[0])+bfac[0];
    dry_terms[1] = mb[1]*part.rdata(FHD_realData::dryDiff)*part.rdata(FHD_realData::forcey)/(k_B*T_init[0])+bfac[1];
    dry_terms[2] = mb[2]*part.rdata(FHD_realData::dryDiff)*part.rdata(FHD_realData::forcez)/(k_B*T_init[0])+bfac[2];

    //std::cout << "terms: " << dry_terms[1] << ", " << mobDir[1]*part.rdata(FHD_realData::dryDiff) << std::endl;

    if(dry_move_tog == 1)
    {
        dry_terms[0] = dry_terms[0] + variance_coef_mom*mobDir[0]*part.rdata(FHD_realData::dryDiff);
        dry_terms[1] = dry_terms[1] + variance_coef_mom*mobDir[1]*part.rdata(FHD_realData::dryDiff);
        dry_terms[2] = dry_terms[2] + variance_coef_mom*mobDir[2]*part.rdata(FHD_realData::dryDiff);
    }

}

#endif
