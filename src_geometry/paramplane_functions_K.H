#ifndef _paramplane_functions_K_H_
#define _paramplane_functions_K_H_

#include <AMReX.H>
#include <common_namespace.H>
#include <FhdParticleContainer.H>

void pre_check(FhdParticleContainer::ParticleType& part, const amrex_real delt, const paramPlane* paramplanes, 
                                                const int ns, int* flag, const amrex_real* phi, const amrex_real* plo)
{

    Real proj[3];
    Real box1lo[3];
    Real box2lo[3];
    Real box1hi[3];
    Real box2hi[3];

    for (int d=0; d<AMREX_SPACEDIM; ++d)
    {
        proj[d] = part.pos(d) + delt*part.rdata(FHD_realData::velx + 0);
    }

    for (int d=0; d<AMREX_SPACEDIM; ++d)
    {
        box1lo[d] = plo[d];
        box1hi[d] = phi[d];

        box2lo[d] = plo[d];
        box2hi[d] = phi[d];
    }

#if (BL_SPACEDIM == 3)
    if(  (part.pos(0) < box1hi[0]) && (part.pos(1) < box1hi[1]) && (part.pos(2) < box1hi[2]) && (part.pos(0) > box1lo[0]) && (part.pos(1) > box1lo[1]) && (part.pos(2) > box1lo[2])  )
    {
        if(  (proj[0] < box1hi[0]) && (proj[1] < box1hi[1]) && (proj[2] < box1hi[2]) && (proj[0] > box1lo[0]) && (proj[1] > box1lo[1]) && (proj[2] > box1lo[2]) )
        {
            *flag = 1;
        }
    }else  //started in box 2
    {
        if(  (proj[0] < box2hi[0]) && (proj[1] < box2hi[1]) && (proj[2] < box2hi[2]) && (proj[0] > box2lo[0]) && (proj[1] > box2lo[1]) && (proj[2] > box2lo[2]) )
        {
            *flag = 1;
        }
    }
#endif
}

void find_inter_gpu(FhdParticleContainer::ParticleType& part, const amrex_real delt, const paramPlane* paramplanes, const int ns, int* intsurf,
                        amrex_real* inttime, int* intside, const amrex_real* phi, const amrex_real* plo)
{

    int flag = 0;
    *inttime = delt;
    *intsurf = -1;
    Real uval, vval, tval;
    
    //pre_check(part, delt, paramplanes, ns, &flag, phi, plo);

    if(flag == 0)
    {
        for(int s=0;s<ns;s++)
        {
            const paramPlane* surf = &paramplanes[s];

            Real denominv = 1.0/(part.rdata(FHD_realData::velz)*surf->uy*surf->vx - part.rdata(FHD_realData::vely)*surf->uz*surf->vx - part.rdata(FHD_realData::velz)*surf->ux*surf->vy + part.rdata(FHD_realData::velx)*surf->uz*surf->vy + part.rdata(FHD_realData::vely)*surf->ux*surf->vz - part.rdata(FHD_realData::velx)*surf->uy*surf->vz);

            uval = (part.rdata(FHD_realData::velz)*part.pos(1)*surf->vx - part.rdata(FHD_realData::vely)*part.pos(2)*surf->vx - part.rdata(FHD_realData::velz)*surf->y0*surf->vx + part.rdata(FHD_realData::vely)*surf->z0*surf->vx - part.rdata(FHD_realData::velz)*part.pos(0)*surf->vy + part.rdata(FHD_realData::velx)*part.pos(2)*surf->vy + part.rdata(FHD_realData::velz)*surf->x0*surf->vy - part.rdata(FHD_realData::velx)*surf->z0*surf->vy + part.rdata(FHD_realData::vely)*part.pos(0)*surf->vz - part.rdata(FHD_realData::velx)*part.pos(1)*surf->vz -  part.rdata(FHD_realData::vely)*surf->x0*surf->vz + part.rdata(FHD_realData::velx)*surf->y0*surf->vz)*denominv;

            vval = (-part.rdata(FHD_realData::velz)*part.pos(1)*surf->ux + part.rdata(FHD_realData::vely)*part.pos(2)*surf->ux + part.rdata(FHD_realData::velz)*surf->y0*surf->ux - part.rdata(FHD_realData::vely)*surf->z0*surf->ux + part.rdata(FHD_realData::velz)*part.pos(0)*surf->uy - part.rdata(FHD_realData::velx)*part.pos(2)*surf->uy - part.rdata(FHD_realData::velz)*surf->x0*surf->uy + part.rdata(FHD_realData::velx)*surf->z0*surf->uy - part.rdata(FHD_realData::vely)*part.pos(0)*surf->uz + part.rdata(FHD_realData::velx)*part.pos(1)*surf->uz + part.rdata(FHD_realData::vely)*surf->x0*surf->uz - part.rdata(FHD_realData::velx)*surf->y0*surf->uz)*denominv;

            tval = (-part.pos(2)*surf->uy*surf->vx + surf->z0*surf->uy*surf->vx + part.pos(1)*surf->uz*surf->vx - surf->y0*surf->uz*surf->vx + part.pos(2)*surf->ux*surf->vy - surf->z0*surf->ux*surf->vy - part.pos(0)*surf->uz*surf->vy + surf->x0*surf->uz*surf->vy - part.pos(1)*surf->ux*surf->vz + surf->y0*surf->ux*surf->vz + part.pos(0)*surf->uy*surf->vz - surf->x0*surf->uy*surf->vz)*denominv;

            //Print() << "Testing "<< s+1 << " particle " << part.id() << ", " << uval << ", " << vval << ", " << tval << ", " << denominv << "\n";
            //Print() << "Vel "<< part.rdata(FHD_realData::velx) << ", " << part.rdata(FHD_realData::vely) << ", " << part.rdata(FHD_realData::velz) << "\n";
            //Print() << "surf1 "<< surf->ux << ", " << surf->uy << ", " << surf->uz << "\n";
            //Print() << "surf2 "<< surf->vx << ", " << surf->vy << ", " << surf->vz << "\n";
       
            if(  ((uval > 0) && (uval < surf->uTop)) && ((vval > 0) && (vval < surf->vTop))  &&  ((tval > 0) && (tval < *inttime))   )
            {
              *inttime = tval;
              *intsurf = s;

              Real dotprod = part.rdata(FHD_realData::velx)*surf->lnx + part.rdata(FHD_realData::vely)*surf->lny + part.rdata(FHD_realData::velz)*surf->lnz;

              if (dotprod > 0)
              {
                *intside = 1; //1 for rhs
              }else
              {
                *intside = 0; //0 for lhs
              }

              //std::cout << "Intersection! " << s << " part " << part.id() << ", " << part.rdata(FHD_realData::velx) << ", " << part.rdata(FHD_realData::vely) << ", " << part.rdata(FHD_realData::velz) << ", " << *intside << "\n";

            }
        }

    }

}


#endif
