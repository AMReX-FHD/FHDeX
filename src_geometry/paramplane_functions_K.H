#ifndef _paramplane_functions_K_H_
#define _paramplane_functions_K_H_

#include <AMReX.H>
#include <common_namespace.H>
#include <FhdParticleContainer.H>
#include <math.h>

void pre_check_gpu(FhdParticleContainer::ParticleType& part, const amrex_real delt, const paramPlane* paramplanes, 
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
    
    pre_check_gpu(part, delt, paramplanes, ns, &flag, phi, plo);

    if(flag == 0)
    {
        for(int s=1;s<=ns;s++)
        {
            const paramPlane* surf = &paramplanes[s-1];

            Real denominv = 1.0/(part.rdata(FHD_realData::velz)*surf->uy*surf->vx - part.rdata(FHD_realData::vely)*surf->uz*surf->vx - part.rdata(FHD_realData::velz)*surf->ux*surf->vy + part.rdata(FHD_realData::velx)*surf->uz*surf->vy + part.rdata(FHD_realData::vely)*surf->ux*surf->vz - part.rdata(FHD_realData::velx)*surf->uy*surf->vz);

            uval = (part.rdata(FHD_realData::velz)*part.pos(1)*surf->vx - part.rdata(FHD_realData::vely)*part.pos(2)*surf->vx - part.rdata(FHD_realData::velz)*surf->y0*surf->vx + part.rdata(FHD_realData::vely)*surf->z0*surf->vx - part.rdata(FHD_realData::velz)*part.pos(0)*surf->vy + part.rdata(FHD_realData::velx)*part.pos(2)*surf->vy + part.rdata(FHD_realData::velz)*surf->x0*surf->vy - part.rdata(FHD_realData::velx)*surf->z0*surf->vy + part.rdata(FHD_realData::vely)*part.pos(0)*surf->vz - part.rdata(FHD_realData::velx)*part.pos(1)*surf->vz -  part.rdata(FHD_realData::vely)*surf->x0*surf->vz + part.rdata(FHD_realData::velx)*surf->y0*surf->vz)*denominv;

            vval = (-part.rdata(FHD_realData::velz)*part.pos(1)*surf->ux + part.rdata(FHD_realData::vely)*part.pos(2)*surf->ux + part.rdata(FHD_realData::velz)*surf->y0*surf->ux - part.rdata(FHD_realData::vely)*surf->z0*surf->ux + part.rdata(FHD_realData::velz)*part.pos(0)*surf->uy - part.rdata(FHD_realData::velx)*part.pos(2)*surf->uy - part.rdata(FHD_realData::velz)*surf->x0*surf->uy + part.rdata(FHD_realData::velx)*surf->z0*surf->uy - part.rdata(FHD_realData::vely)*part.pos(0)*surf->uz + part.rdata(FHD_realData::velx)*part.pos(1)*surf->uz + part.rdata(FHD_realData::vely)*surf->x0*surf->uz - part.rdata(FHD_realData::velx)*surf->y0*surf->uz)*denominv;

            tval = (-part.pos(2)*surf->uy*surf->vx + surf->z0*surf->uy*surf->vx + part.pos(1)*surf->uz*surf->vx - surf->y0*surf->uz*surf->vx + part.pos(2)*surf->ux*surf->vy - surf->z0*surf->ux*surf->vy - part.pos(0)*surf->uz*surf->vy + surf->x0*surf->uz*surf->vy - part.pos(1)*surf->ux*surf->vz + surf->y0*surf->ux*surf->vz + part.pos(0)*surf->uy*surf->vz - surf->x0*surf->uy*surf->vz)*denominv;

//            Print() << "Testing "<< s+1 << " particle " << part.id() << ", " << uval << ", " << vval << ", " << tval << ", " << denominv << "\n";
//            Print() << "Vel "<< part.rdata(FHD_realData::velx) << ", " << part.rdata(FHD_realData::vely) << ", " << part.rdata(FHD_realData::velz) << "\n";
//            Print() << "surf1 "<< surf->ux << ", " << surf->uy << ", " << surf->uz << "\n";
//            Print() << "surf2 "<< surf->vx << ", " << surf->vy << ", " << surf->vz << "\n";
       
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

void rotation(Real costheta, Real sintheta, Real cosphi, Real sinphi, Real *cx, Real *cy, Real *cz)
{

    Real cxnew = (*cx)*costheta*cosphi + *cz*cosphi*sintheta - *cy*sinphi;
    *cy = *cy*cosphi + *cx*costheta*sinphi + *cz*sintheta*sinphi;
    *cz = *cz*costheta - *cx*sintheta;

    *cx = cxnew;

 }

void randomhemisphere(Real costheta, Real sintheta, Real cosphi, Real sinphi, Real *cx, Real *cy, Real *cz)
{
    Real mag, costhetanew, sinthetanew, cosphinew, sinphinew, phi;

    mag = sqrt((*cx)*(*cx) + (*cy)*(*cy) + (*cz)*(*cz));

//    call get_half_angles(costhetanew, sinthetanew, cosphinew, sinphinew)

    costhetanew = amrex::Random();
    sinthetanew = sqrt(1.0 - costhetanew*costhetanew);

    phi = amrex::Random()*2.0*M_PI;

    cosphinew = cos(phi);
    sinphinew = sin(phi);

    *cx = mag*sinthetanew*cosphinew;
    *cy = mag*sinthetanew*sinphinew;
    *cz = mag*costhetanew;

    rotation(costheta, sintheta, cosphi, sinphi, cx, cy, cz);

}

void app_bc_gpu(const paramPlane* surf, FhdParticleContainer::ParticleType& part, int intside, Real* domsize, int *push, Real time, Real inttime)
{
    
    
    if(intside == 1)
    {
      if(amrex::Random() > surf->porosityRight)
      {
        *push = 0;

//        surf%fxright = surf%fxright + part%mass*part%vel(1)
//        surf%fyright = surf%fyright + part%mass*part%vel(2)
//        surf%fzright = surf%fzright + part%mass*part%vel(3)

        if(get_uniform_func() < surf->specularityRight)
        {

          Real dotprod = part.rdata(FHD_realData::velx)*surf->rnx + part.rdata(FHD_realData::vely)*surf->rny + part.rdata(FHD_realData::velz)*surf->rnz;

          part.rdata(FHD_realData::velx) = -2.0*dotprod*surf->rnx + part.rdata(FHD_realData::velx);
          part.rdata(FHD_realData::vely) = -2.0*dotprod*surf->rny + part.rdata(FHD_realData::vely);
          part.rdata(FHD_realData::velz) = -2.0*dotprod*surf->rnz + part.rdata(FHD_realData::velz);
        }else
        {
          
          Real srt = sqrt(part.rdata(FHD_realData::R)*surf->temperatureRight);

          part.rdata(FHD_realData::velx) = srt*amrex::RandomNormal(0.,1.);
          part.rdata(FHD_realData::vely) = srt*amrex::RandomNormal(0.,1.);
          part.rdata(FHD_realData::velz) = 1.414213562*srt*sqrt(-log(amrex::Random()));
        
          rotation(surf->cosThetaRight, surf->sinThetaRight, surf->cosPhiRight, surf->sinPhiRight, &part.rdata(FHD_realData::velx), &part.rdata(FHD_realData::vely), &part.rdata(FHD_realData::velz));


        }
      }else if(amrex::Random() < surf->periodicity)
      {
        *push = 0;

        if(surf->boundary == 1)
        {
          part.pos(0) = part.pos(0) + 0.9999999*domsize[0];
        }
        else if(surf->boundary == 2)
        {
          part.pos(0) = part.pos(0) - 0.9999999*domsize[0];
        }
        else if(surf->boundary == 3)
        {
          part.pos(1) = part.pos(1) + 0.9999999*domsize[1];
        }
        else if(surf->boundary == 4)
        {
          part.pos(1) = part.pos(1) - 0.9999999*domsize[1];
        }
        else if(surf->boundary == 5)
        {
          part.pos(2) = part.pos(2) + 0.9999999*domsize[2];
        }
        else if(surf->boundary == 6)
        {
          part.pos(2) = part.pos(2) - 0.9999999*domsize[2];
        }

      }else
      {
        *push = 1;

        if(amrex::Random() > surf->momentumConsRight)
        {
          randomhemisphere(surf->cosThetaRight, surf->sinThetaRight, surf->cosPhiRight, surf->sinPhiRight, &part.rdata(FHD_realData::velx), &part.rdata(FHD_realData::vely), &part.rdata(FHD_realData::velz));
        }
      }

    }else
    {

//      if(get_uniform_func() > surf%porosityleft)
//      {
//        *push = 0;

//        surf%fxleft = surf%fxleft + part%mass*part%vel(1)
//        surf%fyleft = surf%fyleft + part%mass*part%vel(1)
//        surf%fzleft = surf%fzleft + part%mass*part%vel(1)

//        if(get_uniform_func() < surf%specularityleft)
//        {
//          dotprod = part%vel(1)*surf%lnx + part%vel(2)*surf%lny + part%vel(3)*surf%lnz

//          part%vel(1) = -2d0*dotprod*surf%lnx + part%vel(1)
//          part%vel(2) = -2d0*dotprod*surf%lny + part%vel(2)
//          part%vel(3) = -2d0*dotprod*surf%lnz + part%vel(3)
//        }
//        else
//        { 
//          srt = sqrt(part%r*surf%temperatureleft)

//          part%vel(1) = srt*get_particle_normal_func()
//          part%vel(2) = srt*get_particle_normal_func()
//          part%vel(3) = 1.414213562*srt*sqrt(-log(get_uniform_func()))

//          call rotation(surf%costhetaleft, surf%sinthetaleft, surf%cosphileft, surf%sinphileft, part%vel(1), part%vel(2), part%vel(3))
//        }
       
//      }
//      else if(amrex::Random() < surf->periodicity)
//      {
//        push = 0

//        if(surf%boundary .eq. 1) then

//          part%pos(1) = part%pos(1) + 0.9999999*domsize(1)

//        elseif(surf%boundary .eq. 2) then

//          part%pos(1) = part%pos(1) - 0.9999999*domsize(1)

//        elseif(surf%boundary .eq. 3) then

//          part%pos(2) = part%pos(2) + 0.9999999*domsize(2)

//        elseif(surf%boundary .eq. 4) then

//          part%pos(2) = part%pos(2) - 0.9999999*domsize(2)

//#if (BL_SPACEDIM == 3)
//        elseif(surf%boundary .eq. 5) then

//          part%pos(3) = part%pos(3) + 0.9999999*domsize(3)

//        elseif(surf%boundary .eq. 6) then

//          part%pos(3) = part%pos(3) - 0.9999999*domsize(3)
//#endif
//        endif

//      else
//        push = 1      

//        if(get_uniform_func() > surf%momentumleft) then
//#if (BL_SPACEDIM == 3)
//          call randomhemisphere(surf%costhetaright, surf%sinthetaright, surf%cosphiright, surf%sinphiright, part%vel(1), part%vel(2), part%vel(3))
//#endif
//#if (BL_SPACEDIM == 2)
//          call randomhemisphere(surf%costhetaright, surf%sinthetaright, part%vel(1), part%vel(2), part%vel(3))
//#endif

//        endif

//      endif
  
   }
//   if(graphene_tog .eq. 1) then
//      if(surf%boundary .eq. 6) then

//      !print *, "new: ", part%vel(3)
//      call surf_velocity(surf, part, time, oldvel, inttime)
//   endif
//   endif
        
}

#endif