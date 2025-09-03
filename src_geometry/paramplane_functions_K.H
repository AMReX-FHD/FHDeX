#ifndef _paramplane_functions_K_H_
#define _paramplane_functions_K_H_

#include <AMReX.H>
#include <common_namespace.H>
#include <math.h>

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void pre_check_gpu(FhdParticleContainer::ParticleType& part, const Real delt, const paramPlane* paramplanes,
                                                const int ns, int* flag, const GpuArray<Real, 3>& phi, const GpuArray<Real, 3>& plo, Real* inttime)
{
    Real proj[3];
    Real box1lo[3];
    Real box2lo[3];
    Real box1hi[3];
    Real box2hi[3];

    for (int d=0; d<(AMREX_SPACEDIM); ++d)
    {
        proj[d] = part.pos(d) + delt*part.rdata(FHD_realData::velx + d);
    }

    for (int d=0; d<AMREX_SPACEDIM; ++d)
    {
        box1lo[d] = plo[d];
        box1hi[d] = phi[d];

    }

    if(  (proj[0] < box1hi[0]) && (proj[1] < box1hi[1]) &&
        (proj[2] < box1hi[2]) && (proj[0] > box1lo[0]) &&
        (proj[1] > box1lo[1]) && (proj[2] > box1lo[2]) )
    {
        *flag = 1;
    }
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void dummy(Real* inttime)
{
    printf("delt dummy %e\n",*inttime);
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void find_inter_gpu(FhdParticleContainer::ParticleType& part, const Real delt, const paramPlane* paramplanes, const int ns, int* intsurf,
                        Real* inttime, int* intside, const GpuArray<Real, 3>& phi, const GpuArray<Real, 3>& plo)
{

    //BL_PROFILE_VAR("find_inter_gpu()", find_inter_gpu);
    int flag = 0;
    *inttime = delt;
    *intsurf = -1;
    Real uval, vval, tval;

    //pre_check_gpu(part, delt, paramplanes, ns, &flag, phi, plo, inttime);

    //Complete
    if(flag == 0)
    {
        for(int s=1;s<=ns;s++)
        {
            //if((s != 2) && (s != 3))
            {
                const paramPlane* surf = &paramplanes[s-1];

                Real denominv = 1.0/(part.rdata(FHD_realData::velz)*surf->uy*surf->vx - part.rdata(FHD_realData::vely)*surf->uz*surf->vx - part.rdata(FHD_realData::velz)*surf->ux*surf->vy + part.rdata(FHD_realData::velx)*surf->uz*surf->vy + part.rdata(FHD_realData::vely)*surf->ux*surf->vz - part.rdata(FHD_realData::velx)*surf->uy*surf->vz);

                uval = (part.rdata(FHD_realData::velz)*part.pos(1)*surf->vx - part.rdata(FHD_realData::vely)*part.pos(2)*surf->vx - part.rdata(FHD_realData::velz)*surf->y0*surf->vx + part.rdata(FHD_realData::vely)*surf->z0*surf->vx - part.rdata(FHD_realData::velz)*part.pos(0)*surf->vy + part.rdata(FHD_realData::velx)*part.pos(2)*surf->vy + part.rdata(FHD_realData::velz)*surf->x0*surf->vy - part.rdata(FHD_realData::velx)*surf->z0*surf->vy + part.rdata(FHD_realData::vely)*part.pos(0)*surf->vz - part.rdata(FHD_realData::velx)*part.pos(1)*surf->vz -  part.rdata(FHD_realData::vely)*surf->x0*surf->vz + part.rdata(FHD_realData::velx)*surf->y0*surf->vz)*denominv;

                vval = (-part.rdata(FHD_realData::velz)*part.pos(1)*surf->ux + part.rdata(FHD_realData::vely)*part.pos(2)*surf->ux + part.rdata(FHD_realData::velz)*surf->y0*surf->ux - part.rdata(FHD_realData::vely)*surf->z0*surf->ux + part.rdata(FHD_realData::velz)*part.pos(0)*surf->uy - part.rdata(FHD_realData::velx)*part.pos(2)*surf->uy - part.rdata(FHD_realData::velz)*surf->x0*surf->uy + part.rdata(FHD_realData::velx)*surf->z0*surf->uy - part.rdata(FHD_realData::vely)*part.pos(0)*surf->uz + part.rdata(FHD_realData::velx)*part.pos(1)*surf->uz + part.rdata(FHD_realData::vely)*surf->x0*surf->uz - part.rdata(FHD_realData::velx)*surf->y0*surf->uz)*denominv;

                tval = (-part.pos(2)*surf->uy*surf->vx + surf->z0*surf->uy*surf->vx + part.pos(1)*surf->uz*surf->vx - surf->y0*surf->uz*surf->vx + part.pos(2)*surf->ux*surf->vy - surf->z0*surf->ux*surf->vy - part.pos(0)*surf->uz*surf->vy + surf->x0*surf->uz*surf->vy - part.pos(1)*surf->ux*surf->vz + surf->y0*surf->ux*surf->vz + part.pos(0)*surf->uy*surf->vz - surf->x0*surf->uy*surf->vz)*denominv;

                //Print() << "Checking particle " << part.id() << ", wall " << s << endl;
                //Print() << uval/surf->uTop << ", " << vval/surf->vTop << ", " << tval/delt << endl;

                if(  ((uval > 0) && (uval < surf->uTop)) && ((vval > 0) && (vval < surf->vTop))  &&  ((tval > 0) && (tval < *inttime))   )
                {
                    *inttime = tval;
                    *intsurf = s;

                    Real dotprod = part.rdata(FHD_realData::velx)*surf->lnx +
                        part.rdata(FHD_realData::vely)*surf->lny + part.rdata(FHD_realData::velz)*surf->lnz;

                    if (dotprod > 0)
                    {
                        *intside = 1; //1 for rhs
                    }
                    else
                    {
                        *intside = 0; //0 for lhs
                    }
                }
            }
        }
    }

    //No z velocity
/*    if(flag == 0)*/
/*    {*/
/*        for(int s=1;s<=ns;s++)*/
/*        {*/
/*            const paramPlane* surf = &paramplanes[s-1];*/

/*            Real denominv = 1.0/( - part.rdata(FHD_realData::vely)*surf->uz*surf->vx + part.rdata(FHD_realData::velx)*surf->uz*surf->vy + part.rdata(FHD_realData::vely)*surf->ux*surf->vz - part.rdata(FHD_realData::velx)*surf->uy*surf->vz);*/

/*            uval = (- part.rdata(FHD_realData::vely)*part.pos(2)*surf->vx + part.rdata(FHD_realData::vely)*surf->z0*surf->vx + part.rdata(FHD_realData::velx)*part.pos(2)*surf->vy - part.rdata(FHD_realData::velx)*surf->z0*surf->vy + part.rdata(FHD_realData::vely)*part.pos(0)*surf->vz - part.rdata(FHD_realData::velx)*part.pos(1)*surf->vz -  part.rdata(FHD_realData::vely)*surf->x0*surf->vz + part.rdata(FHD_realData::velx)*surf->y0*surf->vz)*denominv;*/

/*            vval = (part.rdata(FHD_realData::vely)*part.pos(2)*surf->ux - part.rdata(FHD_realData::vely)*surf->z0*surf->ux - part.rdata(FHD_realData::velx)*part.pos(2)*surf->uy + part.rdata(FHD_realData::velx)*surf->z0*surf->uy - part.rdata(FHD_realData::vely)*part.pos(0)*surf->uz + part.rdata(FHD_realData::velx)*part.pos(1)*surf->uz + part.rdata(FHD_realData::vely)*surf->x0*surf->uz - part.rdata(FHD_realData::velx)*surf->y0*surf->uz)*denominv;*/

/*            tval = (-part.pos(2)*surf->uy*surf->vx + surf->z0*surf->uy*surf->vx + part.pos(1)*surf->uz*surf->vx - surf->y0*surf->uz*surf->vx + part.pos(2)*surf->ux*surf->vy - surf->z0*surf->ux*surf->vy - part.pos(0)*surf->uz*surf->vy + surf->x0*surf->uz*surf->vy - part.pos(1)*surf->ux*surf->vz + surf->y0*surf->ux*surf->vz + part.pos(0)*surf->uy*surf->vz - surf->x0*surf->uy*surf->vz)*denominv;*/

/*                       //Print() << "Checking particle " << part.id() << ", wall " << s << endl;*/
/*                       //Print() << uval/surf->uTop << ", " << vval/surf->vTop << ", " << tval/delt << endl;*/
/*       */
/*            if(  ((uval > 0) && (uval < surf->uTop)) && ((vval > 0) && (vval < surf->vTop))  &&  ((tval > 0) && (tval < *inttime))   )*/
/*            {*/
/*                *inttime = tval;*/
/*                *intsurf = s;*/

/*                Real dotprod = part.rdata(FHD_realData::velx)*surf->lnx + */
/*                    part.rdata(FHD_realData::vely)*surf->lny;*/

/*                if (dotprod > 0)*/
/*                {*/
/*                    *intside = 1; //1 for rhs*/
/*                }*/
/*                else*/
/*                {*/
/*                    *intside = 0; //0 for lhs*/
/*                }*/
/*            }*/
/*        }*/
/*    }*/

    //No y or z velocity
/*    if(flag == 0)*/
/*    {*/
/*        for(int s=1;s<=ns;s++)*/
/*        {*/
/*            const paramPlane* surf = &paramplanes[s-1];*/

/*            Real denominv = 1.0/( part.rdata(FHD_realData::velx)*surf->uz*surf->vy - part.rdata(FHD_realData::velx)*surf->uy*surf->vz);*/

/*            uval = (part.rdata(FHD_realData::velx)*part.pos(2)*surf->vy - part.rdata(FHD_realData::velx)*surf->z0*surf->vy - part.rdata(FHD_realData::velx)*part.pos(1)*surf->vz  + part.rdata(FHD_realData::velx)*surf->y0*surf->vz)*denominv;*/

/*            vval = (- part.rdata(FHD_realData::velx)*part.pos(2)*surf->uy + part.rdata(FHD_realData::velx)*surf->z0*surf->uy + part.rdata(FHD_realData::velx)*part.pos(1)*surf->uz - part.rdata(FHD_realData::velx)*surf->y0*surf->uz)*denominv;*/

/*            tval = (-part.pos(2)*surf->uy*surf->vx + surf->z0*surf->uy*surf->vx + part.pos(1)*surf->uz*surf->vx - surf->y0*surf->uz*surf->vx + part.pos(2)*surf->ux*surf->vy - surf->z0*surf->ux*surf->vy - part.pos(0)*surf->uz*surf->vy + surf->x0*surf->uz*surf->vy - part.pos(1)*surf->ux*surf->vz + surf->y0*surf->ux*surf->vz + part.pos(0)*surf->uy*surf->vz - surf->x0*surf->uy*surf->vz)*denominv;*/

/*                       //Print() << "Checking particle " << part.id() << ", wall " << s << endl;*/
/*                       //Print() << uval/surf->uTop << ", " << vval/surf->vTop << ", " << tval/delt << endl;*/
/*       */
/*            if(  ((uval > 0) && (uval < surf->uTop)) && ((vval > 0) && (vval < surf->vTop))  &&  ((tval > 0) && (tval < *inttime))   )*/
/*            {*/
/*                *inttime = tval;*/
/*                *intsurf = s;*/

/*                Real dotprod = part.rdata(FHD_realData::velx)*surf->lnx;*/

/*                if (dotprod > 0)*/
/*                {*/
/*                    *intside = 1; //1 for rhs*/
/*                }*/
/*                else*/
/*                {*/
/*                    *intside = 0; //0 for lhs*/
/*                }*/
/*            }*/
/*        }*/
/*    }*/
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rotation(Real costheta, Real sintheta, Real cosphi, Real sinphi, Real *cx, Real *cy, Real *cz)
{
    Real cxnew = (*cx)*costheta*cosphi + *cz*cosphi*sintheta - *cy*sinphi;
    *cy = *cy*cosphi + *cx*costheta*sinphi + *cz*sintheta*sinphi;
    *cz = *cz*costheta - *cx*sintheta;
    *cx = cxnew;
}

// uncomment this once this routine is on a GPU device
// then the amrex::Random() calls will need to be fixed
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void randomhemisphere(Real costheta, Real sintheta, Real cosphi, Real sinphi, Real *cx, Real *cy, Real *cz, amrex::RandomEngine const& engine)
{
    Real mag, costhetanew, sinthetanew, cosphinew, sinphinew, phi;

    mag = sqrt((*cx)*(*cx) + (*cy)*(*cy) + (*cz)*(*cz));

    //    call get_half_angles(costhetanew, sinthetanew, cosphinew, sinphinew)

    costhetanew = amrex::Random(engine);
    sinthetanew = sqrt(1.0 - costhetanew*costhetanew);

    phi = amrex::Random(engine)*2.0*M_PI;

    cosphinew = cos(phi);
    sinphinew = sin(phi);

    *cx = mag*sinthetanew*cosphinew;
    *cy = mag*sinthetanew*sinphinew;
    *cz = mag*costhetanew;

    rotation(costheta, sintheta, cosphi, sinphi, cx, cy, cz);
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void randomSphere(Real *cx, Real *cy, Real *cz, amrex::RandomEngine const& engine)
{
    Real mag, costhetanew, sinthetanew, cosphinew, sinphinew, phi;

    mag = sqrt((*cx)*(*cx) + (*cy)*(*cy) + (*cz)*(*cz));

    costhetanew = 2.0*amrex::Random(engine)-1.0;
    sinthetanew = sqrt(1.0 - costhetanew*costhetanew);

    phi = amrex::Random(engine)*2.0*M_PI;

    cosphinew = cos(phi);
    sinphinew = sin(phi);

    *cx = mag*sinthetanew*cosphinew;
    *cy = mag*sinthetanew*sinphinew;
    *cz = mag*costhetanew;
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void cosineLawHemisphere(Real costheta, Real sintheta, Real cosphi, Real sinphi, Real *cx, Real *cy, Real *cz, Real mag, amrex::RandomEngine const& engine)
{
    Real costhetanew, sinthetanew, cosphinew, sinphinew, phi;

    //    call get_half_angles(costhetanew, sinthetanew, cosphinew, sinphinew)

    sinthetanew = sqrt(amrex::Random(engine));
    costhetanew = sqrt(1.0 - sinthetanew*sinthetanew);

    phi = amrex::Random(engine)*2.0*M_PI;

    cosphinew = cos(phi);
    sinphinew = sin(phi);

    *cx = mag*sinthetanew*cosphinew;
    *cy = mag*sinthetanew*sinphinew;
    *cz = mag*costhetanew;

    rotation(costheta, sintheta, cosphi, sinphi, cx, cy, cz);
}


AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real pn(int n)
{
    return 90.0/pow(M_PI*n,4);
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real erlang4Dist(int n, Real temp, amrex::RandomEngine const& engine)
{
    return -log(amrex::Random(engine)*amrex::Random(engine)*amrex::Random(engine)*amrex::Random(engine))*k_B*temp/(h_bar*n);
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real plankDist(Real temp, amrex::RandomEngine const& engine)
{
    Real sample = 360*k_B*temp*1.037/(h_bar*pow(M_PI,4));//Mean value incase sampling fails

    Real pnTotal = 0;
    Real test = amrex::Random(engine);

    for(int n=1;n<=50;n++)
    {
        pnTotal = pnTotal + pn(n);
        if(test < pnTotal)
        {
            sample = erlang4Dist(n, temp, engine);
            break;
        }
    }

    return sample;
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real truncatedFluxDist(Real vel, Real trunc, Real srt, amrex::RandomEngine const& engine)
{
    vel = -10e10;
    while(vel < trunc)
    {
        if(amrex::Random(engine) < 0.5)
        {
            vel = sqrt(2)*srt*sqrt(-log(amrex::Random()));

        }else
        {
            vel = -sqrt(2)*srt*sqrt(-log(amrex::Random()));
        }
    }

    return vel;
}


AMREX_GPU_HOST_DEVICE AMREX_INLINE
void reduceMassFlux(paramPlane* paramPlaneList, const int paramplanes)
{
    for(int i=0; i<paramplanes; i++)
    {
        for(int l=0;l<nspecies;l++)
        {
            ParallelDescriptor::ReduceRealSum(paramPlaneList[i].massFluxRight[l]);
            ParallelDescriptor::ReduceRealSum(paramPlaneList[i].massFluxLeft[l]);
            //Print() << "Reduced mass flux right: " << paramPlaneList[i].massFluxRight[l] << endl;
            //Print() << "Reduced mass flux left: " << paramPlaneList[i].massFluxLeft[l] << endl;
        }

    }
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void zeroMassFlux(paramPlane* paramPlaneList, const int paramplanes)
{
    for(int i=0; i<paramplanes; i++)
    {
        for(int l=0;l<nspecies;l++)
        {
            paramPlaneList[i].massFluxRight[l]  = 0;
            paramPlaneList[i].massFluxLeft[l] = 0;
            paramPlaneList[i].xMomFluxRight[l]  = 0;
            paramPlaneList[i].yMomFluxRight[l]  = 0;
            paramPlaneList[i].zMomFluxRight[l]  = 0;
            paramPlaneList[i].xMomFluxLeft[l]  = 0;
            paramPlaneList[i].yMomFluxLeft[l]  = 0;
            paramPlaneList[i].zMomFluxLeft[l]  = 0;

            paramPlaneList[i].energyFluxRight[l]  = 0;
            paramPlaneList[i].energyFluxLeft[l] = 0;

        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void app_bc_gpu(paramPlane* surf, FhdParticleContainer::ParticleType& part, int intside, Real* domsize, int *push, Real *runtime, Real inttime, amrex::RandomEngine const& engine)
{

    //BL_PROFILE_VAR("app_bc_gpu()", app_bc_gpu);
    if(intside == 1)
    {
        if(amrex::Random(engine) < surf->sinkRight)
        {
            *runtime = 0;

            part.id() = -1;

            int spec = part.idata(FHD_intData::species);
            surf->massFluxRight[spec] += part.rdata(FHD_realData::mass);
            surf->xMomFluxRight[spec] += part.rdata(FHD_realData::mass)*part.rdata(FHD_realData::velx);
            surf->yMomFluxRight[spec] += part.rdata(FHD_realData::mass)*part.rdata(FHD_realData::vely);
            surf->zMomFluxRight[spec] += part.rdata(FHD_realData::mass)*part.rdata(FHD_realData::velz);

        }
        else if(amrex::Random(engine) > surf->porosityRight)
        {
            *push = 0;

            if(amrex::Random(engine) < surf->specularityRight)
            {

                if(surf->species[0] > -1)
                {
                    Real rand = amrex::Random(engine);
                    //Real massOld = part.idata(FHD_realData::mass);
                    //cout << "Old: " << part.rdata(FHD_realData::velx) << ", " << part.rdata(FHD_realData::vely) << ", " << part.rdata(FHD_realData::velz) << endl;
                    if(rand < surf->species[0])
                    {
                        part.idata(FHD_intData::newSpecies) = 0;
                    }else
                    {
                        part.idata(FHD_intData::newSpecies) = 1;
                    }

/*               Real dotprod = part.rdata(FHD_realData::velx)*surf->rnx + part.rdata(FHD_realData::vely)*surf->rny*/
/*                  + part.rdata(FHD_realData::velz)*surf->rnz;*/
/*                  part.rdata(FHD_realData::velx) = -2.0*dotprod*surf->rnx + part.rdata(FHD_realData::velx);*/
/*                  part.rdata(FHD_realData::vely) = -2.0*dotprod*surf->rny + part.rdata(FHD_realData::vely);*/
/*                  part.rdata(FHD_realData::velz) = -2.0*dotprod*surf->rnz + part.rdata(FHD_realData::velz);     */

/*               part.rdata(FHD_realData::mass) = mass[part.idata(FHD_intData::newSpecies)];*/
/*               part.rdata(FHD_realData::R) = k_B/mass[part.idata(FHD_intData::newSpecies)];*/
/*               part.rdata(FHD_realData::radius) = diameter[part.idata(FHD_intData::newSpecies)]/2.0;    */

                    Real srt = sqrt(part.rdata(FHD_realData::R)*surf->temperatureRight);
                    part.rdata(FHD_realData::velx) = srt*amrex::RandomNormal(0.,1.,engine);
                    part.rdata(FHD_realData::vely) = srt*amrex::RandomNormal(0.,1.,engine);
                    part.rdata(FHD_realData::velz) = 1.414213562*srt*sqrt(-log(amrex::Random(engine)));

                    rotation(surf->cosThetaRight, surf->sinThetaRight,
                             surf->cosPhiRight, surf->sinPhiRight, &part.rdata(FHD_realData::velx), &part.rdata(FHD_realData::vely), &part.rdata(FHD_realData::velz));


/*               if(part.idata(FHD_intData::newSpecies) == part.idata(FHD_intData::species))*/
/*               {*/
/*                  Real dotprod = part.rdata(FHD_realData::velx)*surf->rnx + part.rdata(FHD_realData::vely)*surf->rny*/
/*                      + part.rdata(FHD_realData::velz)*surf->rnz;*/

/*                  part.rdata(FHD_realData::velx) = -2.0*dotprod*surf->rnx + part.rdata(FHD_realData::velx);*/
/*                  part.rdata(FHD_realData::vely) = -2.0*dotprod*surf->rny + part.rdata(FHD_realData::vely);*/
/*                  part.rdata(FHD_realData::velz) = -2.0*dotprod*surf->rnz + part.rdata(FHD_realData::velz);*/
/*               */
/*               }else*/
/*               {*/
/*                  */
/*                  Real velx = -part.rdata(FHD_realData::velx)*massOld/massNew;*/
/*                  Real vely, velz;*/
/*                  if(amrex::Random(engine) < 0.5)*/
/*                  {*/
/*                    vely = part.rdata(FHD_realData::vely)*massOld/massNew;*/
/*                    velz = sqrt(-massOld*massOld*part.rdata(FHD_realData::velx)*part.rdata(FHD_realData::velx) + massOld*massNew*part.rdata(FHD_realData::velx)*part.rdata(FHD_realData::velx) - massOld*massOld*part.rdata(FHD_realData::vely)*part.rdata(FHD_realData::vely)*/
/*                                 + massOld*massNew*part.rdata(FHD_realData::vely)*part.rdata(FHD_realData::vely)  +  massOld*massNew*part.rdata(FHD_realData::velz)*part.rdata(FHD_realData::velz)  )/massNew;*/
/*                    if(amrex::Random(engine) < 0.5)*/
/*                    {*/
/*                        velz = -velz;*/
/*                    }*/
/*                  }else*/
/*                  {*/
/*                    velz = part.rdata(FHD_realData::velz)*massOld/massNew;*/
/*                    vely = sqrt(-massOld*massOld*part.rdata(FHD_realData::velx)*part.rdata(FHD_realData::velx) + massOld*massNew*part.rdata(FHD_realData::velx)*part.rdata(FHD_realData::velx) - massOld*massOld*part.rdata(FHD_realData::velz)*part.rdata(FHD_realData::velz)*/
/*                                 + massOld*massNew*part.rdata(FHD_realData::velz)*part.rdata(FHD_realData::velz)  +  massOld*massNew*part.rdata(FHD_realData::vely)*part.rdata(FHD_realData::vely)  )/massNew;*/
/*                    if(amrex::Random(engine) < 0.5)*/
/*                    {*/
/*                        vely = -vely;*/
/*                    }*/
/*                  */
/*                  }*/
/*                  */
/*                  part.rdata(FHD_realData::velx) = velx;*/
/*                  part.rdata(FHD_realData::vely) = vely;*/
/*                  part.rdata(FHD_realData::velz) = velz;*/
/*            */
/*               */
/*               }*/
/*               cout << "New: " << part.rdata(FHD_realData::velx) << ", " << part.rdata(FHD_realData::vely) << ", " << part.rdata(FHD_realData::velz) << endl;*/

                }else
                {

                    Real dotprod = part.rdata(FHD_realData::velx)*surf->rnx + part.rdata(FHD_realData::vely)*surf->rny
                                 + part.rdata(FHD_realData::velz)*surf->rnz;
                    part.rdata(FHD_realData::velx) = -2.0*dotprod*surf->rnx + part.rdata(FHD_realData::velx);
                    part.rdata(FHD_realData::vely) = -2.0*dotprod*surf->rny + part.rdata(FHD_realData::vely);
                    part.rdata(FHD_realData::velz) = -2.0*dotprod*surf->rnz + part.rdata(FHD_realData::velz);

                }

            }
            else
            {
                Real srt = sqrt(part.rdata(FHD_realData::R)*surf->temperatureRight);
                part.rdata(FHD_realData::velx) = srt*amrex::RandomNormal(0.,1.,engine);
                part.rdata(FHD_realData::vely) = srt*amrex::RandomNormal(0.,1.,engine);
                part.rdata(FHD_realData::velz) = 1.414213562*srt*sqrt(-log(amrex::Random(engine)));

                rotation(surf->cosThetaRight, surf->sinThetaRight,
                         surf->cosPhiRight, surf->sinPhiRight,
                         &part.rdata(FHD_realData::velx), &part.rdata(FHD_realData::vely), &part.rdata(FHD_realData::velz));

            }
        }
        else if(amrex::Random(engine) < surf->periodicity)
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

        }
        else
        {
            *push = 1;
            if(amrex::Random(engine) > surf->momentumConsRight)
            {
                randomhemisphere(surf->cosThetaRight, surf->sinThetaRight,
                                 surf->cosPhiRight, surf->sinPhiRight,
                                 &part.rdata(FHD_realData::velx),&part.rdata(FHD_realData::vely), &part.rdata(FHD_realData::velz),
                                 engine);
            }
        }

    }
    // Left boundary
    else
    {

        if(amrex::Random(engine) < surf->sinkLeft)
        {
            *runtime = 0;
            part.id() = -1;
            int spec = part.idata(FHD_intData::species);

            surf->massFluxLeft[spec] += part.rdata(FHD_realData::mass);
            surf->xMomFluxLeft[spec] += part.rdata(FHD_realData::mass)*part.rdata(FHD_realData::velx);
            surf->yMomFluxLeft[spec] += part.rdata(FHD_realData::mass)*part.rdata(FHD_realData::vely);
            surf->zMomFluxLeft[spec] += part.rdata(FHD_realData::mass)*part.rdata(FHD_realData::velz);

        }
        else if(amrex::Random(engine) > surf->porosityLeft)
        {
            *push = 0;

            if(amrex::Random(engine) < surf->specularityLeft)
            {

                if(surf->species[0] > -1)
                {
                    Real rand = amrex::Random(engine);
/*                  Real massOld = part.idata(FHD_realData::mass);*/
                    if(rand < surf->species[0])
                    {
                        part.idata(FHD_intData::newSpecies) = 0;
                    }else
                    {
                        part.idata(FHD_intData::newSpecies) = 1;
                    }

/*               Real dotprod = part.rdata(FHD_realData::velx)*surf->lnx + part.rdata(FHD_realData::vely)*surf->lny*/
/*                      + part.rdata(FHD_realData::velz)*surf->lnz;*/

/*               part.rdata(FHD_realData::velx) = -2.0*dotprod*surf->lnx + part.rdata(FHD_realData::velx);*/
/*               part.rdata(FHD_realData::vely) = -2.0*dotprod*surf->lny + part.rdata(FHD_realData::vely);*/
/*               part.rdata(FHD_realData::velz) = -2.0*dotprod*surf->lnz + part.rdata(FHD_realData::velz);*/

/*               Real massNew = mass[part.idata(FHD_intData::newSpecies)];*/
/*               */
/*               part.rdata(FHD_realData::mass) = mass[part.idata(FHD_intData::newSpecies)];*/
/*               part.rdata(FHD_realData::R) = k_B/mass[part.idata(FHD_intData::newSpecies)];*/
/*               part.rdata(FHD_realData::radius) = diameter[part.idata(FHD_intData::newSpecies)]/2.0;    */

                    Real srt = sqrt(part.rdata(FHD_realData::R)*surf->temperatureLeft);
                    part.rdata(FHD_realData::velx) = srt*amrex::RandomNormal(0.,1.,engine);
                    part.rdata(FHD_realData::vely) = srt*amrex::RandomNormal(0.,1.,engine);
                    part.rdata(FHD_realData::velz) = 1.414213562*srt*sqrt(-log(amrex::Random(engine)));

                    rotation(surf->cosThetaLeft, surf->sinThetaLeft,
                             surf->cosPhiLeft, surf->sinPhiLeft, &part.rdata(FHD_realData::velx), &part.rdata(FHD_realData::vely), &part.rdata(FHD_realData::velz));

/*               if(part.idata(FHD_intData::newSpecies) == part.idata(FHD_intData::species))*/
/*               {*/
/*                  Real dotprod = part.rdata(FHD_realData::velx)*surf->lnx + part.rdata(FHD_realData::vely)*surf->lny*/
/*                      + part.rdata(FHD_realData::velz)*surf->lnz;*/

/*                  part.rdata(FHD_realData::velx) = -2.0*dotprod*surf->lnx + part.rdata(FHD_realData::velx);*/
/*                  part.rdata(FHD_realData::vely) = -2.0*dotprod*surf->lny + part.rdata(FHD_realData::vely);*/
/*                  part.rdata(FHD_realData::velz) = -2.0*dotprod*surf->lnz + part.rdata(FHD_realData::velz);*/
/*               */
/*               }else*/
/*               {*/
/*                  */
/*                  Real velx = -part.rdata(FHD_realData::velx)*massOld/massNew;*/
/*                  Real vely, velz;*/
/*                  if(amrex::Random(engine) < 0.5)*/
/*                  {*/
/*                    vely = part.rdata(FHD_realData::vely)*massOld/massNew;*/
/*                    velz = sqrt(-massOld*massOld*part.rdata(FHD_realData::velx)*part.rdata(FHD_realData::velx) + massOld*massNew*part.rdata(FHD_realData::velx)*part.rdata(FHD_realData::velx) - massOld*massOld*part.rdata(FHD_realData::vely)*part.rdata(FHD_realData::vely)*/
/*                                 + massOld*massNew*part.rdata(FHD_realData::vely)*part.rdata(FHD_realData::vely)  +  massOld*massNew*part.rdata(FHD_realData::velz)*part.rdata(FHD_realData::velz)  )/massNew;*/
/*                    if(amrex::Random(engine) < 0.5)*/
/*                    {*/
/*                        velz = -velz;*/
/*                    }*/
/*                  }else*/
/*                  {*/
/*                    velz = part.rdata(FHD_realData::velz)*massOld/massNew;*/
/*                    vely = sqrt(-massOld*massOld*part.rdata(FHD_realData::velx)*part.rdata(FHD_realData::velx) + massOld*massNew*part.rdata(FHD_realData::velx)*part.rdata(FHD_realData::velx) - massOld*massOld*part.rdata(FHD_realData::velz)*part.rdata(FHD_realData::velz)*/
/*                                 + massOld*massNew*part.rdata(FHD_realData::velz)*part.rdata(FHD_realData::velz)  +  massOld*massNew*part.rdata(FHD_realData::vely)*part.rdata(FHD_realData::vely)  )/massNew;*/
/*                    if(amrex::Random(engine) < 0.5)*/
/*                    {*/
/*                        vely = -vely;*/
/*                    }*/
/*                  */
/*                  }*/
/*                  */
/*                  part.rdata(FHD_realData::velx) = velx;*/
/*                  part.rdata(FHD_realData::vely) = vely;*/
/*                  part.rdata(FHD_realData::velz) = velz;*/
/*                  part.rdata(FHD_realData::mass) = massNew;*/
/*                  part.rdata(FHD_realData::R) = k_B/massNew;*/
/*                  part.rdata(FHD_realData::radius) = diameter[part.idata(FHD_intData::newSpecies)]/2.0;                */
/*               */
/*               }*/



                }else
                {
                    Real dotprod = part.rdata(FHD_realData::velx)*surf->lnx + part.rdata(FHD_realData::vely)*surf->lny
                                 + part.rdata(FHD_realData::velz)*surf->lnz;

                    part.rdata(FHD_realData::velx) = -2.0*dotprod*surf->lnx + part.rdata(FHD_realData::velx);
                    part.rdata(FHD_realData::vely) = -2.0*dotprod*surf->lny + part.rdata(FHD_realData::vely);
                    part.rdata(FHD_realData::velz) = -2.0*dotprod*surf->lnz + part.rdata(FHD_realData::velz);
                }
            }
            else
            {
                Real srt = sqrt(part.rdata(FHD_realData::R)*surf->temperatureLeft);
                part.rdata(FHD_realData::velx) = srt*amrex::RandomNormal(0.,1.,engine);
                part.rdata(FHD_realData::vely) = srt*amrex::RandomNormal(0.,1.,engine);
                part.rdata(FHD_realData::velz) = sqrt(2)*srt*sqrt(-log(amrex::Random(engine)));

                rotation(surf->cosThetaLeft, surf->sinThetaLeft,
                         surf->cosPhiLeft, surf->sinPhiLeft,
                         &part.rdata(FHD_realData::velx), &part.rdata(FHD_realData::vely), &part.rdata(FHD_realData::velz));
            }

        }
        else if(amrex::Random(engine) < surf->periodicity)
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

        }
        else
        {
            *push = 1;
            if(amrex::Random(engine) > surf->momentumConsLeft)
            {
                randomhemisphere(surf->cosThetaRight, surf->sinThetaRight,
                                 surf->cosPhiRight, surf->sinPhiRight,
                                 &part.rdata(FHD_realData::velx), &part.rdata(FHD_realData::vely), &part.rdata(FHD_realData::velz),
                                 engine);
            }
        }
    }
//   if(graphene_tog .eq. 1) then
//      if(surf%boundary .eq. 6) then

//      !print *, "new: ", part%vel(3)
//      call surf_velocity(surf, part, time, oldvel, inttime)
//   endif
//   endif

}


#endif
