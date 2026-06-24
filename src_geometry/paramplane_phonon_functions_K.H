#ifndef _paramplane_phonon_functions_K_H_
#define _paramplane_phonon_functions_K_H_

#include <AMReX.H>
#include <common_namespace.H>
#include <math.h>

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void app_bc_phonon_gpu(const paramPlane* surf, FhdParticleContainer::ParticleType& part, int intside, Real* domsize, int *push, Real *runtime, const int step, int *count, int *specCount, amrex::RandomEngine const& engine)
{

    if(intside == 1)
    {
        if(amrex::Random(engine) < surf->momentumConsRight)
        {
            //           std::string plotfilename = std::to_string(surf->boundary) + "_" + std::to_string(ParallelDescriptor::MyProc()) + amrex::Concatenate("_particles_right_",step,12);
            //           string plotfilename = "test";
            /*           std::ofstream ofs(plotfilename, std::ios::app);*/
            /*           ofs << part.pos(0) << " " << part.pos(1) << " " << part.pos(2) << " " << part.rdata(FHD_realData::velx) << " " << part.rdata(FHD_realData::vely) << " " << part.rdata(FHD_realData::velz) << " " << part.rdata(FHD_realData::omega) << std::endl;*/
            /*           ofs.close();*/
            part.idata(FHD_intData::fluxRec) = surf->boundary;
            //         char str[] = std::to_string(surf->boundary);
        }

        if(amrex::Random(engine) < surf->sinkRight)
        {
            *runtime = 0;

            part.id() = -1;
        }
        else if(amrex::Random(engine) > surf->porosityRight)
        {
            *push = 0;

            //(*count)++;
            amrex::Gpu::Atomic::Add(count, 1);
            Real dotprod = part.rdata(FHD_realData::velx)*surf->rnx + part.rdata(FHD_realData::vely)*surf->rny + part.rdata(FHD_realData::velz)*surf->rnz;
            Real speed = sqrt(part.rdata(FHD_realData::velx)*part.rdata(FHD_realData::velx)+part.rdata(FHD_realData::vely)*part.rdata(FHD_realData::vely)+part.rdata(FHD_realData::velz)*part.rdata(FHD_realData::velz));
            Real cosAlpha = dotprod/speed;

            Real probSpec = exp(-16*M_PI*M_PI*surf->specularityRight*surf->specularityRight*cosAlpha*cosAlpha/(part.rdata(FHD_realData::lambda)*part.rdata(FHD_realData::lambda)));

            //cout << "SPEC " << surf->boundary << ": " << part.id() << ", " << surf->specularityRight << ", " << part.rdata(FHD_realData::lambda) << ", " << probSpec << endl;

            if(amrex::Random(engine) < probSpec)
            {
                //Specular
                //(*specCount)++;
                amrex::Gpu::Atomic::Add(specCount, 1);
                dotprod = part.rdata(FHD_realData::velx)*surf->rnx + part.rdata(FHD_realData::vely)*surf->rny
                         + part.rdata(FHD_realData::velz)*surf->rnz;

                part.rdata(FHD_realData::velx) = -2.0*dotprod*surf->rnx + part.rdata(FHD_realData::velx);
                part.rdata(FHD_realData::vely) = -2.0*dotprod*surf->rny + part.rdata(FHD_realData::vely);
                part.rdata(FHD_realData::velz) = -2.0*dotprod*surf->rnz + part.rdata(FHD_realData::velz);

            }
            else
            {
                //Lambert cosine law

                /*          randomhemisphere(surf->cosThetaRight, surf->sinThetaRight, surf->cosPhiRight, surf->sinPhiRight,*/
                /*                            &part.rdata(FHD_realData::velx),&part.rdata(FHD_realData::vely), &part.rdata(FHD_realData::velz), */
                /*                            engine);        */

                cosineLawHemisphere(surf->cosThetaRight, surf->sinThetaRight, surf->cosPhiRight, surf->sinPhiRight,
                                    &part.rdata(FHD_realData::velx),&part.rdata(FHD_realData::vely), &part.rdata(FHD_realData::velz), speed,
                                    engine);

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
        }


    }
    // Left boundary
    else
    {
        if(amrex::Random(engine) < surf->momentumConsLeft)
        {
            /*           std::string plotfilename = std::to_string(surf->boundary) + "_" + std::to_string(ParallelDescriptor::MyProc()) + amrex::Concatenate("_particles_left_",step,12);*/
            /*           std::ofstream ofs(plotfilename, std::ios::app);*/
            /*           ofs << part.pos(0) << " " << part.pos(1) << " " << part.pos(2) << " " << part.rdata(FHD_realData::velx) << " " << part.rdata(FHD_realData::vely) << " " << part.rdata(FHD_realData::velz) << " " << part.rdata(FHD_realData::omega) << std::endl;*/
            /*           ofs.close();*/
            // cout << "BC: " << surf->boundary << endl;
            part.idata(FHD_intData::fluxRec) = -surf->boundary;
            //   cout << "BC2: " << part.idata(FHD_intData::fluxRec) << endl;
        }

        if(amrex::Random(engine) < surf->sinkLeft)
        {
            *runtime = 0;
            part.id() = -1;
        }
        else if(amrex::Random(engine) > surf->porosityLeft)
        {
            *push = 0;
            //(*count)++;
            amrex::Gpu::Atomic::Add(count, 1);
            Real dotprod = part.rdata(FHD_realData::velx)*surf->lnx + part.rdata(FHD_realData::vely)*surf->lny + part.rdata(FHD_realData::velz)*surf->lnz;
            Real speed = sqrt(part.rdata(FHD_realData::velx)*part.rdata(FHD_realData::velx)+part.rdata(FHD_realData::vely)*part.rdata(FHD_realData::vely)+part.rdata(FHD_realData::velz)*part.rdata(FHD_realData::velz));
            Real cosAlpha = dotprod/speed;

            Real probSpec = exp(-16*M_PI*M_PI*surf->specularityLeft*surf->specularityLeft*cosAlpha*cosAlpha/(part.rdata(FHD_realData::lambda)*part.rdata(FHD_realData::lambda)));

            if(amrex::Random(engine) < probSpec)
            {
                //Specular
                //(*specCount)++;
                amrex::Gpu::Atomic::Add(specCount, 1);
                dotprod = part.rdata(FHD_realData::velx)*surf->lnx + part.rdata(FHD_realData::vely)*surf->lny
                         + part.rdata(FHD_realData::velz)*surf->lnz;

                part.rdata(FHD_realData::velx) = -2.0*dotprod*surf->lnx + part.rdata(FHD_realData::velx);
                part.rdata(FHD_realData::vely) = -2.0*dotprod*surf->lny + part.rdata(FHD_realData::vely);
                part.rdata(FHD_realData::velz) = -2.0*dotprod*surf->lnz + part.rdata(FHD_realData::velz);

            }
            else
            {
                //Lambert cosine law
                cosineLawHemisphere(surf->cosThetaLeft, surf->sinThetaLeft, surf->cosPhiLeft, surf->sinPhiLeft,
                                    &part.rdata(FHD_realData::velx),&part.rdata(FHD_realData::vely), &part.rdata(FHD_realData::velz), speed,
                                    engine);

                /*            randomhemisphere(surf->cosThetaLeft, surf->sinThetaLeft, surf->cosPhiLeft, surf->sinPhiLeft,*/
                /*                &part.rdata(FHD_realData::velx),&part.rdata(FHD_realData::vely), &part.rdata(FHD_realData::velz), */
                /*                engine);       */
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
        }
    }


}


#endif
