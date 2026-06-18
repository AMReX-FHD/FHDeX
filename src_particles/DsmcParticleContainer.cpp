#include "DsmcParticleContainer.H"

// #include "particle_functions_K.H"
#include "paramplane_functions_K.H"
#include "paramplane_phonon_functions_K.H"
#include <math.h>
using namespace std;
FhdParticleContainer::FhdParticleContainer(const Geometry & geom, const DistributionMapping & dmap,
    const BoxArray & ba, int ncells)
    : NeighborParticleContainer<FHD_realData::count, FHD_intData::count> (geom, dmap, ba, ncells)
{
    BL_PROFILE_VAR("FhdParticleContainer()",FhdParticleContainer);

    realParticles = 0;
    simParticles = 0;

    totalCollisionCells = n_cells[0]*n_cells[1]*n_cells[2];
    domainVol = (prob_hi[0] - prob_lo[0])*(prob_hi[1] - prob_lo[1])*(prob_hi[2] - prob_lo[2]);

    collisionCellVol = domainVol/totalCollisionCells;
    ocollisionCellVol = 1.0/collisionCellVol;

    bool rho_defined = true;
    if(rho0<0)
    {
        rho0=0.;
        rho_defined = false;
    }

    // total = simulated total (total*particle_neff = real number of particles)
    // n0 = real number density
    for(int i=0;i<nspecies;i++)
    {
        properties[i].mass = mass[i];
        properties[i].radius = diameter[i]/2.0;
        properties[i].partVol = pow(diameter[i],3.0)*pi_usr/6.0;
        properties[i].part2cellVol = properties[i].partVol*ocollisionCellVol;
        properties[i].Neff = particle_neff; //   real / simulated
        properties[i].R = k_B/properties[i].mass;

        if(particle_count[i]>=0)
        {
            //Print() << "pcnt\n";
            properties[i].total = particle_count[i];
            properties[i].n0 = particle_neff*properties[i].total/domainVol;
            particle_n0[i] = properties[i].n0;

            amrex::Print() <<  "Species "<< i << " count " << properties[i].total << "\n";
            amrex::Print() <<  "Species "<< i << " n0 " << properties[i].n0 << "\n";
            amrex::Print() <<  "Species " << i << " rho0 " << properties[i].n0*properties[i].mass << "\n";

            rho0 += properties[i].n0*properties[i].mass;
            Yk0[i] = properties[i].n0*properties[i].mass;
        }
        else if(phi_domain[i]>=0)
        {
            //Print() << "phi\n";
            properties[i].total = (int)amrex::Math::ceil(
                (phi_domain[i]*domainVol)/(properties[i].partVol*particle_neff) );
            properties[i].n0 = particle_neff*properties[i].total/domainVol;
            particle_n0[i] = properties[i].n0;

            amrex::Print() <<  "Species "<< i << " count " << properties[i].total << "\n";
            amrex::Print() <<  "Species "<< i << " n0 " << properties[i].n0 << "\n";
            amrex::Print() <<  "Species " << i << " rho0 " << properties[i].n0*properties[i].mass << "\n";

            rho0 += properties[i].n0*properties[i].mass;
            Yk0[i] = properties[i].n0*properties[i].mass;
        }
        else if(rho_defined > 0)
        {
            //Print() << "rho\n";
            Real Yktot = 0.0;
            for(int j=0; j<nspecies; j++) { Yktot += Yk0[j]; }
            for(int j=0; j<nspecies; j++) {    Yk0[j] /= Yktot; }

            Real rhop = properties[i].mass/domainVol;
            amrex::Print() <<  "Species "<< i << " rhop " << rhop << "\n";
            amrex::Print() <<  "Species "<< i << " Yk0 " << Yk0[i] << "\n";
            amrex::Print() <<  "Species "<< i << " rho0 " << rho0 << "\n";
            amrex::Print() <<  "Species "<< i << " neff " << particle_neff << "\n";
            properties[i].total = std::ceil((rho0*Yk0[i])/(rhop*particle_neff));
            amrex::Print() <<  "Species "<< i << " total " << properties[i].total << "\n";
            properties[i].n0 = particle_neff*properties[i].total/domainVol;
            particle_n0[i] = properties[i].n0;

    //            Real specRho = rho0*Yk0[i];
    //            properties[i].total = std::ceil((specRho/properties[i].mass)*domainVol/particle_neff);
    //            properties[i].n0 = particle_neff*properties[i].total/properties[i].mass;
    //            particle_n0[i] = properties[i].n0;


            amrex::Print() <<  "Species "<< i << " count " << properties[i].total << "\n";
            amrex::Print() <<  "Species "<< i << " n0 " << properties[i].n0 << "\n";
            amrex::Print() <<  "Species " << i << " rho0 " << properties[i].n0*properties[i].mass << ", from " << rho0*Yk0[i] << "\n";
            for(int j=0; j<nspecies; j++)
            {
                Yk0[j] *= rho0;
            }
        }
        else
        {
            //Print() << "n0\n";
            amrex::Print() << "n0: " << particle_n0[i] << "\n";
            properties[i].total = (int)amrex::Math::ceil(particle_n0[i]
                *domainVol/particle_neff);
            properties[i].n0 = particle_neff*properties[i].total/domainVol;
            particle_n0[i] = properties[i].n0;

            amrex::Print() <<  "Species " << i << " count " << properties[i].total << "\n";
            amrex::Print() <<  "Species " << i << " n0 " << properties[i].n0 << "\n";
            amrex::Print() <<  "Species " << i << " rho0 " << properties[i].n0*properties[i].mass << "\n";

            rho0 += properties[i].n0*properties[i].mass;
            Yk0[i] = properties[i].n0*properties[i].mass;
        }
        realParticles = realParticles + properties[i].total*particle_neff;
        simParticles = simParticles + properties[i].total;
    //    amrex::Print() << "Particles per cell for species " << i
    //        << " is " << properties[i].total/totalCollisionCells << "\n";
    }
    amrex::Print() <<  "Rho0: " << rho0 << "\n";
    amrex::Print() <<  "Total n0: " << realParticles/domainVol << "\n";

    for(int i=0; i<nspecies ; i++)
    {
        Yk0[i] = Yk0[i]/rho0;
    }

    for (int d=0; d<AMREX_SPACEDIM; ++d)
    {
        domSize[d] = prob_hi[d] - prob_lo[d];
    }

    //    const int lev=0;
    //    for(MFIter mfi = MakeMFIter(lev, false); mfi.isValid(); ++mfi)
    //    {
    //        const Box& box = mfi.validbox();
    //        const int grid_id = mfi.index();
    //
    //        for(int i=0;i<nspecies;i++)
    //        {
    //            m_cell_vectors[i][grid_id].resize(box.numPts());
    //        }
    //    }

    }

    void FhdParticleContainer::MoveParticlesCPP(const Real dt, paramPlane* paramPlaneList, const int paramPlaneCount)
    {
        BL_PROFILE_VAR("MoveParticlesCPP()", MoveParticlesCPP);

        const int lev = 0;
        const GpuArray<Real, 3> dx = Geom(lev).CellSizeArray();
        const GpuArray<Real, 3> plo = Geom(lev).ProbLoArray();
        const GpuArray<Real, 3> phi = Geom(lev).ProbHiArray();

        int np_tile = 0, np_proc = 0;

        GpuArray<int, MAX_SPECIES> delCountExt;
        for(int i = 0; i<MAX_SPECIES;i++)
        {
            delCountExt[i]=0;
        }
        int* delCount = delCountExt.data();

        Gpu::ManagedVector<Real> domsize(3);
        Real* pdomsize = domsize.data();

        for (int d=0; d<AMREX_SPACEDIM; ++d)
        {
            pdomsize[d] = phi[d]-plo[d];
        }

        Gpu::ManagedVector<paramPlane> paramPlaneListTmp;
        paramPlaneListTmp.resize(paramPlaneCount);
        for(int i=0;i<paramPlaneCount;i++)
        {
            paramPlaneListTmp[i]=paramPlaneList[i];

        }
        paramPlane* paramPlaneListPtr = paramPlaneListTmp.data();
        //paramPlane* paramPlaneListPtr = paramPlaneList;

        int totalParts = 0;
        //amrex::RandomEngine engine;

        for (FhdParIter pti(* this, lev); pti.isValid(); ++pti)
        {
            const int grid_id = pti.index();
            const int tile_id = pti.LocalTileIndex();
            const Box& tile_box  = pti.tilebox();

            auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];
            auto& aos = particle_tile.GetArrayOfStructs();
            ParticleType* particles = aos().dataPtr();
            const long np = particle_tile.numParticles();

            Box bx  = pti.tilebox();
            IntVect myLo = bx.smallEnd();
            IntVect myHi = bx.bigEnd();

            //cout << "Rank " << ParallelDescriptor::MyProc() << " sees " << np << " particles\n";

            totalParts += np;

    //        for (int i = 0; i < np; i++)
            amrex::ParallelForRNG(np, [=] AMREX_GPU_DEVICE (int i, amrex::RandomEngine const& engine) noexcept
            {
                ParticleType & part = particles[i];
                Real runtime = dt*part.rdata(FHD_realData::timeFrac);
                Real inttime;
                int intsurf, intside, push;
    //            Real adj = 0.99999;
    //            Real adjalt = 2.0*(1.0-adj);
                //printf ("moving %d\n", i);
                //paramPlane paramPlaneListT[3];
                //int paramPlaneCountT = 3;

                while(runtime > 0)
                {
                    find_inter_gpu(part, runtime, paramPlaneListPtr, paramPlaneCount,
                        &intsurf, &inttime, &intside, AMREX_ZFILL(plo), AMREX_ZFILL(phi));

                    for (int d=0; d<(AMREX_SPACEDIM); ++d)
                    //for (int d=0; d<(AMREX_SPACEDIM-2); ++d)
                    {
                        part.pos(d) += inttime * part.rdata(FHD_realData::velx + d)*ADJ;
                    }

                    runtime = runtime - inttime;

                    if(intsurf > 0)
                    {
                        //find_inter indexes from 1 to maintain compatablity with fortran version
                        paramPlane& surf = paramPlaneListPtr[intsurf-1];

                        Real posAlt[3];

                        for (int d=0; d<(AMREX_SPACEDIM); ++d)
                        //for (int d=0; d<(AMREX_SPACEDIM-2); ++d)
                        {
                            posAlt[d] = inttime * part.rdata(FHD_realData::velx + d)*ADJALT;
                        }

                        Real dummy = 1;
                        //cout << "particle " << part.id() << " intersected " << intsurf << " with vel " << part.rdata(FHD_realData::velx) << ", " << part.rdata(FHD_realData::vely) << ", " << part.rdata(FHD_realData::velz) << endl;
                        app_bc_gpu(&surf, part, intside, pdomsize, &push, &runtime, dummy, engine);
                        if(part.id() == -1)
                        {
                            delCount[part.idata(FHD_intData::species)]++;
                        }
                        if(push == 1)
                        {
                            //for (int d=0; d<(AMREX_SPACEDIM-2); ++d)
                            for (int d=0; d<(AMREX_SPACEDIM); ++d)
                            {
                                part.pos(d) += part.pos(d) + posAlt[d];
                            }
                        }
                    }

                }

                part.idata(FHD_intData::sorted) = -1;

                part.rdata(FHD_realData::timeFrac) = 1;

    //            IntVect iv(part.idata(FHD_intData::i), part.idata(FHD_intData::j), part.idata(FHD_intData::k));
    //            long imap = tile_box.index(iv);
    //            m_cell_vectors[part.idata(FHD_intData::species)][grid_id][imap].pop_back();

    //            int cell[3];
    //            cell[0] = (int)floor((part.pos(0)-plo[0])/dx[0]);
    //            cell[1] = (int)floor((part.pos(1)-plo[1])/dx[1]);
    //            cell[2] = (int)floor((part.pos(2)-plo[2])/dx[2]);

    //            if((part.idata(FHD_intData::i) != cell[0]) || (part.idata(FHD_intData::j) != cell[1]) ||
    //                (part.idata(FHD_intData::k) != cell[2]) || part.id() < 0 || part.idata(FHD_intData::newSpecies) != -1)
    //            {
    //                IntVect iv(part.idata(FHD_intData::i), part.idata(FHD_intData::j), part.idata(FHD_intData::k));
    //                long imap = tile_box.index(iv);

    //                int lastIndex = m_cell_vectors[part.idata(FHD_intData::species)][pti.index()][imap].size() - 1;
    //                int lastPart = m_cell_vectors[part.idata(FHD_intData::species)][pti.index()][imap][lastIndex];
    //                int newIndex = part.idata(FHD_intData::sorted);

    //                m_cell_vectors[part.idata(FHD_intData::species)][pti.index()][imap][newIndex] = lastPart;
    //                m_cell_vectors[part.idata(FHD_intData::species)][pti.index()][imap].pop_back();

    //                particles[lastPart].idata(FHD_intData::sorted) = newIndex;

    //                part.idata(FHD_intData::sorted) = -1;
    //                if(part.idata(FHD_intData::newSpecies) != -1)
    //                {
    //                    part.idata(FHD_intData::species) = part.idata(FHD_intData::newSpecies);
    //                    part.idata(FHD_intData::newSpecies) = -1;
    //                }
    //            }

                if(part.idata(FHD_intData::newSpecies) != -1)
                {
                    part.idata(FHD_intData::species) = part.idata(FHD_intData::newSpecies);
                    part.idata(FHD_intData::newSpecies) = -1;
                }

            });
        }

    //    for(int i = 0; i<MAX_SPECIES;i++)
    //    {
    //        int temp = delCount[i];
    //        ParallelDescriptor::ReduceIntSum(temp);
    //        if(temp != 0)
    //        {
    //            Print() << "Deleted " << temp << " of species " << i << "\n";
    //        }
    //    }

        ParallelDescriptor::ReduceIntSum(totalParts);
        simParticles = totalParts;
        //Print() << "Total particles: " << totalParts << "\n";
        Redistribute();
        //SortParticles();
        SortParticlesDB();
    }

    void FhdParticleContainer::MovePhononsCPP(const Real dt, paramPlane* paramPlaneList, const int paramPlaneCount, const int step, const int istep, iMultiFab& bCell)
    {
        BL_PROFILE_VAR("MoveParticlesCPP()", MoveParticlesCPP);

        const int lev = 0;
        const GpuArray<Real, 3> dx = Geom(lev).CellSizeArray();
        const GpuArray<Real, 3> plo = Geom(lev).ProbLoArray();
        const GpuArray<Real, 3> phi = Geom(lev).ProbHiArray();

        int np_tile = 0, np_proc = 0, scatterCount = 0, count = 0, specCount = 0;

        GpuArray<int, MAX_SPECIES> delCountExt;
        for(int i = 0; i<MAX_SPECIES;i++)
        {
            delCountExt[i]=0;
        }
        int* delCount = delCountExt.data();

        Gpu::ManagedVector<Real> domsize(3);
        Real* pdomsize = domsize.data();

        for (int d=0; d<AMREX_SPACEDIM; ++d)
        {
            pdomsize[d] = phi[d]-plo[d];
        }

        Gpu::ManagedVector<paramPlane> paramPlaneListTmp;
        paramPlaneListTmp.resize(paramPlaneCount);
        for(int i=0;i<paramPlaneCount;i++)
        {
            paramPlaneListTmp[i]=paramPlaneList[i];

        }
        paramPlane* paramPlaneListPtr = paramPlaneListTmp.data();
        //paramPlane* paramPlaneListPtr = &paramPlaneList[0];

        int totalParts = 0;

        //int* delCountPtr = delCount;

        Gpu::ManagedVector<int> scatterCountVec(1);
        int * scatterCountPtr = scatterCountVec.data();
        scatterCountVec[0] = 0;

        Gpu::ManagedVector<int> specCountVec(1);
        int * specCountPtr = specCountVec.data();
        specCountVec[0] = 0;

        Gpu::ManagedVector<int> countVec(1);
        int * countPtr = countVec.data();
        countVec[0] = 0;

        //cout << "Rank " << ParallelDescriptor::MyProc() << " start move\n";

        for (FhdParIter pti(* this, lev); pti.isValid(); ++pti)
        {
            const int grid_id = pti.index();
            const int tile_id = pti.LocalTileIndex();
            const Box& tile_box  = pti.tilebox();

            auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];
            auto& aos = particle_tile.GetArrayOfStructs();
            ParticleType* particles = aos().dataPtr();
            const long np = particle_tile.numParticles();

            Array4<int> bCellArr  = bCell[pti].array();

            Box bx  = pti.tilebox();
            IntVect myLo = bx.smallEnd();
            IntVect myHi = bx.bigEnd();

            Real tau_i_p = tau_i;
            Real tau_ta_p = tau_ta;
            Real tau_la_p = tau_la;

            //cout << "Rank " << ParallelDescriptor::MyProc() << " sees " << np << " particles. dt: " << dt << endl;

            totalParts += np;

            //for (int i = 0; i < np; i++)
            //{
            amrex::ParallelForRNG(np, [=] AMREX_GPU_DEVICE (int i, amrex::RandomEngine const& engine) noexcept
            {
                ParticleType & part = particles[i];
                Real runtime = dt*part.rdata(FHD_realData::timeFrac);
                int intsurf, intside, push;
                Real inttime = dt;


                part.idata(FHD_intData::fluxRec) = 0;


                while(runtime > 0)
                {
                    //Print() << "Pre " << part.id() << ": " << part.rdata(FHD_realData::velx + 0) << ", " << part.rdata(FHD_realData::velx + 1) << ", " << part.rdata(FHD_realData::velx + 2) << endl;
    //                printf("DT: %e\n", dt);
    //                cout << "DT: " << dt << endl;

                    int cell[AMREX_SPACEDIM];
                    for (int d=0; d<AMREX_SPACEDIM; ++d)
                    {
                        cell[d] = (int)floor((part.pos(d)-prob_lo[d])/dx[d]);
                    }

                    if(bCellArr(cell[0],cell[1],cell[2]) != 1)
                    {

                        find_inter_gpu(part, runtime, paramPlaneListPtr, paramPlaneCount,
                            &intsurf, &inttime, &intside, AMREX_ZFILL(plo), AMREX_ZFILL(phi));
                    }


                    Real tauImpurityInv = pow(part.rdata(FHD_realData::omega),4)/tau_i_p;
                    Real tauTAInv = part.rdata(FHD_realData::omega)*pow(T_init[0],4)/tau_ta_p;
                    Real tauLAInv = pow(part.rdata(FHD_realData::omega),2)*pow(T_init[0],3)/tau_la_p;
                    Real tauNormalInv = (2.0*tauTAInv+tauLAInv)/3.0;
                    Real tauInv = tauImpurityInv + tauNormalInv;

                    Real scatterTime = -log(amrex::Random(engine))/tauInv;

                    if(scatterTime > inttime)
                    {

                      for (int d=0; d<AMREX_SPACEDIM; ++d)
                      {
                            part.pos(d) += inttime * part.rdata(FHD_realData::velx + d)*ADJ;
                      }

                        runtime = runtime - inttime;
                        part.rdata(FHD_realData::travelTime) += inttime;

                        if(intsurf > 0)
                        {
                            //find_inter indexes from 1 to maintain compatablity with fortran version
                            const paramPlane& surf = paramPlaneListPtr[intsurf-1];
                            //Print() << "Hitting " << intsurf << endl;
                            Real posAlt[3];

                            for (int d=0; d<AMREX_SPACEDIM; ++d)
                            {
                                posAlt[d] = inttime * part.rdata(FHD_realData::velx + d)*ADJALT;
                            }

                            app_bc_phonon_gpu(&surf, part, intside, pdomsize, &push, &runtime, step, countPtr, specCountPtr, engine);
                            //app_bc_gpu(&surf, part, intside, pdomsize, &push, &runtime, dummy, engine);
                            //Print() << "Post " << part.id() << ": " << part.rdata(FHD_realData::velx + 0) << ", " << part.rdata(FHD_realData::velx + 1) << ", " << part.rdata(FHD_realData::velx + 2) << endl;
                            if(part.id() == -1)
                            {
                                //delCountPtr[part.idata(FHD_intData::species)]++;
                            }
                            if(push == 1)
                            {
                                for (int d=0; d<AMREX_SPACEDIM; ++d)
                                {
                                    part.pos(d) += part.pos(d) + posAlt[d];
                                }
                            }
                        }


                    }else
                    {
                        runtime = runtime - scatterTime;
                        //scatterCountPtr[0]++;
                        amrex::Gpu::Atomic::Add(scatterCountPtr, 1);

                        for (int d=0; d<AMREX_SPACEDIM; ++d)
                        {
                            part.pos(d) += scatterTime * part.rdata(FHD_realData::velx + d)*ADJ;
                        }

                        randomSphere(&part.rdata(FHD_realData::velx),&part.rdata(FHD_realData::vely), &part.rdata(FHD_realData::velz), engine);

                    }


                }

                part.idata(FHD_intData::sorted) = -1;

                part.rdata(FHD_realData::timeFrac) = 1.0;

    //            if(step%2==0)
    //            {
    //                part.rdata(FHD_realData::velz) = -part.rdata(FHD_realData::velz);
    //            }

                if(part.idata(FHD_intData::newSpecies) != -1)
                {
                    part.idata(FHD_intData::species) = part.idata(FHD_intData::newSpecies);
                    part.idata(FHD_intData::newSpecies) = -1;
                }


            });


            //Print() << "Pre buffer size: " << paramPlaneList[1].recCountRight << endl;
            for (int i = 0; i < np; i++)
            {
                ParticleType & part = particles[i];
                if(part.idata(FHD_intData::fluxRec) > 0)
                {

                    int surfnum = part.idata(FHD_intData::fluxRec)-1;
                    int surfcount = paramPlaneList[surfnum].recCountRight;

                    paramPlaneList[surfnum].xVelRecRight[surfcount] = part.rdata(FHD_realData::velx);
                    paramPlaneList[surfnum].yVelRecRight[surfcount] = part.rdata(FHD_realData::vely);
                    paramPlaneList[surfnum].zVelRecRight[surfcount] = part.rdata(FHD_realData::velz);

                    paramPlaneList[surfnum].xPosRecRight[surfcount] = part.pos(0);
                    paramPlaneList[surfnum].yPosRecRight[surfcount] = part.pos(1);
                    paramPlaneList[surfnum].zPosRecRight[surfcount] = part.pos(2);

                    //Print() << "Buffering " << paramPlaneList[surfnum].yPosRecRight[surfcount] << " to surf " << surfnum << endl;

                    paramPlaneList[surfnum].freqRecRight[surfcount] = part.rdata(FHD_realData::omega);

                    paramPlaneList[surfnum].timeRecRight[surfcount] = part.rdata(FHD_realData::travelTime);

                    paramPlaneList[surfnum].recCountRight++;

                    if(paramPlaneList[surfnum].recCountRight == WRITE_BUFFER)
                    {
                        //Print() << "Writing on step " << istep << "\n";
                        std::string plotfilename = std::to_string(surfnum+1) + "_" + std::to_string(ParallelDescriptor::MyProc()) + amrex::Concatenate("_particles_right_",step,12);
                        std::ofstream ofs(plotfilename, std::ios::app);

                        for(int j=0; j <WRITE_BUFFER;j++)
                        {
                            ofs << paramPlaneList[surfnum].xPosRecRight[j] << " " << paramPlaneList[surfnum].yPosRecRight[j] << " " << paramPlaneList[surfnum].zPosRecRight[j]
                                << " " << paramPlaneList[surfnum].xVelRecRight[j] << " " << paramPlaneList[surfnum].yVelRecRight[j] << " " << paramPlaneList[surfnum].zVelRecRight[j]
                                << " " << paramPlaneList[surfnum].freqRecRight[j] << " " << paramPlaneList[surfnum].timeRecRight[j] << std::endl;
                         }
                         ofs.close();
                         paramPlaneList[surfnum].recCountRight = 0;

                    }

                }else if(part.idata(FHD_intData::fluxRec) < 0)
                {

                    int surfnum = -part.idata(FHD_intData::fluxRec)-1;
                    int surfcount = paramPlaneList[surfnum].recCountLeft;

                    paramPlaneList[surfnum].xVelRecLeft[surfcount] = part.rdata(FHD_realData::velx);
                    paramPlaneList[surfnum].yVelRecLeft[surfcount] = part.rdata(FHD_realData::vely);
                    paramPlaneList[surfnum].zVelRecLeft[surfcount] = part.rdata(FHD_realData::velz);

                    paramPlaneList[surfnum].xPosRecLeft[surfcount] = part.pos(0);
                    paramPlaneList[surfnum].yPosRecLeft[surfcount] = part.pos(1);
                    paramPlaneList[surfnum].zPosRecLeft[surfcount] = part.pos(2);

                    paramPlaneList[surfnum].freqRecLeft[surfcount] = part.rdata(FHD_realData::omega);

                    paramPlaneList[surfnum].timeRecLeft[surfcount] = part.rdata(FHD_realData::travelTime);

                    paramPlaneList[surfnum].recCountLeft++;

                    if(paramPlaneList[surfnum].recCountLeft == WRITE_BUFFER)
                    {
                        std::string plotfilename = std::to_string(surfnum+1) + "_" + std::to_string(ParallelDescriptor::MyProc()) + amrex::Concatenate("_particles_left_",step,12);
                        std::ofstream ofs(plotfilename, std::ios::app);
                        for(int j=0; j <WRITE_BUFFER;j++)
                        {

                            ofs << paramPlaneList[surfnum].xPosRecLeft[j] << " " << paramPlaneList[surfnum].yPosRecLeft[j] << " " << paramPlaneList[surfnum].zPosRecLeft[j]
                                << " " << paramPlaneList[surfnum].xVelRecLeft[j] << " " << paramPlaneList[surfnum].yVelRecLeft[j] << " " << paramPlaneList[surfnum].zVelRecLeft[j]
                                << " " << paramPlaneList[surfnum].freqRecLeft[j] << " " << paramPlaneList[surfnum].timeRecLeft[j] << std::endl;

                         }
                         ofs.close();
                         paramPlaneList[surfnum].recCountLeft = 0;

                    }
                }
            }
            //Print() << "Post buffer size: " << paramPlaneList[1].recCountRight << endl;


        }

    //    for(int i = 0; i<MAX_SPECIES;i++)
    //    {
    //        int temp = delCount[i];
    //        ParallelDescriptor::ReduceIntSum(temp);
    //        if(temp != 0)
    //        {
    //            Print() << "Deleted " << temp << " of species " << i << "\n";
    //        }
    //    }

        int scatterCountInt = scatterCountPtr[0];
        int countInt = countPtr[0];
        int specCountInt = specCountPtr[0];
        ParallelDescriptor::ReduceIntSum(scatterCountInt);
        ParallelDescriptor::ReduceIntSum(totalParts);
        ParallelDescriptor::ReduceIntSum(countInt);
        ParallelDescriptor::ReduceIntSum(specCountInt);

        if(istep%100==0)
        {
            Print() << "Total particles: " << totalParts << "\n";
            Print() << "Internal scattering events: " << scatterCountInt << "\n";
            if(countInt != 0)
            {
                Print() << "Fraction of boundary interactions which were specular: " << (double)specCountInt/((double)countInt) << "\n";
            }
        }

        Redistribute();
        SortParticlesDB();
        //SortParticles();

        if (ParallelDescriptor::MyProc() == 0){
            std::string plotfilename = "totalparts";
            std::fstream ofs;
            ofs.open(plotfilename, std::ios::app);
            if (!ofs){
                std::ofstream ofs(plotfilename, std::ios::app);
                ofs << totalParts << std::endl;
                ofs.close();
            }
            else{
                ofs << totalParts << std::endl;
                ofs.close();
            }
        }
    }

    void FhdParticleContainer::SortParticles()
    {
        int lev = 0;

        const Real* dx = Geom(lev).CellSize();
        const Real* plo = Geom(lev).ProbLo();
        const Real* phi = Geom(lev).ProbHi();

        for (FhdParIter pti(* this, lev); pti.isValid(); ++pti)
        {
            const int grid_id = pti.index();
            const int tile_id = pti.LocalTileIndex();
            const Box& tile_box  = pti.tilebox();

            auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];
            auto& particles = particle_tile.GetArrayOfStructs();
            const long np = particles.numParticles();

            //cout << "Sorting started\n";

            for (int i = 0; i < np; ++ i)
            {
                ParticleType & part = particles[i];

                IntVect iv ={0,0,0};
                if(part.idata(FHD_intData::sorted) == -1)
                {
    //                if(part.id() == 229592)
    //                {
    //                    //cout << "Sorting1!\n";
    //                }
                    iv[0] = (int)floor((part.pos(0)-plo[0])/dx[0]);
                    iv[1] = (int)floor((part.pos(1)-plo[1])/dx[1]);
                    iv[2] = (int)floor((part.pos(2)-plo[2])/dx[2]);

    //                if(part.id() == 229592)
    //                {
    //                    //cout << "Sorting2!\n";
    //                }

                    part.idata(FHD_intData::i) = iv[0];
                    part.idata(FHD_intData::j) = iv[1];
                    part.idata(FHD_intData::k) = iv[2];

    //                if(part.id() == 229592)
    //                {
    //                   // cout << "Cell: " << part.idata(FHD_intData::i) << ", " << part.idata(FHD_intData::j) << ", " << part.idata(FHD_intData::k) << endl;
    //                   // cout << "Species: " << part.idata(FHD_intData::species) << endl;
    //                }

                    long imap = tile_box.index(iv);
                    part.idata(FHD_intData::sorted) = m_cell_vectors[part.idata(FHD_intData::species)][pti.index()][imap].size();
    //                if(part.id() == 229592)
    //                {
    //                   // cout << "Sorting3!\n";
    //                }

                    m_cell_vectors[part.idata(FHD_intData::species)][pti.index()][imap].push_back(i);
                    //if(part.id() == 229592)
                   // {
                      //  cout << "Sorting4!\n";
                    //}
                }
                else
                {
                    iv[0] = part.idata(FHD_intData::i);
                    iv[1] = part.idata(FHD_intData::j);
                    iv[2] = part.idata(FHD_intData::k);
                    long imap = tile_box.index(iv);
                    m_cell_vectors[part.idata(FHD_intData::species)][pti.index()][imap][part.idata(FHD_intData::sorted)] = i;
                }
            }

            //cout << "Sorting finished\n";
        }
    }


    void FhdParticleContainer::SortParticlesDB()
    {

        int lev = 0;

        const Real* dx = Geom(lev).CellSize();
        const GpuArray<Real, 3> dxInv = Geom(lev).InvCellSizeArray();
        //const Real* plo = Geom(lev).ProbLo();
        //const Real* phi = Geom(lev).ProbHi();
        const GpuArray<Real, 3> plo = Geom(lev).ProbLoArray();
        const GpuArray<Real, 3> phi = Geom(lev).ProbHiArray();

        for (FhdParIter pti(* this, lev); pti.isValid(); ++pti)
        {
            const int grid_id = pti.index();
            const int tile_id = pti.LocalTileIndex();
            const Box& tile_box  = pti.tilebox();
            const Box& tile_box2  = pti.tilebox();

            auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];
            auto& particles = particle_tile.GetArrayOfStructs();
            const long np = particles.numParticles();
            auto pstruct_ptr = particles().dataPtr();

            auto cellLo = tile_box.smallEnd();
            auto cellHi = tile_box.bigEnd();


            int ncells = (cellHi[0]-cellLo[0]+1)*(cellHi[1]-cellLo[1]+1)*(cellHi[2]-cellLo[2]+1);
            int nbins = ncells*nspecies;
    //        int specCount[nspecies];
    //
    //        for (int i = 0; i < np; ++ i)
    //        {
    //            int spec = particles[i].idata(FHD_intData::species);
    //            specCount[spec]++;
    //        }

    //        if(ParallelDescriptor::MyProc()==1)
    //        {
    //            AllPrint() << "lo: " << cellLo << ", hi: " << cellHi << ", ncells: " << ncells << endl;
    //            AllPrint() << "lo index: " << tile_box.index(cellLo) << endl;
    //
    //        }
    //        Print() << "np: " << np << ", nbins: " << nbins << ", ncells: " << ncells << endl;
        m_bins.build(np, pstruct_ptr, nbins, getBin{plo, dxInv, tile_box, ncells});
            auto inds = m_bins.permutationPtr();
            auto offs = m_bins.offsetsPtr();
    //
            if(ParallelDescriptor::MyProc()==1)
            {
    //            AllPrint() << "lo: " << cellLo << ", hi: " << cellHi << ", ncells: " << ncells << endl;
    //
    //            //auto iv = {4
    //            for(int i=0;i<5;i++)
    //            {
    //                AllPrint() << "off: " << offs[i] << endl;
    //
    //            }
    //
    //            //AllPrint() << "test cell: " << getParticleCell(particles[inds[i]], plo, dxInv, tile_box) << " map: " << tile_box.index(getParticleCell(particles[inds[i]], plo, dxInv, tile_box)) << endl;
    //            for(int i=0;i<5;i++)
    //            {
    //                AllPrint() << "part: " << inds[i] << ", order: " << particles[inds[i]].idata(FHD_intData::sorted) << ", spec: " << particles[inds[i]].idata(FHD_intData::species)
    //                        << ", cell: " << particles[inds[i]].idata(FHD_intData::i) << ", " <<particles[inds[i]].idata(FHD_intData::j) << ", " << particles[inds[i]].idata(FHD_intData::k) << endl;
    //                AllPrint() << "pos: " << particles[inds[i]].pos(0) << ", " << particles[inds[i]].pos(1) << ", " << particles[inds[i]].pos(2) << ", " << endl;
    //                AllPrint() << "cell: " << getPartCell(particles[inds[i]], plo, dxInv, tile_box) << " imap: " << tile_box.index(getPartCell(particles[inds[i]], plo, dxInv, tile_box)) << endl;
    //                AllPrint() << "bin: " << getBin{plo, dxInv, tile_box, ncells}(particles[inds[i]]) << endl;
    //            }
    //
    //            unsigned int* cellList = getCellList(inds,offs,cellLo,1,tile_box);
    //            unsigned int listSize = getBinSize(offs,cellLo,1,tile_box);
    //            AllPrint() << "size: " << listSize << " list: ";
    //            for(int i=0; i<listSize;i++)
    //            {
    //                AllPrint() << cellList[i] << " ";
    //            }
    //            AllPrint() << endl;
    //            long np_spec = m_cell_vectors[1][grid_id][tile_box.index(cellLo)].size();
    //            AllPrint() << "old size: " << np_spec << endl;
            }


        }


    }

    //void FhdParticleContainer::SpecChange(FhdParticleContainer::ParticleType& part) {
    //    int lev = 0;
    //    bool proc_enter = true;
    //
    //    const Real* dx = Geom(lev).CellSize();
    //    Real smallNumber = dx[0];
    //    if(dx[1] < smallNumber){smallNumber = dx[1];}
    //    if(dx[2] < smallNumber){smallNumber = dx[2];}
    //    smallNumber = smallNumber*0.000001;
    //
    //    for (MFIter mfi = MakeMFIter(lev, true); mfi.isValid(); ++mfi)
    //    {
    //        if(proc_enter && ParallelDescriptor::MyProc()==0)
    //        {
    //            proc_enter = false;//Make sure this runs only once incase of tiling
    //
    //            const int grid_id = mfi.index();
    //            const int tile_id = mfi.LocalTileIndex();
    //            auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];
    //
    //            for(int i = 0; i< paramPlaneCount; i++)
    //            {
    //
    //                if(paramPlaneList[i].sourceRight == 1)
    //                {
    //                    for(int j=nspecies-1; j>=0; j--)
    //                    {
    //                        Real density = paramPlaneList[i].densityRight[j];
    //                        Real temp = paramPlaneList[i].temperatureRight;
    //                        //Real area = paramPlaneList[i].area/ParallelDescriptor::NProcs();
    //                        Real area = paramPlaneList[i].area;
    //                        Real fluxMean = density*area*sqrt(properties[j].R*temp/(2.0*M_PI))/particle_neff;
    //
    //                        //cout << "R: " << properties[j].R << endl;
    //                        //cout << "fluxmean right species " << j << " surface " << i << ": " << fluxMean << "\n";
    //
    //                        //Real elapsedTime = -log(amrex::Random())/fluxMean;
    //                        //int totalFluxInt = 0;
    //                        //while(elapsedTime < dt)
    //                        //{
    //                        //    totalFluxInt++;
    //                        //    elapsedTime += -log(amrex::Random())/fluxMean;
    //                        //}
    //
    //                        int totalFluxInt = amrex::RandomPoisson(dt*fluxMean);
    //
    //                        //Real fluxMean = density*area*sqrt(properties[j].R*temp/(2.0*M_PI))/particle_neff;
    //                        //Real fluxVar = density*area*sqrt(properties[j].R*temp/(2.0*M_PI))/particle_neff;
    //
    //                        //Real totalFlux = dt*fluxMean + sqrt(dt*fluxVar)*amrex::RandomNormal(0.,1.);
    //                        //totalFlux = std::max(totalFlux,0.);
    //
    //                        //int totalFluxInt =  (int)floor(totalFlux);
    //                        //Real totalFluxLeftOver = totalFlux - totalFluxInt;
    //
    //                        //if(amrex::Random() < totalFluxLeftOver)
    //                        //{
    //                        //    totalFluxInt++;
    //                        //}
    //                        //Print() << "Surface " << i << " generating " << totalFluxInt << " of species " << j << " on the right.\n";
    //
    //                        for(int k=0;k<totalFluxInt;k++)
    //                        {
    //                            Real uCoord = amrex::Random()*paramPlaneList[i].uTop;
    //                            Real vCoord = amrex::Random()*paramPlaneList[i].vTop;
    //
    //                            ParticleType p;
    //                            p.id() = ParticleType::NextID();
    //
    //                            p.cpu() = ParallelDescriptor::MyProc();
    //                            p.idata(FHD_intData::sorted) = -1;
    //
    //                            p.idata(FHD_intData::species) = j;
    //
    //                            p.pos(0) = paramPlaneList[i].x0 + paramPlaneList[i].ux*uCoord + paramPlaneList[i].vx*vCoord;
    //                            p.pos(1) = paramPlaneList[i].y0 + paramPlaneList[i].uy*uCoord + paramPlaneList[i].vy*vCoord;
    //                            p.pos(2) = paramPlaneList[i].z0 + paramPlaneList[i].uz*uCoord + paramPlaneList[i].vz*vCoord;
    //
    //                            p.pos(0) = p.pos(0) + smallNumber*paramPlaneList[i].rnx;
    //                            p.pos(1) = p.pos(1) + smallNumber*paramPlaneList[i].rny;
    //                            p.pos(2) = p.pos(2) + smallNumber*paramPlaneList[i].rnz;
    //
    //                            p.rdata(FHD_realData::boostx) = 0;
    //                            p.rdata(FHD_realData::boosty) = 0;
    //                            p.rdata(FHD_realData::boostz) = 0;
    //
    //                            p.idata(FHD_intData::i) = -100;
    //                            p.idata(FHD_intData::j) = -100;
    //                            p.idata(FHD_intData::k) = -100;
    //
    //                            p.rdata(FHD_realData::R) = properties[j].R;
    //                            p.rdata(FHD_realData::timeFrac) = amrex::Random();
    //
    //                            Real srt = sqrt(p.rdata(FHD_realData::R)*temp);
    //                            p.rdata(FHD_realData::velx) = srt*amrex::RandomNormal(0.,1.);
    //                            p.rdata(FHD_realData::vely) = srt*amrex::RandomNormal(0.,1.);
    //                            p.rdata(FHD_realData::velz) = sqrt(2)*srt*sqrt(-log(amrex::Random()));
    //
    //                            const paramPlane surf = paramPlaneList[i];
    //
    //                            rotation(surf.cosThetaRight, surf.sinThetaRight, surf.cosPhiRight, surf.sinPhiRight,
    //                                &p.rdata(FHD_realData::velx), &p.rdata(FHD_realData::vely), &p.rdata(FHD_realData::velz));
    //
    //    //                        cout << p.rdata(FHD_realData::velx) << "," << p.rdata(FHD_realData::vely) << ", " << p.rdata(FHD_realData::velz) << "\n";
    //    //                        cout << p.pos(0) << "," << p.pos(1) << ", " << p.pos(2) << "\n";
    //
    //                            particle_tile.push_back(p);
    //                        }
    //    //                    ParallelDescriptor::ReduceIntSum(totalFluxInt);
    //    //                    if(totalFluxInt != 0)
    //    //                    {
    //    //                        Print() << "Surface " << i << " generated " << totalFluxInt << " of species " << j << "\n";
    //    //                    }
    //                    }
    //                }else if(paramPlaneList[i].sourceLeft == 1)
    //                {
    //                    for(int j=nspecies-1; j>=0; j--)
    //                    {
    //                        Real density = paramPlaneList[i].densityLeft[j];
    //                        Real temp = paramPlaneList[i].temperatureLeft;
    //                        //Real area = paramPlaneList[i].area/ParallelDescriptor::NProcs();
    //                        Real area = paramPlaneList[i].area;
    //                        Real fluxMean = density*area*sqrt(properties[j].R*temp/(2.0*M_PI))/particle_neff;
    //
    //                        //cout << "R: " << properties[j].R << endl;
    //                        //cout << "fluxmean left species " << j << " surface " << i << ": " << fluxMean << "\n";
    //
    //                        //Real elapsedTime = -log(amrex::Random())/fluxMean;
    //                        //int totalFluxInt = 0;
    //                        //while(elapsedTime < dt)
    //                        //{
    //                        //    totalFluxInt++;
    //                        //    elapsedTime += -log(amrex::Random())/fluxMean;
    //                        //}
    //
    //                        int totalFluxInt = amrex::RandomPoisson(dt*fluxMean);
    //
    //                        //Real fluxMean = density*area*sqrt(properties[j].R*temp/(2.0*M_PI))/particle_neff;
    //                        //Real fluxVar = density*area*sqrt(properties[j].R*temp/(2.0*M_PI))/particle_neff;
    //
    //                        //Real totalFlux = dt*fluxMean + sqrt(dt*fluxVar)*amrex::RandomNormal(0.,1.);
    //                        //totalFlux = std::max(totalFlux,0.);
    //
    //                        //int totalFluxInt =  (int)floor(totalFlux);
    //                        //Real totalFluxLeftOver = totalFlux - totalFluxInt;
    //
    //                        //if(amrex::Random() < totalFluxLeftOver)
    //                        //{
    //                        //    totalFluxInt++;
    //                        //}
    //
    //                        for(int k=0;k<totalFluxInt;k++)
    //                        {
    //                            Real uCoord = amrex::Random()*paramPlaneList[i].uTop;
    //                            Real vCoord = amrex::Random()*paramPlaneList[i].vTop;
    //
    //                            ParticleType p;
    //                            p.id() = ParticleType::NextID();
    //
    //                            p.cpu() = ParallelDescriptor::MyProc();
    //                            p.idata(FHD_intData::sorted) = -1;
    //
    //                            p.idata(FHD_intData::species) = j;
    //
    //                            p.pos(0) = paramPlaneList[i].x0 + paramPlaneList[i].ux*uCoord + paramPlaneList[i].vx*vCoord;
    //                            p.pos(1) = paramPlaneList[i].y0 + paramPlaneList[i].uy*uCoord + paramPlaneList[i].vy*vCoord;
    //                            p.pos(2) = paramPlaneList[i].z0 + paramPlaneList[i].uz*uCoord + paramPlaneList[i].vz*vCoord;
    //
    //                            //move the particle slightly off the surface so it doesn't intersect it when it moves
    //                            p.pos(0) = p.pos(0) + smallNumber*paramPlaneList[i].lnx;
    //                            p.pos(1) = p.pos(1) + smallNumber*paramPlaneList[i].lny;
    //                            p.pos(2) = p.pos(2) + smallNumber*paramPlaneList[i].lnz;
    //
    //                            p.rdata(FHD_realData::boostx) = 0;
    //                            p.rdata(FHD_realData::boosty) = 0;
    //                            p.rdata(FHD_realData::boostz) = 0;
    //
    //                            p.idata(FHD_intData::i) = -100;
    //                            p.idata(FHD_intData::j) = -100;
    //                            p.idata(FHD_intData::k) = -100;
    //
    //                            p.rdata(FHD_realData::R) = properties[j].R;
    //                            p.rdata(FHD_realData::timeFrac) = amrex::Random();
    //
    //                            Real srt = sqrt(p.rdata(FHD_realData::R)*temp);
    //                            p.rdata(FHD_realData::velx) = srt*amrex::RandomNormal(0.,1.);
    //                            p.rdata(FHD_realData::vely) = srt*amrex::RandomNormal(0.,1.);
    //                            p.rdata(FHD_realData::velz) = sqrt(2)*srt*sqrt(-log(amrex::Random()));
    //
    //                            const paramPlane surf = paramPlaneList[i];
    //
    //
    //                            rotation(surf.cosThetaLeft, surf.sinThetaLeft, surf.cosPhiLeft, surf.sinPhiLeft,
    //                                                &p.rdata(FHD_realData::velx),&p.rdata(FHD_realData::vely), &p.rdata(FHD_realData::velz));
    //
    //
    //                            particle_tile.push_back(p);
    //                        }
    //    //                    ParallelDescriptor::ReduceIntSum(totalFluxInt);
    //    //                    if(totalFluxInt != 0)
    //    //                    {
    //    //                        Print() << "Surface " << i << " generated " << totalFluxInt << " of species " << j << "\n";
    //    //                    }
    //
    //                    }
    //                }
    //            }
    //        }
    //    }
    //    Redistribute();
    //    SortParticles();
    //}


    void FhdParticleContainer::Source(const Real dt, paramPlane* paramPlaneList, const int paramPlaneCount, MultiFab& mfcuInst) {
        int lev = 0;
        bool proc_enter = true;

        const Real* dx = Geom(lev).CellSize();
        Real smallNumber = dx[0];
        if(dx[1] < smallNumber){smallNumber = dx[1];}
        if(dx[2] < smallNumber){smallNumber = dx[2];}
        smallNumber = smallNumber*0.000001;

        for (MFIter mfi = MakeMFIter(lev, true); mfi.isValid(); ++mfi)
        {
            if(proc_enter)
            {
                proc_enter = false;//Make sure this runs only once incase of tiling

                const int grid_id = mfi.index();
                const int tile_id = mfi.LocalTileIndex();
                auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];

                Array4<Real> cuInst  = mfcuInst[mfi].array();

                for(int i = 0; i< paramPlaneCount; i++)
                {

                    if(paramPlaneList[i].sourceRight == 1)
                    {
                        for(int j=nspecies-1; j>=0; j--)
                        {
                            Real density = paramPlaneList[i].densityRight[j]*rho_lo[0]/properties[j].mass;
                            //Real density = paramPlaneList[i].densityRight[j];
                            Real temp = paramPlaneList[i].temperatureRight;
                            Real area = paramPlaneList[i].area/ParallelDescriptor::NProcs();
                            Real fluxMean = density*area*sqrt(properties[j].R*temp/(2.0*M_PI))/particle_neff;

                            int totalFluxInt = amrex::RandomPoisson(dt*fluxMean);

                            for(int k=0;k<totalFluxInt;k++)
                            {
                                Real uCoord = amrex::Random()*paramPlaneList[i].uTop;
                                Real vCoord = amrex::Random()*paramPlaneList[i].vTop;

                                ParticleType p;
                                p.id() = ParticleType::NextID();

                                p.cpu() = ParallelDescriptor::MyProc();
                                p.idata(FHD_intData::sorted) = -1;

                                p.idata(FHD_intData::species) = j;
                                p.idata(FHD_intData::newSpecies) = -1;

                                p.pos(0) = paramPlaneList[i].x0 + paramPlaneList[i].ux*uCoord + paramPlaneList[i].vx*vCoord;
                                p.pos(1) = paramPlaneList[i].y0 + paramPlaneList[i].uy*uCoord + paramPlaneList[i].vy*vCoord;
                                p.pos(2) = paramPlaneList[i].z0 + paramPlaneList[i].uz*uCoord + paramPlaneList[i].vz*vCoord;

                                p.pos(0) = p.pos(0) + smallNumber*paramPlaneList[i].rnx;
                                p.pos(1) = p.pos(1) + smallNumber*paramPlaneList[i].rny;
                                p.pos(2) = p.pos(2) + smallNumber*paramPlaneList[i].rnz;

                                p.rdata(FHD_realData::boostx) = 0;
                                p.rdata(FHD_realData::boosty) = 0;
                                p.rdata(FHD_realData::boostz) = 0;

                                p.rdata(FHD_realData::mass) = mass[j];

                                p.idata(FHD_intData::i) = -100;
                                p.idata(FHD_intData::j) = -100;
                                p.idata(FHD_intData::k) = -100;

                                p.rdata(FHD_realData::R) = properties[j].R;
                                p.rdata(FHD_realData::timeFrac) = amrex::Random();

                                Real srt = sqrt(p.rdata(FHD_realData::R)*temp);
                                p.rdata(FHD_realData::velx) = srt*amrex::RandomNormal(0.,1.);
                                p.rdata(FHD_realData::vely) = srt*amrex::RandomNormal(0.,1.);
                                p.rdata(FHD_realData::velz) = sqrt(2)*srt*sqrt(-log(amrex::Random()));

                                const paramPlane surf = paramPlaneList[i];

                                rotation(surf.cosThetaRight, surf.sinThetaRight, surf.cosPhiRight, surf.sinPhiRight,
                                    &p.rdata(FHD_realData::velx), &p.rdata(FHD_realData::vely), &p.rdata(FHD_realData::velz));

    //                            cout << p.rdata(FHD_realData::velx) << "," << p.rdata(FHD_realData::vely) << ", " << p.rdata(FHD_realData::velz) << "\n";
    //                            cout << p.pos(0) << "," << p.pos(1) << ", " << p.pos(2) << "\n";

                                particle_tile.push_back(p);
                            }

                        }
                    }else if(paramPlaneList[i].sourceRight == 3)
                    {
                        Real srtV = 0;
                        for(int j=nspecies-1; j>=0; j--)
                        {
                            Real temp = paramPlaneList[i].temperatureRight;
                            Real srtVTemp = sqrt(properties[j].R*temp);
                            if(srtVTemp > srtV){srtV = srtVTemp;}
                        }
                        Real depth = 6.0*srtV*dt;
                        Real area = paramPlaneList[i].area/ParallelDescriptor::NProcs();
                        Real vol = depth*area;

                        for(int j=nspecies-1; j>=0; j--)
                        {
                            Real density = paramPlaneList[i].densityRight[j]*rho_lo[0]/properties[j].mass;
                                            //cout << "n: " << density << endl;
                            //Real density = paramPlaneList[i].densityRight[j];

                            Real xMom = paramPlaneList[i].xMomFluxRight[j];
                            Real yMom = paramPlaneList[i].yMomFluxRight[j];
                            Real zMom = paramPlaneList[i].zMomFluxRight[j];

                            int totalAttemptsInt = amrex::RandomPoisson(density*vol);

                            Real totalMass = totalAttemptsInt*properties[j].mass;

                            Real vel = 0.5*(xMom*paramPlaneList[i].rnx + yMom*paramPlaneList[i].rny + zMom*paramPlaneList[i].rnz)/totalMass;

                            //Print() << "Attempting " << density << ", " << sqrt(area) << ", " << depth << endl;

                            int generated = 0;
                            for(int k=0;k<totalAttemptsInt;k++)
                            {
                                ParticleType p;
                                Real srt = sqrt(properties[j].R*paramPlaneList[i].temperatureRight);
                                p.rdata(FHD_realData::velz) = srt*amrex::RandomNormal(0.,1.);
                                Real wCoord = -amrex::Random()*depth;
                                Real projPos = wCoord + dt*p.rdata(FHD_realData::velz);

                                if(projPos > 0)
                                {
                                    Real uCoord = amrex::Random()*paramPlaneList[i].uTop;
                                    Real vCoord = amrex::Random()*paramPlaneList[i].vTop;

                                    p.id() = ParticleType::NextID();

                                    p.cpu() = ParallelDescriptor::MyProc();
                                    p.idata(FHD_intData::sorted) = -1;

                                    p.idata(FHD_intData::species) = j;
                                    p.idata(FHD_intData::newSpecies) = -1;

                                    p.pos(0) = paramPlaneList[i].x0 + paramPlaneList[i].ux*uCoord + paramPlaneList[i].vx*vCoord;
                                    p.pos(1) = paramPlaneList[i].y0 + paramPlaneList[i].uy*uCoord + paramPlaneList[i].vy*vCoord;
                                    p.pos(2) = paramPlaneList[i].z0 + paramPlaneList[i].uz*uCoord + paramPlaneList[i].vz*vCoord;

                                    p.pos(0) = p.pos(0) + smallNumber*paramPlaneList[i].rnx;
                                    p.pos(1) = p.pos(1) + smallNumber*paramPlaneList[i].rny;
                                    p.pos(2) = p.pos(2) + smallNumber*paramPlaneList[i].rnz;

                                    p.rdata(FHD_realData::boostx) = 0;
                                    p.rdata(FHD_realData::boosty) = 0;
                                    p.rdata(FHD_realData::boostz) = 0;

                                    p.rdata(FHD_realData::mass) = mass[j];

                                    p.idata(FHD_intData::i) = -100;
                                    p.idata(FHD_intData::j) = -100;
                                    p.idata(FHD_intData::k) = -100;

                                    p.rdata(FHD_realData::R) = properties[j].R;
                                    p.rdata(FHD_realData::timeFrac) = (dt + wCoord/p.rdata(FHD_realData::velz))/dt;

                                    p.rdata(FHD_realData::velx) = srt*amrex::RandomNormal(0.,1.);
                                    p.rdata(FHD_realData::vely) = srt*amrex::RandomNormal(0.,1.);


                                    const paramPlane surf = paramPlaneList[i];

                                    rotation(surf.cosThetaRight, surf.sinThetaRight, surf.cosPhiRight, surf.sinPhiRight,
                                        &p.rdata(FHD_realData::velx), &p.rdata(FHD_realData::vely), &p.rdata(FHD_realData::velz));

                                    particle_tile.push_back(p);
                                    generated++;
                                }
                            }
                            //Print() << "Generated " << generated << endl;

                        }
                    }else if(paramPlaneList[i].sourceRight == 2)
                    {
                        Real totalMass = 0;
                        int rankTotalParticles = 0;
                        Real residual = 0;
                        int particleFluxInt[nspecies];
                        //Real particleFluxReal[nspecies];
                        Real massFluxReal[nspecies];

                        for(int j=nspecies-1; j>=0; j--)
                        {
                            totalMass += paramPlaneList[i].massFluxRight[j];
                        }

                        //Print() << "Total mass out: " << totalMass << endl;

                        for(int j=nspecies-1; j>=0; j--)
                        {
                            massFluxReal[j] = bc_Yk_x_lo[j]*totalMass/((double)ParallelDescriptor::NProcs());
                            //particleFluxReal[j] = massFluxReal[j]/mass[j];
                            particleFluxInt[j] = (int)floor(massFluxReal[j]/mass[j]);
                            rankTotalParticles += particleFluxInt[j];
                            residual += massFluxReal[j] - particleFluxInt[j]*mass[j];
                        }
                        for(int j=nspecies-1; j>=0; j--)
                        {
                            paramPlaneList[i].massFluxRight[j] = bc_Yk_x_lo[j]*residual;
                            //paramPlaneList[i].massFluxRight[j] = 0;
                        }

                        ParticleType partList[rankTotalParticles];
                        int partCount = 0;
                        Real generatedMass = 0;
                        for(int j=nspecies-1; j>=0; j--)
                        {

                            Real temp = paramPlaneList[i].temperatureRight;

                            for(int k=0;k<particleFluxInt[j];k++)
                            {

                                Real uCoord = amrex::Random()*paramPlaneList[i].uTop;
                                Real vCoord = amrex::Random()*paramPlaneList[i].vTop;

                                partList[partCount].id() = ParticleType::NextID();

                                partList[partCount].cpu() = ParallelDescriptor::MyProc();
                                partList[partCount].idata(FHD_intData::sorted) = -1;

                                partList[partCount].idata(FHD_intData::species) = j;
                                partList[partCount].idata(FHD_intData::newSpecies) = -1;

                                partList[partCount].pos(0) = paramPlaneList[i].x0 + paramPlaneList[i].ux*uCoord + paramPlaneList[i].vx*vCoord;
                                partList[partCount].pos(1) = paramPlaneList[i].y0 + paramPlaneList[i].uy*uCoord + paramPlaneList[i].vy*vCoord;
                                partList[partCount].pos(2) = paramPlaneList[i].z0 + paramPlaneList[i].uz*uCoord + paramPlaneList[i].vz*vCoord;

                                partList[partCount].pos(0) = partList[partCount].pos(0) + smallNumber*paramPlaneList[i].rnx;
                                partList[partCount].pos(1) = partList[partCount].pos(1) + smallNumber*paramPlaneList[i].rny;
                                partList[partCount].pos(2) = partList[partCount].pos(2) + smallNumber*paramPlaneList[i].rnz;

                                partList[partCount].rdata(FHD_realData::boostx) = 0;
                                partList[partCount].rdata(FHD_realData::boosty) = 0;
                                partList[partCount].rdata(FHD_realData::boostz) = 0;

                                partList[partCount].rdata(FHD_realData::mass) = mass[j];

                                generatedMass += mass[j];

                                partList[partCount].idata(FHD_intData::i) = -100;
                                partList[partCount].idata(FHD_intData::j) = -100;
                                partList[partCount].idata(FHD_intData::k) = -100;

                                partList[partCount].rdata(FHD_realData::R) = properties[j].R;
                                partList[partCount].rdata(FHD_realData::timeFrac) = amrex::Random();

                                Real srt = sqrt(partList[partCount].rdata(FHD_realData::R)*temp);
                                partList[partCount].rdata(FHD_realData::velx) = srt*amrex::RandomNormal(0.,1.);
                                partList[partCount].rdata(FHD_realData::vely) = srt*amrex::RandomNormal(0.,1.);
                                partList[partCount].rdata(FHD_realData::velz) = sqrt(2)*srt*sqrt(-log(amrex::Random()));

                                const paramPlane surf = paramPlaneList[i];

                                rotation(surf.cosThetaRight, surf.sinThetaRight, surf.cosPhiRight, surf.sinPhiRight,
                                    &partList[partCount].rdata(FHD_realData::velx), &partList[partCount].rdata(FHD_realData::vely), &partList[partCount].rdata(FHD_realData::velz));

                                partCount++;

                            }

                            ParallelDescriptor::ReduceRealSum(paramPlaneList[i].massFluxRight[j]);
                        }


                        ParallelDescriptor::ReduceRealSum(generatedMass);
                        //Print() << "Total mass in: " << generatedMass << endl;


                        for(int l=0; l<rankTotalParticles; l++)
                        {
                            particle_tile.push_back(partList[l]);
                        }

                    }

                    if(paramPlaneList[i].sourceLeft == 1)
                    {
                        //Print() << "Type 1\n";
                        for(int j=nspecies-1; j>=0; j--)
                        {
                            Real density = paramPlaneList[i].densityLeft[j]*rho_lo[0]/properties[j].mass;
                            //Real density = paramPlaneList[i].densityLeft[j];
                            Real temp = paramPlaneList[i].temperatureLeft;
                            Real area = paramPlaneList[i].area/ParallelDescriptor::NProcs();
                            Real fluxMean = density*area*sqrt(properties[j].R*temp/(2.0*M_PI))/particle_neff;

                            int totalFluxInt = amrex::RandomPoisson(dt*fluxMean);

                            for(int k=0;k<totalFluxInt;k++)
                            {
                                Real uCoord = amrex::Random()*paramPlaneList[i].uTop;
                                Real vCoord = amrex::Random()*paramPlaneList[i].vTop;

                                ParticleType p;
                                p.id() = ParticleType::NextID();

                                p.cpu() = ParallelDescriptor::MyProc();
                                p.idata(FHD_intData::sorted) = -1;

                                p.idata(FHD_intData::species) = j;
                                p.idata(FHD_intData::newSpecies) = -1;

                                p.pos(0) = paramPlaneList[i].x0 + paramPlaneList[i].ux*uCoord + paramPlaneList[i].vx*vCoord;
                                p.pos(1) = paramPlaneList[i].y0 + paramPlaneList[i].uy*uCoord + paramPlaneList[i].vy*vCoord;
                                p.pos(2) = paramPlaneList[i].z0 + paramPlaneList[i].uz*uCoord + paramPlaneList[i].vz*vCoord;

                                //move the particle slightly off the surface so it doesn't intersect it when it moves
                                p.pos(0) = p.pos(0) + smallNumber*paramPlaneList[i].lnx;
                                p.pos(1) = p.pos(1) + smallNumber*paramPlaneList[i].lny;
                                p.pos(2) = p.pos(2) + smallNumber*paramPlaneList[i].lnz;

                                p.rdata(FHD_realData::boostx) = 0;
                                p.rdata(FHD_realData::boosty) = 0;
                                p.rdata(FHD_realData::boostz) = 0;

                                p.rdata(FHD_realData::mass) = mass[j];

                                p.idata(FHD_intData::i) = -100;
                                p.idata(FHD_intData::j) = -100;
                                p.idata(FHD_intData::k) = -100;

                                p.rdata(FHD_realData::R) = properties[j].R;
                                p.rdata(FHD_realData::timeFrac) = amrex::Random();

                                Real srt = sqrt(p.rdata(FHD_realData::R)*temp);
                                p.rdata(FHD_realData::velx) = srt*amrex::RandomNormal(0.,1.);
                                p.rdata(FHD_realData::vely) = srt*amrex::RandomNormal(0.,1.);
                                p.rdata(FHD_realData::velz) = sqrt(2)*srt*sqrt(-log(amrex::Random()));

                                const paramPlane surf = paramPlaneList[i];


                                rotation(surf.cosThetaLeft, surf.sinThetaLeft, surf.cosPhiLeft, surf.sinPhiLeft,
                                                    &p.rdata(FHD_realData::velx),&p.rdata(FHD_realData::vely), &p.rdata(FHD_realData::velz));


                                particle_tile.push_back(p);
                            }


                        }
                    }else if(paramPlaneList[i].sourceLeft == 3)
                    {
                        Real srtV = 0;
                        for(int j=nspecies-1; j>=0; j--)
                        {
                            Real temp = paramPlaneList[i].temperatureLeft;
                            Real srtVTemp = sqrt(properties[j].R*temp);
                            if(srtVTemp > srtV){srtV = srtVTemp;}
                        }
                        Real depth = 6.0*srtV*dt;
                        Real area = paramPlaneList[i].area/ParallelDescriptor::NProcs();
                        Real vol = depth*area;

                        for(int j=nspecies-1; j>=0; j--)
                        {
                            Real density = paramPlaneList[i].densityLeft[j]*rho_lo[0]/properties[j].mass;
                            //Real density = paramPlaneList[i].densityLeft[j];

                            int totalAttemptsInt = amrex::RandomPoisson(density*vol);

                            for(int k=0;k<totalAttemptsInt;k++)
                            {
                                ParticleType p;
                                Real srt = sqrt(properties[j].R*paramPlaneList[i].temperatureRight);
                                p.rdata(FHD_realData::velz) = srt*amrex::RandomNormal(0.,1.);
                                Real wCoord = -amrex::Random()*depth;
                                Real projPos = wCoord + dt*p.rdata(FHD_realData::velz);

                                if(projPos > 0)
                                {
                                    Real uCoord = amrex::Random()*paramPlaneList[i].uTop;
                                    Real vCoord = amrex::Random()*paramPlaneList[i].vTop;

                                    p.id() = ParticleType::NextID();

                                    p.cpu() = ParallelDescriptor::MyProc();
                                    p.idata(FHD_intData::sorted) = -1;

                                    p.idata(FHD_intData::species) = j;
                                    p.idata(FHD_intData::newSpecies) = -1;

                                    p.pos(0) = paramPlaneList[i].x0 + paramPlaneList[i].ux*uCoord + paramPlaneList[i].vx*vCoord;
                                    p.pos(1) = paramPlaneList[i].y0 + paramPlaneList[i].uy*uCoord + paramPlaneList[i].vy*vCoord;
                                    p.pos(2) = paramPlaneList[i].z0 + paramPlaneList[i].uz*uCoord + paramPlaneList[i].vz*vCoord;

                                    p.pos(0) = p.pos(0) + smallNumber*paramPlaneList[i].lnx;
                                    p.pos(1) = p.pos(1) + smallNumber*paramPlaneList[i].lny;
                                    p.pos(2) = p.pos(2) + smallNumber*paramPlaneList[i].lnz;

                                    p.rdata(FHD_realData::boostx) = 0;
                                    p.rdata(FHD_realData::boosty) = 0;
                                    p.rdata(FHD_realData::boostz) = 0;

                                    p.rdata(FHD_realData::mass) = mass[j];

                                    p.idata(FHD_intData::i) = -100;
                                    p.idata(FHD_intData::j) = -100;
                                    p.idata(FHD_intData::k) = -100;

                                    p.rdata(FHD_realData::R) = properties[j].R;
                                    p.rdata(FHD_realData::timeFrac) = (dt + wCoord/p.rdata(FHD_realData::velz))/dt;

                                    p.rdata(FHD_realData::velx) = srt*amrex::RandomNormal(0.,1.);
                                    p.rdata(FHD_realData::vely) = srt*amrex::RandomNormal(0.,1.);


                                    const paramPlane surf = paramPlaneList[i];

                                    rotation(surf.cosThetaLeft, surf.sinThetaLeft, surf.cosPhiLeft, surf.sinPhiLeft,
                                        &p.rdata(FHD_realData::velx), &p.rdata(FHD_realData::vely), &p.rdata(FHD_realData::velz));

                                    particle_tile.push_back(p);
                                }
                            }

                        }
                    }else if(paramPlaneList[i].sourceLeft == 2)
                    {
                        Real totalMass = 0;
                        int rankTotalParticles = 0;
                        Real residual = 0;
                        int particleFluxInt[nspecies];
                        //Real particleFluxReal[nspecies];
                        Real massFluxReal[nspecies];

                        for(int j=nspecies-1; j>=0; j--)
                        {
                            totalMass += paramPlaneList[i].massFluxLeft[j];
                        }

                        for(int j=nspecies-1; j>=0; j--)
                        {
                            massFluxReal[j] = bc_Yk_x_hi[j]*totalMass/((double)ParallelDescriptor::NProcs());
                            //particleFluxReal[j] = massFluxReal[j]/mass[j];
                            particleFluxInt[j] = (int)floor(massFluxReal[j]/mass[j]);
                            rankTotalParticles += particleFluxInt[j];
                            residual += massFluxReal[j] - particleFluxInt[j]*mass[j];
                        }
                        for(int j=nspecies-1; j>=0; j--)
                        {
                            paramPlaneList[i].massFluxLeft[j] = bc_Yk_x_hi[j]*residual;
                            //paramPlaneList[i].massFluxLeft[j] = 0;
                        }

                        ParticleType partList[rankTotalParticles];
                        int partCount = 0;

                        for(int j=nspecies-1; j>=0; j--)
                        {

                            Real temp = paramPlaneList[i].temperatureLeft;

                            for(int k=0;k<particleFluxInt[j];k++)
                            {
                                Real uCoord = amrex::Random()*paramPlaneList[i].uTop;
                                Real vCoord = amrex::Random()*paramPlaneList[i].vTop;

                                partList[partCount].id() = ParticleType::NextID();

                                partList[partCount].cpu() = ParallelDescriptor::MyProc();
                                partList[partCount].idata(FHD_intData::sorted) = -1;

                                partList[partCount].idata(FHD_intData::species) = j;
                                partList[partCount].idata(FHD_intData::newSpecies) = -1;

                                partList[partCount].pos(0) = paramPlaneList[i].x0 + paramPlaneList[i].ux*uCoord + paramPlaneList[i].vx*vCoord;
                                partList[partCount].pos(1) = paramPlaneList[i].y0 + paramPlaneList[i].uy*uCoord + paramPlaneList[i].vy*vCoord;
                                partList[partCount].pos(2) = paramPlaneList[i].z0 + paramPlaneList[i].uz*uCoord + paramPlaneList[i].vz*vCoord;

                                //move the particle slightly off the surface so it doesn't intersect it when it moves
                                partList[partCount].pos(0) = partList[partCount].pos(0) + smallNumber*paramPlaneList[i].lnx;
                                partList[partCount].pos(1) = partList[partCount].pos(1) + smallNumber*paramPlaneList[i].lny;
                                partList[partCount].pos(2) = partList[partCount].pos(2) + smallNumber*paramPlaneList[i].lnz;

                                partList[partCount].rdata(FHD_realData::boostx) = 0;
                                partList[partCount].rdata(FHD_realData::boosty) = 0;
                                partList[partCount].rdata(FHD_realData::boostz) = 0;

                                partList[partCount].rdata(FHD_realData::mass) = mass[j];

                                partList[partCount].idata(FHD_intData::i) = -100;
                                partList[partCount].idata(FHD_intData::j) = -100;
                                partList[partCount].idata(FHD_intData::k) = -100;

                                partList[partCount].rdata(FHD_realData::R) = properties[j].R;
                                partList[partCount].rdata(FHD_realData::timeFrac) = amrex::Random();

                                Real srt = sqrt(partList[partCount].rdata(FHD_realData::R)*temp);
                                partList[partCount].rdata(FHD_realData::velx) = srt*amrex::RandomNormal(0.,1.);
                                partList[partCount].rdata(FHD_realData::vely) = srt*amrex::RandomNormal(0.,1.);
                                partList[partCount].rdata(FHD_realData::velz) = sqrt(2)*srt*sqrt(-log(amrex::Random()));

                                const paramPlane surf = paramPlaneList[i];


                                rotation(surf.cosThetaLeft, surf.sinThetaLeft, surf.cosPhiLeft, surf.sinPhiLeft,
                                                    &partList[partCount].rdata(FHD_realData::velx),&partList[partCount].rdata(FHD_realData::vely), &partList[partCount].rdata(FHD_realData::velz));

                                partCount++;
                            }

                            ParallelDescriptor::ReduceRealSum(paramPlaneList[i].massFluxLeft[j]);

                        }


                        for(int l=0; l<rankTotalParticles; l++)
                        {
                            particle_tile.push_back(partList[l]);

                            //cout << "Pushing back " << partList[l].rdata(FHD_realData::velx) << ", " << partList[l].rdata(FHD_realData::mass) << endl;
                        }


                    }
                }
            }
        }
        Redistribute();
        //SortParticles();
        SortParticlesDB();
    }

    void FhdParticleContainer::SourcePhonons(const Real dt, const paramPlane* paramPlaneList, const int paramPlaneCount) {
        int lev = 0;
        bool proc_enter = true;

        const Real* dx = Geom(lev).CellSize();
        Real smallNumber = dx[0];
        if(dx[1] < smallNumber){smallNumber = dx[1];}
        if(dx[2] < smallNumber){smallNumber = dx[2];}
        smallNumber = smallNumber*0.00000001;

        int procID = ParallelDescriptor::MyProc();
        amrex::RandomEngine engine;

        for(MFIter mfi = MakeMFIter(lev, true); mfi.isValid(); ++mfi)
        {
            if(proc_enter)
            {
                proc_enter = false;//Make sure this runs only once incase of tiling

                const int grid_id = mfi.index();
                const int tile_id = mfi.LocalTileIndex();
                auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];
                auto& aos = particle_tile.GetArrayOfStructs();
                ParticleType* particles = aos().dataPtr();

                Real pSpeed = phonon_sound_speed;

                Gpu::ManagedVector<paramPlane> paramPlaneListTmp;
                paramPlaneListTmp.resize(paramPlaneCount);
                for(int i=0;i<paramPlaneCount;i++)
                {
                    paramPlaneListTmp[i]=paramPlaneList[i];

                }
                paramPlane* paramPlaneListPtr = paramPlaneListTmp.data();

                for(int i = 0; i< paramPlaneCount; i++)
                {
                    if(paramPlaneList[i].sourceLeft == 1)
                    {

                        for(int j = 0; j< nspecies; j++)
                        {
                            Real density = paramPlaneList[i].densityLeft[j];
                            Real temp = paramPlaneList[i].temperatureLeft;
                            Real area = paramPlaneList[i].area/ParallelDescriptor::NProcs();
                            Real fluxMean = density*area;

                            Real totalFlux = dt*fluxMean;

                            int totalFluxInt =  (int)floor(totalFlux);
                            Real totalFluxLeftOver = totalFlux - totalFluxInt;

                            if(amrex::Random() < totalFluxLeftOver)
                            {
                                totalFluxInt++;
                            }
                            //Print() << "Flux: " << totalFluxInt << endl;
                            //auto old_size = particle_tile.GetArrayOfStructs().size();
                            //auto new_size = old_size + totalFluxInt;
                            //particle_tile.resize(new_size);
                            //ParticleType* pstruct = particle_tile.GetArrayOfStructs()().data();

                            //amrex::ParallelForRNG(totalFluxInt, [=] AMREX_GPU_DEVICE (int k, amrex::RandomEngine const& engine) noexcept
                            for(int k=0;k<totalFluxInt;k++)
                            {
                                Real uCoord = amrex::Random(engine)*paramPlaneListPtr[i].uTop;
                                Real vCoord = amrex::Random(engine)*paramPlaneListPtr[i].vTop;

                                ParticleType p;
    //                            ParticleType& p = pstruct[old_size+k];

                                p.id() = ParticleType::NextID();
                                //p.id() = old_size+k;
    //                            p.cpu() = ParallelDescriptor::MyProc();
                                p.cpu() = procID;
                                p.idata(FHD_intData::sorted) = -1;

                                p.idata(FHD_intData::species) = j;
                                p.idata(FHD_intData::newSpecies) = -1;

                                p.pos(0) = paramPlaneListPtr[i].x0 + paramPlaneListPtr[i].ux*uCoord + paramPlaneListPtr[i].vx*vCoord;
                                p.pos(1) = paramPlaneListPtr[i].y0 + paramPlaneListPtr[i].uy*uCoord + paramPlaneListPtr[i].vy*vCoord;
                                p.pos(2) = paramPlaneListPtr[i].z0 + paramPlaneListPtr[i].uz*uCoord + paramPlaneListPtr[i].vz*vCoord;

                                //move the particle slightly off the surface so it doesn't intersect it when it moves
                                p.pos(0) = p.pos(0) + smallNumber*paramPlaneListPtr[i].lnx;
                                p.pos(1) = p.pos(1) + smallNumber*paramPlaneListPtr[i].lny;
                                p.pos(2) = p.pos(2) + smallNumber*paramPlaneListPtr[i].lnz;


                                p.idata(FHD_intData::i) = -100;
                                p.idata(FHD_intData::j) = -100;
                                p.idata(FHD_intData::k) = -100;

                                p.idata(FHD_intData::fluxRec) = 0;

                                if(toggleTimeFrac = 1)
                                {
                                    p.rdata(FHD_realData::timeFrac) = amrex::Random(engine);            //boundary 'reservoir' source
                                }else
                                {
                                    p.rdata(FHD_realData::timeFrac) = 1;
                                }
                                p.rdata(FHD_realData::travelTime) = 0;

                                const paramPlane surf = paramPlaneListPtr[i];

                                if(toggleTimeFrac = 1)
                                {
                                    cosineLawHemisphere(surf.cosThetaLeft, surf.sinThetaLeft, surf.cosPhiLeft, surf.sinPhiLeft,
                                                    &p.rdata(FHD_realData::velx),&p.rdata(FHD_realData::vely), &p.rdata(FHD_realData::velz), pSpeed, engine);
                                }else
                                {
                                    randomhemisphere(surf.cosThetaLeft, surf.sinThetaLeft, surf.cosPhiLeft, surf.sinPhiLeft,
                                                    &p.rdata(FHD_realData::velx),&p.rdata(FHD_realData::vely), &p.rdata(FHD_realData::velz), engine);
                                }


                                p.rdata(FHD_realData::omega) = plankDist(surf.temperatureLeft, engine);
                                //I hope this is right?
                                p.rdata(FHD_realData::lambda) = pSpeed*2.0*M_PI/p.rdata(FHD_realData::omega);

                                //Print() << "Gen " << p.id() << ": " << p.rdata(FHD_realData::velx + 0) << ", " << p.rdata(FHD_realData::velx + 1) << ", " << p.rdata(FHD_realData::velx + 2) << endl;

                                //cout << "velx: " << p.rdata(FHD_realData::velx) << "\n";


                                particle_tile.push_back(p);

                            }


                        }
                    }
                    else if(paramPlaneList[i].sourceRight == 1)
                    {
                        for(int j=0; j< nspecies; j++)
                        {
                            Real density = paramPlaneList[i].densityRight[j];
                            Real temp = paramPlaneList[i].temperatureRight;
                            Real area = paramPlaneList[i].area/ParallelDescriptor::NProcs();
                            Real fluxMean = density*area;

                            Real totalFlux = dt*fluxMean;

                            int totalFluxInt =  (int)floor(totalFlux);
                            Real totalFluxLeftOver = totalFlux - totalFluxInt;

                            if(amrex::Random() < totalFluxLeftOver)
                            {
                                totalFluxInt++;
                            }

                            //auto old_size = particle_tile.GetArrayOfStructs().size();
                            //auto new_size = old_size + totalFluxInt;
                            //particle_tile.resize(new_size);
                            //ParticleType* pstruct = particle_tile.GetArrayOfStructs()().data();

                            //amrex::ParallelForRNG(totalFluxInt, [=] AMREX_GPU_DEVICE (int k, amrex::RandomEngine const& engine) noexcept
                            for(int k=0;k<totalFluxInt;k++)
                            {
                                Real uCoord = amrex::Random(engine)*paramPlaneListPtr[i].uTop;
                                Real vCoord = amrex::Random(engine)*paramPlaneListPtr[i].vTop;

                                ParticleType p;
                                //ParticleType& p = pstruct[old_size+k];
                                p.id() = ParticleType::NextID();
                                //p.id() = old_size+k;

                                p.cpu() = procID;
                                p.idata(FHD_intData::sorted) = -1;

                                p.idata(FHD_intData::species) = j;
                                p.idata(FHD_intData::newSpecies) = -1;

                                p.pos(0) = paramPlaneList[i].x0 + paramPlaneListPtr[i].ux*uCoord + paramPlaneListPtr[i].vx*vCoord;
                                p.pos(1) = paramPlaneList[i].y0 + paramPlaneListPtr[i].uy*uCoord + paramPlaneListPtr[i].vy*vCoord;
                                p.pos(2) = paramPlaneList[i].z0 + paramPlaneListPtr[i].uz*uCoord + paramPlaneListPtr[i].vz*vCoord;

                                p.pos(0) = p.pos(0) + smallNumber*paramPlaneListPtr[i].rnx;
                                p.pos(1) = p.pos(1) + smallNumber*paramPlaneListPtr[i].rny;
                                p.pos(2) = p.pos(2) + smallNumber*paramPlaneListPtr[i].rnz;

                                p.rdata(FHD_realData::boostx) = 0;
                                p.rdata(FHD_realData::boosty) = 0;
                                p.rdata(FHD_realData::boostz) = 0;

                                p.idata(FHD_intData::i) = -100;
                                p.idata(FHD_intData::j) = -100;
                                p.idata(FHD_intData::k) = -100;

                                p.idata(FHD_intData::fluxRec) = 0;

                                if(toggleTimeFrac = 1)
                                {
                                    p.rdata(FHD_realData::timeFrac) = amrex::Random(engine);            //boundary 'reservoir' source
                                }else
                                {
                                    p.rdata(FHD_realData::timeFrac) = 1;
                                }
                                p.rdata(FHD_realData::travelTime) = 0;


                                const paramPlane surf = paramPlaneListPtr[i];

                                if(toggleTimeFrac = 1)
                                {
                                    cosineLawHemisphere(surf.cosThetaRight, surf.sinThetaRight, surf.cosPhiRight, surf.sinPhiRight,
                                                    &p.rdata(FHD_realData::velx),&p.rdata(FHD_realData::vely), &p.rdata(FHD_realData::velz), pSpeed, engine);
                                }else
                                {
                                    randomhemisphere(surf.cosThetaRight, surf.sinThetaRight, surf.cosPhiRight, surf.sinPhiRight,
                                                    &p.rdata(FHD_realData::velx),&p.rdata(FHD_realData::vely), &p.rdata(FHD_realData::velz), engine);
                                }

                                p.rdata(FHD_realData::omega) = plankDist(surf.temperatureRight, engine);
                                //I hope this is right?
                                p.rdata(FHD_realData::lambda) = pSpeed*2.0*M_PI/p.rdata(FHD_realData::omega);

                                particle_tile.push_back(p);
                            }
//                            ParallelDescriptor::ReduceIntSum(totalFluxInt);
//                            if(totalFluxInt != 0)
//                            {
//                                Print() << "Surface " << i << " generated " << totalFluxInt << " of species " << j << "\n";
//                            }
                        }
                    }
                }
            }
        }
        Redistribute();
        //SortParticles();
        SortParticlesDB();
    }

    void FhdParticleContainer::zeroCells()
    {
        BL_PROFILE_VAR("zeroCells()",zeroCells);
        int lev = 0;
        for(MFIter mfi(mfvrmax); mfi.isValid(); ++mfi)
        {
            const Box& tile_box  = mfi.tilebox();
            const int grid_id = mfi.index();
            const int tile_id = mfi.LocalTileIndex();
            auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];
            auto& particles = particle_tile.GetArrayOfStructs();

            const Array4<Real> & arrvrmax = mfvrmax.array(mfi);
            const Array4<Real> & arrselect = mfselect.array(mfi);

            auto inds = m_bins.permutationPtr();
            auto offs = m_bins.offsetsPtr();
            //const long np = particles.numParticles();
            //amrex::ParallelForRNG(tile_box,
            //    [=] AMREX_GPU_DEVICE (int i, int j, int k, amrex::RandomEngine const& engine) noexcept {
            IntVect smallEnd = tile_box.smallEnd();
            IntVect bigEnd = tile_box.bigEnd();

            for (int i = smallEnd[0]; i <= bigEnd[0]; i++) {
            for (int j = smallEnd[1]; j <= bigEnd[1]; j++) {
            for (int k = smallEnd[2]; k <= bigEnd[2]; k++) {

                const IntVect& iv = {i,j,k};
                long imap = tile_box.index(iv);

                unsigned int np[nspecies];

                Real uComVel=0;
                Real vComVel=0;
                Real wComVel=0;

                Real cellMass = 0;
                int specTotal = 0;

                for(int i_spec = 0; i_spec<nspecies; i_spec++)
                {
                    unsigned int* cellList = getCellList(inds,offs,iv,i_spec,tile_box);
                    np[i_spec] = getBinSize(offs,iv,i_spec,tile_box);
                    for(int i_cell = 0; i_cell<np[i_spec]; i_cell++)
                    {
                        int pindex = cellList[i_cell];
                        ParticleType & p = particles[pindex];

                        uComVel += p.rdata(FHD_realData::mass)*p.rdata(FHD_realData::velx);
                        vComVel += p.rdata(FHD_realData::mass)*p.rdata(FHD_realData::vely);
                        wComVel += p.rdata(FHD_realData::mass)*p.rdata(FHD_realData::velz);

                        cellMass += p.rdata(FHD_realData::mass);
                        specTotal++;

                    }

                }

                uComVel /= cellMass;
                vComVel /= cellMass;
                wComVel /= cellMass;

                Real comTemp = 0;

                for(int i_spec = 0; i_spec<nspecies; i_spec++)
                {
                    unsigned int* cellList = getCellList(inds,offs,iv,i_spec,tile_box);
                    for(int i_cell = 0; i_cell<np[i_spec]; i_cell++)
                    {
                        int pindex = cellList[i_cell];
                        ParticleType & p = particles[pindex];

                        p.rdata(FHD_realData::velx) -= uComVel;
                        p.rdata(FHD_realData::vely) -= vComVel;
                        p.rdata(FHD_realData::velz) -= wComVel;

                    }

                }
                for(int i_spec = 0; i_spec<nspecies; i_spec++)
                {
                    unsigned int* cellList = getCellList(inds,offs,iv,i_spec,tile_box);
                    for(int i_cell = 0; i_cell<np[i_spec]; i_cell++)
                    {
                        int pindex = cellList[i_cell];
                        ParticleType & p = particles[pindex];

                        comTemp += p.rdata(FHD_realData::mass)*(p.rdata(FHD_realData::velx)*p.rdata(FHD_realData::velx) + p.rdata(FHD_realData::vely)*p.rdata(FHD_realData::vely) + p.rdata(FHD_realData::velz)*p.rdata(FHD_realData::velz));
                    }

                }
                comTemp /= (3.0*k_B*specTotal);

                Real tRatio = sqrt(T_init[0]/comTemp);
                Real cellAvMass = cellMass/specTotal;

                for(int i_spec = 0; i_spec<nspecies; i_spec++)
                {
                    unsigned int* cellList = getCellList(inds,offs,iv,i_spec,tile_box);
                    for(int i_cell = 0; i_cell<np[i_spec]; i_cell++)
                    {
                        int pindex = cellList[i_cell];
                        ParticleType & p = particles[pindex];

                        p.rdata(FHD_realData::velx) *= (tRatio*cellAvMass/p.rdata(FHD_realData::mass));
                        p.rdata(FHD_realData::vely) *= (tRatio*cellAvMass/p.rdata(FHD_realData::mass));
                        p.rdata(FHD_realData::velz) *= (tRatio*cellAvMass/p.rdata(FHD_realData::mass));

                    }

                }


            }
            }
            }
        }
    }
