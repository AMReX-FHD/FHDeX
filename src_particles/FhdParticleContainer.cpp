#include "FhdParticleContainer.H"
#include <filesystem>
#include "particle_functions_K.H"
#include "paramplane_functions_K.H"
#include <math.h>

bool FhdParticleContainer::use_neighbor_list  {true};
bool FhdParticleContainer::sort_neighbor_list {false};

FhdParticleContainer::FhdParticleContainer(const Geometry & geom,
                                           const Geometry & geomF,
                                           const DistributionMapping & dmap,
                                           const BoxArray & ba,
                                           const BoxArray & baF,
                                           int n_nbhd,
                                           int ngF)
    : IBMarkerContainerBase<FHD_realData, FHD_intData>(geom, geomF, dmap, ba, baF, n_nbhd, ngF), n_list(0)
{
    BL_PROFILE_VAR("FhdParticleContainer()",FhdParticleContainer);

    InitInternals(n_nbhd);
    nghost = n_nbhd;

    double domx, domy, domz;

    domx = (prob_hi[0] - prob_lo[0]);
    domy = (prob_hi[1] - prob_lo[1]);
    domz = (prob_hi[2] - prob_lo[2]);

    if (radialdist_int > 0 || cartdist_int > 0) {
        if (binSize == 0. || searchDist == 0.) {
            Abort("Must set binSize and searchDist if computing g(r)");
        }
    }

    if (searchDist > 0.5*domx || searchDist > 0.5*domy || searchDist > 0.5*domz) {
        Abort("searchDist is greater than half the domain length");
    }

    if (radialdist_int > 0 || cartdist_int > 0) {

        // create enough bins to look within a sphere with radius equal to "half" of the domain
        totalBins = (int)amrex::Math::floor((searchDist)/((double)binSize)) - 1;

        Print() << "Bin size for pair correlation function: " << binSize << std::endl;
        Print() << "Number of pair correlation bins: " << totalBins << std::endl;

        // compute the volume of each bin
        binVolRadial = new Real[totalBins]();
        for(int i=0;i<totalBins;i++) {
            binVolRadial[i]= (4.0/3.0)*3.14159265359*(pow((i+1)*binSize,3) - pow((i)*binSize,3));
        }

        binVolCartesian =  2.*binSize * 2.*searchDist * 2.*searchDist;

        // how many snapshots
        radialStatsCount = 0;
        cartesianStatsCount = 0;
    }

    // storage for mean nearest neighbour counts
    if (radialdist_int > 0) {
        nearestN    = new Real[nspecies*nspecies + 1]();
        for(int i=0;i<nspecies*nspecies + 1;i++)
        {
            nearestN[i] = 0;
        }
    }

    // storage for mean radial distribution
    if (radialdist_int > 0) {
        meanRadialDistribution    = new Real[totalBins]();
        meanRadialDistribution_pp = new Real[totalBins]();
        meanRadialDistribution_pm = new Real[totalBins]();
        meanRadialDistribution_mm = new Real[totalBins]();
    }

    // storage for mean Cartesian distributions
    if (cartdist_int > 0) {
        meanXDistribution    = new Real[totalBins]();
        meanXDistribution_pp = new Real[totalBins]();
        meanXDistribution_pm = new Real[totalBins]();
        meanXDistribution_mm = new Real[totalBins]();

        meanYDistribution    = new Real[totalBins]();
        meanYDistribution_pp = new Real[totalBins]();
        meanYDistribution_pm = new Real[totalBins]();
        meanYDistribution_mm = new Real[totalBins]();

        meanZDistribution    = new Real[totalBins]();
        meanZDistribution_pp = new Real[totalBins]();
        meanZDistribution_pm = new Real[totalBins]();
        meanZDistribution_mm = new Real[totalBins]();
    }

    //Remove files that we will be appending to.
    remove("msdEst");
    for(int i=0;i<nspecies;i++)
    {
        string filename = Concatenate("msdEst_",i);
        //        int n = filename.lenght();
        //        char char_array[n+1];
        //        strcpy
        remove(filename.c_str());
    }

    remove("velOut");
    remove("matOut");
    remove("conductivityEst");
    remove("currentEst");
    remove("nearestNeighbour");

    remove("part1X");
    remove("part1Y");
    remove("part1Z");
    remove("part2X");
    remove("part2Y");
    remove("part2Z");

    Real dr = threepmRange/threepmBins;

    threepmVals[0] = 0;
    threepmMin[0] = 0;
    threepmMax[0] = 0;
    threepmPoints[0] = 0;

    if(sr_tog == 4)
    {
        std::ifstream bottomFile("bottomWall.dat");
        std::ifstream topFile("topWall.dat");

        if(bottomFile.good() && topFile.good())
        {
            bottomFile >> bottomListLength;
            bottomList = new Triplet[bottomListLength];
            for(int i=0;i<bottomListLength;i++)
            {
                bottomFile >> bottomList[i].x;
                bottomFile >> bottomList[i].y;
                bottomFile >> bottomList[i].z;

            }

            topFile >> topListLength;
            topList = new Triplet[topListLength];
            for(int i=0;i<topListLength;i++)
            {
                topFile >> topList[i].x;
                topFile >> topList[i].y;
                topFile >> topList[i].z;
            }

            bottomFile.close();
            topFile.close();
            Print() << "Top/bottom particle files read.\n";

        }else
        {
            Abort("Couldn't read bottom/top particle input file!\n");
        }

        bottomFile.close();
        topFile.close();
    }

    for(int i=1;i<threepmBins;i++)
    {
       threepmPoints[i] = i*dr;

       threepmMax[i] = 0;
       threepmMin[i] = 10000000;
    }

    doRedist = 1;
}

void FhdParticleContainer::forceFunction(Real dt)
{
   const int lev = 0;

   //Real k = 5e2;
   Real maxU = 0;
   Real maxD = 0;
   int pinCheck = 0;

   for (FhdParIter pti(*this, lev, MFItInfo().SetDynamic(false)); pti.isValid(); ++pti) {

        PairIndex index(pti.index(), pti.LocalTileIndex());
        AoS& particles = pti.GetArrayOfStructs();
        int np = pti.numParticles();

        const Box& tile_box  = pti.tilebox();

        Real maxUtile = 0;
        Real maxDtile = 0;

         for (int i = 0; i < np; ++i) {

            ParticleType & part = particles[i];

            if(part.rdata(FHD_realData::spring) != 0)
            {
                Real radVec[3];
                radVec[0] = part.pos(0);
                radVec[1] = part.pos(1);
                radVec[2] = part.pos(2);

                //radVec[0] = part.rdata(FHD_realData::ax);
                //radVec[1] = part.rdata(FHD_realData::ay);
                //radVec[2] = part.rdata(FHD_realData::az);

                //Real kFac = 6*M_PI*part.rdata(FHD_realData::radius)*visc_coef/dt;

                //Print() << "k: " << kFac << endl;

                part.rdata(FHD_realData::forcex) = part.rdata(FHD_realData::forcex) - part.rdata(FHD_realData::spring)*radVec[0];
                part.rdata(FHD_realData::forcey) = part.rdata(FHD_realData::forcey) - part.rdata(FHD_realData::spring)*radVec[1];
                part.rdata(FHD_realData::forcez) = part.rdata(FHD_realData::forcez) - part.rdata(FHD_realData::spring)*radVec[2];

                //                part.rdata(FHD_realData::forcex) = part.rdata(FHD_realData::forcex) - kFac*radVec[0];
                //                part.rdata(FHD_realData::forcey) = part.rdata(FHD_realData::forcey) - kFac*radVec[1];
                //                part.rdata(FHD_realData::forcez) = part.rdata(FHD_realData::forcez) - kFac*radVec[2];

                Real dSqr = (pow(radVec[0],2) + pow(radVec[1],2) + pow(radVec[2],2));

                part.rdata(FHD_realData::potential) = 0.5*part.rdata(FHD_realData::spring)*dSqr;

                if(part.rdata(FHD_realData::potential) > maxUtile)
                {
                    maxUtile = part.rdata(FHD_realData::potential);
                }

                if((dSqr/(0.5*part.rdata(FHD_realData::sigma))) > maxDtile)
                {
                    maxDtile = dSqr/(0.5*part.rdata(FHD_realData::sigma));
                }
                pinCheck = 1;
                //Print() << radVec[0] << endl;
            }
         }
        maxU = amrex::max(maxUtile, maxU);
        maxD = amrex::max(maxDtile, maxD);
    }

    ParallelDescriptor::ReduceRealMax(maxU);
    ParallelDescriptor::ReduceRealMax(maxD);
    ParallelDescriptor::ReduceIntMax(pinCheck);
    //Print() << "Max potential: " << maxU << std::endl;
    if(pinCheck != 0)
    {
        Print() << "Maximum observed pinned particle displacement (fraction of radius): " << sqrt(maxD) << std::endl;
    }
}

void FhdParticleContainer::pinForce()
{
   const int lev = 0;

   for (FhdParIter pti(*this, lev, MFItInfo().SetDynamic(false)); pti.isValid(); ++pti) {

        PairIndex index(pti.index(), pti.LocalTileIndex());
        AoS& particles = pti.GetArrayOfStructs();
        int np = pti.numParticles();

        const Box& tile_box  = pti.tilebox();

        Real maxUtile = 0;
        Real maxDtile = 0;

         for (int i = 0; i < np; ++i) {

            ParticleType & part = particles[i];

            if(part.idata(FHD_intData::pinned) != 0)
            {
                part.rdata(FHD_realData::forcex) += -1.0*k_B*T_init[0]*part.rdata(FHD_realData::velx)/part.rdata(FHD_realData::wetDiff);
                part.rdata(FHD_realData::forcey) += -1.0*k_B*T_init[0]*part.rdata(FHD_realData::vely)/part.rdata(FHD_realData::wetDiff);
                part.rdata(FHD_realData::forcez) += -1.0*k_B*T_init[0]*part.rdata(FHD_realData::velz)/part.rdata(FHD_realData::wetDiff);

                Print() << part.rdata(FHD_realData::forcex) << ", " << part.rdata(FHD_realData::velx) << endl;
            }
         }
    }
}

void FhdParticleContainer::velNorm()
{
   const int lev = 0;

   //Real k = 5e2;
   Real velN = 0;
   Real maxS = 0;
   int pinCheck;

   for (FhdParIter pti(*this, lev, MFItInfo().SetDynamic(false)); pti.isValid(); ++pti) {

        PairIndex index(pti.index(), pti.LocalTileIndex());
        AoS& particles = pti.GetArrayOfStructs();
        int np = pti.numParticles();

        const Box& tile_box  = pti.tilebox();

        Real velNtile = 0;
        Real maxStile = 0;

         for (int i = 0; i < np; ++i) {

            ParticleType & part = particles[i];

            if(part.idata(FHD_intData::pinned) != 0)
            {

                Real velSqr = (pow(part.rdata(FHD_realData::velx),2) + pow(part.rdata(FHD_realData::vely),2) + pow(part.rdata(FHD_realData::velz),2));

                velNtile += velSqr;
                //Print() << "actual vel: " << part.rdata(FHD_realData::velx) << ", " << part.rdata(FHD_realData::vely) << ", " << part.rdata(FHD_realData::velz) << std::endl;

                if(velSqr > maxStile)
                {
                    maxStile = velSqr;
                }
                pinCheck = 1;
            }
         }
        maxS = amrex::max(maxStile, maxS);
        velN += velNtile;
    }

    ParallelDescriptor::ReduceRealMax(maxS);
    ParallelDescriptor::ReduceRealSum(velN);

    if(pinCheck != 0)
    {
        Print() << setprecision(15) << "Pinned particle velocity norm: " << sqrt(velN) << ", max speed: " << sqrt(maxS) << std::endl;
    }
}

void FhdParticleContainer::computeForcesNLGPU(const MultiFab& charge, const MultiFab& coords, const Real* dx) {

    BL_PROFILE_VAR("computeForcesNL()",computeForcesNL);

    Real rcount = 0;
    Real rdcount = 0;
    Real recount = 0;
    Real recountI = 0;
    const int lev = 0;

#ifdef _OPENMP
#pragma omp parallel
#endif

    if(doRedist != 0)
    {
        fillNeighbors();

        buildNeighborList(CHECK_PAIR{});
    }

   for (FhdParIter pti(*this, lev, MFItInfo().SetDynamic(false)); pti.isValid(); ++pti)
   {
        PairIndex index(pti.index(), pti.LocalTileIndex());
        AoS& particles = pti.GetArrayOfStructs();
        int Np = pti.numParticles();
        int Nn = pti.numNeighborParticles();
        int size = neighbor_list[lev][index].size();

        const Box& tile_box  = pti.tilebox();

        if (sr_tog != 0)
        {

            compute_forces_nl_gpu(particles, Np, Nn,
                              m_neighbor_list[lev][index], topList, bottomList, topListLength, bottomListLength, rcount, rdcount);
        }

        if (es_tog==3)
        {
            compute_p3m_sr_correction_nl_gpu(particles, Np, Nn,
                                        m_neighbor_list[lev][index], dx, recount, recountI);
        }
    }

    if(sr_tog != 0)
    {
            ParallelDescriptor::ReduceRealSum(rcount);
            ParallelDescriptor::ReduceRealSum(rdcount);

            Print() << rcount/2 << " close range interactions.\n";
            Print() << rdcount << " wall interactions.\n";
    }
    if(es_tog==3)
    {
            ParallelDescriptor::ReduceRealSum(recount);
            ParallelDescriptor::ReduceRealSum(recountI);

            Print() << recount/2 << " p3m interactions.\n";
            Print() << recountI << " image charge interactions.\n";
    }
}

void FhdParticleContainer::computeForcesCoulombGPU(long totalParticles) {

    BL_PROFILE_VAR("computeForcesCoulomb()",computeForcesCoulomb);

    using namespace amrex;
    using common::permittivity;
    using common::images;

    Real ee = (1.0/(permittivity*4*3.14159265));

    GpuArray<Real, 3> plo = {prob_lo[0], prob_lo[1], prob_lo[2]};
    GpuArray<Real, 3> phi = {prob_hi[0], prob_hi[1], prob_hi[2]};

    const int lev = 0;
    double domx, domy, domz;

    domx = (phi[0] - plo[0]);
    domy = (phi[1] - plo[1]);
    domz = (phi[2] - plo[2]);

//    Real posx[totalParticles];
    Gpu::ManagedVector<Real> posx;
    posx.resize(totalParticles);
    Real * posxPtr = posx.dataPtr();

    Gpu::ManagedVector<Real> posy;
    posy.resize(totalParticles);
    Real * posyPtr = posy.dataPtr();

    Gpu::ManagedVector<Real> posz;
    posz.resize(totalParticles);
    Real * poszPtr = posz.dataPtr();

    Gpu::ManagedVector<Real> charge;
    charge.resize(totalParticles);
    Real * chargePtr = charge.dataPtr();

//    Real posx[totalParticles];
//    Real posy[totalParticles];
//    Real posz[totalParticles];
//    int  species[totalParticles];
//    Real charge[totalParticles];

    Print() << "Calculating Coulomb force for each particle pair\n";

    // collect particle positions onto one processor
    PullDown(0, posxPtr, -1, totalParticles);
    PullDown(0, posyPtr, -2, totalParticles);
    PullDown(0, poszPtr, -3, totalParticles);
    PullDown(0, chargePtr, FHD_realData::q , totalParticles);

    int imag = (images == 0) ? 1 : images;

    Real maxdist = 0.99*amrex::min(imag * domx,
                                   //imag * domy,
                                   imag * domz);

    for (FhdParIter pti(*this, lev); pti.isValid(); ++pti) {

        const int grid_id = pti.index();
        const int tile_id = pti.LocalTileIndex();

        auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];
        auto& particles = particle_tile.GetArrayOfStructs();
        const int np = particles.numParticles();

       auto pstruct = particles().dataPtr();

        // loop over particles
//                for(int i = 0; i < np; i++)
        AMREX_FOR_1D( np, i,
        {

            ParticleType & part = pstruct[i];

            double dr2;
            double rtdr2;
            double dx;
            double dy;
            double dz;

            Real q1 = part.rdata(FHD_realData::q);

            for(int j = 0; j < totalParticles; j++)
            {

                 Real q2 = chargePtr[j];

                // (currently hard-coded for y-wall)
                for(int ii = -images; ii <= images; ii++)
                {
                       for(int kk = -images; kk <= images; kk++)
                       {

                          // get distance between particles
                          dx = part.pos(0)-posxPtr[j] - ii*domx;
                          dy = part.pos(1)-posyPtr[j];// - jj*domy;
                          dz = part.pos(2)-poszPtr[j] - kk*domz;

                          dr2 = dx*dx + dy*dy + dz*dz;
                          rtdr2 = sqrt(dr2);

                            if (rtdr2 < maxdist && rtdr2 > 0.0)
                          {
                              part.rdata(FHD_realData::forcex) += ee*(dx/rtdr2)*q1*q2/dr2;
                              part.rdata(FHD_realData::forcey) += ee*(dy/rtdr2)*q1*q2/dr2;
                              part.rdata(FHD_realData::forcez) += ee*(dz/rtdr2)*q1*q2/dr2;
                          }
                       }
                }
            }
        });
    }

    Print() << "Finished Coulomb calculation\n";
}

void FhdParticleContainer::MoveIonsCPP(const Real dt, const Real* dxFluid, const Real* dxE, const Geometry geomF,
                                    const std::array<MultiFab, AMREX_SPACEDIM>& umac, const std::array<MultiFab, AMREX_SPACEDIM>& efield,
                                    const std::array<MultiFab, AMREX_SPACEDIM>& RealFaceCoords,
                                    std::array<MultiFab, AMREX_SPACEDIM>& source,
                                    std::array<MultiFab, AMREX_SPACEDIM>& sourceTemp,
                                    paramPlane* paramPlaneList, const int paramPlaneCount, int sw)
{
    BL_PROFILE_VAR("MoveIons()",MoveIons);

    //AMREX_GPU_ERROR_CHECK();
    const int lev = 0;
    const GpuArray<Real, 3> dx = Geom(lev).CellSizeArray();
    const GpuArray<Real, 3> plo = Geom(lev).ProbLoArray();
    const GpuArray<Real, 3> phi = Geom(lev).ProbHiArray();

    Gpu::ManagedVector<Real> domsize(3);
    Real* pdomsize = domsize.data();

    for (int d=0; d<AMREX_SPACEDIM; ++d)
    {
        pdomsize[d] = phi[d]-plo[d];
    }

    double kinetic = 0;

    doRedist = 0;

    int        np_tile = 0 ,       np_proc = 0 ; // particle count
    Real rejected_tile = 0., rejected_proc = 0.; // rejected moves in midpoint scheme
    Real    moves_tile = 0.,    moves_proc = 0.; // total moves in midpoint scheme
    Real maxspeed_tile = 0., maxspeed_proc = 0.; // max speed
    Real  maxdist_tile = 0.,  maxdist_proc = 0.; // max displacement (fraction of radius)
    Real diffinst_tile = 0., diffinst_proc = 0.; // average diffusion coefficient

    Real adj = 0.99999;
    Real adjalt = 2.0*(1.0-0.99999);
    //Real runtime, inttime;
    //int intsurf, intside, push;
    //int midpoint = 0;
    //Real posAlt[3];
    Real check;

    if(all_dry != 1)
    {
    InterpolateMarkersGpu(0, dxFluid, umac, RealFaceCoords,check);
    if(move_tog == 2)
    {
        //// Set up reducing operation across gpu (instead of ParallelFor)
        //ReduceOps<ReduceOpSum> reduce_op;
        //ReduceData<int> reduce_data(reduce_op);
        //using ReduceTuple = typename decltype(reduce_data)::Type;

        for (MyIBMarIter pti(* this, lev); pti.isValid(); ++pti) {

            TileIndex index(pti.index(), pti.LocalTileIndex());

            AoS & aos = this->GetParticles(lev).at(index).GetArrayOfStructs();
        ParticleType* particles = aos().dataPtr();
            long np = this->GetParticles(lev).at(index).numParticles();

            // Set up vector to do reduction
        Gpu::DeviceVector<Real> increment_moves_tile(np, 0.);
            Real* pincrement_moves_tile = increment_moves_tile.data();

            //Real posOld[3];
            //Real velOld[3];
            //moves_tile = 0;

        //// Using this loop to set up counter via a reduce sum operation across each thread
        //reduce_op.eval(np, reduce_data, [=] AMREX_GPU_DEVICE (int i) -> ReduceTuple

        // Set up RNG engine with ParallelForRNG, and do reduction using a np-sized vector storing value for each particle
        amrex::ParallelForRNG(np, [=] AMREX_GPU_DEVICE (int i, amrex::RandomEngine const& engine) noexcept
        //for (int i = 0; i < np; ++ i)
        {

                ParticleType & part = particles[i];

                GpuArray<Real, 3> posOld;
                GpuArray<Real, 3> velOld;
                GpuArray<Real, 3> posAlt;

                if(part.idata(FHD_intData::pinned) == 0)
                {
                        //moves_tile++;
                        pincrement_moves_tile[i] = 1.;

                        for (int d=0; d<AMREX_SPACEDIM; ++d)
                        {
                            velOld[d] = part.rdata(FHD_realData::velx + d);
                            part.rdata(FHD_realData::pred_posx + d) = part.pos(d);
                        }

                        Real runtime = 0.5*dt;
                        Real inttime = 0;
                        int midpoint = 0;
                        int intsurf, intside, push;

                        while(runtime > 0)
                        {
                            //find_inter(&part, &runtime, paramPlaneList, &paramPlaneCount, &intsurf, &inttime, &intside, AMREX_ZFILL(plo), AMREX_ZFILL(phi));
                            find_inter_gpu(part, runtime, paramPlaneList, paramPlaneCount, &intsurf, &inttime, &intside, AMREX_ZFILL(plo), AMREX_ZFILL(phi));

                            for (int d=0; d<AMREX_SPACEDIM; ++d)
                            {
                                posAlt[d] = inttime * part.rdata(FHD_realData::velx + d)*adjalt;
                            }
                            for (int d=0; d<AMREX_SPACEDIM; ++d)
                            {
                                part.pos(d) += inttime * part.rdata(FHD_realData::velx + d)*adj;
                            }

                            if(intsurf > 0)
                            {

                                paramPlane& surf = paramPlaneList[intsurf-1]; //find_inter indexes from 1 to maintain compatablity with fortran version

                                if(surf.periodicity == 0)
                                {
                                   Real dummy = 1;
                                   //std::cout << "Pre: " << part.rdata(FHD_realData::velx) << ", " << part.rdata(FHD_realData::vely) << ", " << part.rdata(FHD_realData::velz) << ", " << intside << "\n";
                                   //app_bc(&surf, &part, &intside, domsize, &push, &dummy, &dummy);
                                   //app_bc_gpu(&surf, part, intside, pdomsize, &push, &runtime, dummy, engine);
                                   //std::cout << "Post: " << part.rdata(FHD_realData::velx) << ", " << part.rdata(FHD_realData::vely) << ", " << part.rdata(FHD_realData::velz) << ", " << intside << "\n";

                                   runtime = runtime - inttime;

                                }else
                                {
                                  runtime = runtime - inttime;

                                  for (int d=0; d<AMREX_SPACEDIM; ++d)
                                  {

                                    part.pos(d) += runtime * part.rdata(FHD_realData::velx + d);
                                  }
                                  runtime = 0;

                                }

                            }else
                            {
                               runtime = 0;

                            }
                        }
                }
            //return increment_moves_tile;
            });

            //// Get the value of the reduce sum on one processor
            //moves_tile = amrex::get<0>(reduce_data.value());

            //Gpu::synchronize();
            //AMREX_GPU_ERROR_CHECK();
            moves_proc = Reduce::Sum(np, pincrement_moves_tile);
            //Gpu::synchronize();
            //moves_proc    += moves_tile;
            //Gpu::synchronize();

            //std::cout << "Moves " << moves_tile << std::endl;

        }

        //Need to add midpoint rejecting feature here.
        InterpolateMarkersGpu(0, dxFluid, umac, RealFaceCoords, check);

        for (MyIBMarIter pti(* this, lev); pti.isValid(); ++pti) {

            TileIndex index(pti.index(), pti.LocalTileIndex());

        AoS & aos = this->GetParticles(lev).at(index).GetArrayOfStructs();
        ParticleType* particles = aos().dataPtr();
        long np = this->GetParticles(lev).at(index).numParticles();

        amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE (int i) noexcept
            //for (int i = 0; i < np; ++ i)
        {
                ParticleType & part = particles[i];
                if(part.idata(FHD_intData::pinned) == 0)
                {
                        for (int d=0; d<AMREX_SPACEDIM; ++d)
                        {
                            part.pos(d) = part.rdata(FHD_realData::pred_posx + d);
                        }
                }
            });
        }
    }
    }

    //cin.get();

    if((dry_move_tog == 1) || (dry_move_tog == 2))
    {
        for (MyIBMarIter pti(* this, lev); pti.isValid(); ++pti) {

            TileIndex index(pti.index(), pti.LocalTileIndex());

            AoS & aos = this->GetParticles(lev).at(index).GetArrayOfStructs();
            ParticleType* particles = aos().dataPtr();
            long np = this->GetParticles(lev).at(index).numParticles();

        amrex::ParallelForRNG(np, [=] AMREX_GPU_DEVICE (int i, amrex::RandomEngine const& engine) noexcept
            //for (int i = 0; i < np; ++ i)
            {
                ParticleType & part = particles[i];
                if(part.idata(FHD_intData::pinned) == 0)
                {
                        GpuArray<Real, 3> mb;
                        GpuArray<Real, 3> mbDer;
                        GpuArray<Real, 3> dry_terms;

                        get_explicit_mobility_gpu(mb, mbDer, part, plo, phi);

                        dry_gpu(dt, part,dry_terms, mb, mbDer, engine);

                        for (int d=0; d<AMREX_SPACEDIM; ++d)
                        {
                            part.rdata(FHD_realData::velx + d) += dry_terms[d];
                        }
                }
            });
        }
    }

    //Real maxspeed = 0;
    //Real maxdist = 0;
    //Real totaldist, diffest;
    //Real diffinst = 0;
    int moves = 0;
    int reDist = 0;

    //// Initialize 5 reduce operations, in the order of maxspeed, maxdist, diffest, moves, reDist
    //ReduceOps<ReduceOpMax,ReduceOpMax, ReduceOpSum, ReduceOpSum, ReduceOpSum> reduce_op5;
    //ReduceData<Real, Real, Real, int, int> reduce_data5(reduce_op5);
    //using ReduceTuple = typename decltype(reduce_data5)::Type;

    for (MyIBMarIter pti(* this, lev); pti.isValid(); ++pti) {

        TileIndex index(pti.index(), pti.LocalTileIndex());

        Box bx  = pti.tilebox();
        IntVect myLo = bx.smallEnd();
        IntVect myHi = bx.bigEnd();

        AoS & aos = this->GetParticles(lev).at(index).GetArrayOfStructs();
        ParticleType* particles = aos().dataPtr();
        long np = this->GetParticles(lev).at(index).numParticles();

        np_proc += np;

        // Set up vectors to do reduction
        Gpu::DeviceVector<int> increment_moves(np, 0);
        Gpu::DeviceVector<int> increment_reDist(np, 0);
        Gpu::DeviceVector<Real> increment_maxspeed(np, 0.);
        Gpu::DeviceVector<Real> increment_maxdist(np, 0.);
        Gpu::DeviceVector<Real> increment_diffest(np, 0.);
        int* pincrement_moves = increment_moves.data();
        int* pincrement_reDist = increment_reDist.data();
        Real* pincrement_maxspeed = increment_maxspeed.data();
        Real* pincrement_maxdist = increment_maxdist.data();
        Real* pincrement_diffest = increment_diffest.data();

        //reduce_op5.eval(np, reduce_data5, [=] AMREX_GPU_DEVICE (int i) -> ReduceTuple

        // Set up RNG engine with ParallelForRNG, and do reduction using a np-sized vector storing value for each particle
        amrex::ParallelForRNG(np, [=] AMREX_GPU_DEVICE (int i, amrex::RandomEngine const& engine) noexcept
        //for (int i = 0; i < np; ++ i)
        {

            ParticleType & part = particles[i];

            GpuArray<Real, 3> posAlt;

            if(part.idata(FHD_intData::pinned) == 0)
            {
                Real speed = 0;

                for (int d=0; d<AMREX_SPACEDIM; ++d)
                {
                    speed += part.rdata(FHD_realData::velx + d)*part.rdata(FHD_realData::velx + d);
                }

                if(speed > pincrement_maxspeed[i])
                {
                    pincrement_maxspeed[i] = speed;
                }

                pincrement_moves[i] = 1;

                Real runtime = dt;
                Real inttime = 0;
                int intsurf, intside, push;

                //Real thisMove[3] = {0,0,0};

                while(runtime > 0)
                {
                    //find_inter(&part, &runtime, paramPlaneList, &paramPlaneCount, &intsurf, &inttime, &intside, AMREX_ZFILL(plo), AMREX_ZFILL(phi));
                    find_inter_gpu(part, runtime, paramPlaneList, paramPlaneCount, &intsurf, &inttime, &intside, AMREX_ZFILL(plo), AMREX_ZFILL(phi));
                    //Print() << "PART " << part.id() << ", " << intsurf << "\n";
                    //cin.get();

                    for (int d=0; d<AMREX_SPACEDIM; ++d)
                    {
                        posAlt[d] = inttime * part.rdata(FHD_realData::velx + d)*adjalt;
                    }
                    for (int d=0; d<AMREX_SPACEDIM; ++d)
                    {
                        part.pos(d) += inttime * part.rdata(FHD_realData::velx + d)*adj;
                        part.rdata(FHD_realData::ax +d ) += inttime * part.rdata(FHD_realData::velx + d)*adj;
                        //thisMove[d] += inttime * part.rdata(FHD_realData::velx + d)*adj;
                    }
                    runtime = runtime - inttime;
                    if(intsurf > 0)
                    {
                        paramPlane& surf = paramPlaneList[intsurf-1];//find_inter indexes from 1 to maintain compatablity with fortran version

                        Real dummy = 1;
                        //app_bc(&surf, &part, &intside, domsize, &push, &dummy, &dummy);

                        app_bc_gpu(&surf, part, intside, pdomsize, &push, &runtime, dummy, engine);

                        //std::cout << "Post: " << part.id() << ", " << part.rdata(FHD_realData::velx) << ", " << part.rdata(FHD_realData::vely) << ", " << part.rdata(FHD_realData::velz) << ", " << intsurf << "\n";

                        if(push == 1)
                        {
                            for (int d=0; d<AMREX_SPACEDIM; ++d)
                            {
                                part.pos(d) += part.pos(d) + posAlt[d];
                                part.rdata(FHD_realData::ax + d) += part.rdata(FHD_realData::ax + d) + posAlt[d];
                                //thisMove[d] += part.pos(d) + posAlt[d];
                            }
                        }
                    }
                }

//                for (int d=0; d<AMREX_SPACEDIM; ++d)
//                {
//                    part.rdata(FHD_realData::ax + d) += part.rdata(FHD_realData::velx + d)*dt;
//                }

    //            Print() << part.id() << " vel: " << setprecision(15) << part.rdata(FHD_realData::velx) << " pos: " << part.pos(0) << "\n";

                Real dist = dt*sqrt(part.rdata(FHD_realData::velx)*part.rdata(FHD_realData::velx) + part.rdata(FHD_realData::vely)*part.rdata(FHD_realData::vely) + part.rdata(FHD_realData::velz)*part.rdata(FHD_realData::velz))/part.rdata(FHD_realData::radius);

                Real totaldist = sqrt(part.rdata(FHD_realData::ax)*part.rdata(FHD_realData::ax) + part.rdata(FHD_realData::ay)*part.rdata(FHD_realData::ay) + part.rdata(FHD_realData::az)*part.rdata(FHD_realData::az));

                if(dist > pincrement_maxdist[i])
                {
                    pincrement_maxdist[i] = dist;
                }

                //std::cout << "MAXDIST: " << maxdist << "\n";

                part.rdata(FHD_realData::travelTime) += dt;

                pincrement_diffest[i] = totaldist/(6.0*part.rdata(FHD_realData::travelTime));

                //diffinst += diffest;
            }

            GpuArray<int, 3> cell;
            cell[0] = (int)floor((part.pos(0)-plo[0])/dx[0]);
            cell[1] = (int)floor((part.pos(1)-plo[1])/dx[1]);
            cell[2] = (int)floor((part.pos(2)-plo[2])/dx[2]);

            if((cell[0] < myLo[0]) || (cell[1] < myLo[1]) || (cell[2] < myLo[2]) || (cell[0] > myHi[0]) || (cell[1] > myHi[1]) || (cell[2] > myHi[2]))
            {
                pincrement_reDist[i] = 1;
            }

            //// For ReduceOP, return reduced data at the end.
            //return { maxspeed, maxdist, diffest, increment_moves, increment_reDist };
        });

        //// These are ReduceOp-way of getting the value of reduced data
        ////   Need to update to these once there is a RNG-engine feature of ReduceOp
        //maxspeed_proc = amrex::max(maxspeed_proc, amrex::get<0>(reduce_data5.value()));
        //maxdist_proc = amrex::max(maxdist_proc, amrex::get<1>(reduce_data5.value()));
        //diffinst_proc = amrex::get<2>(reduce_data5.value());
        //moves += amrex::get<3>(reduce_data5.value());
        //reDist += amrex::get<4>(reduce_data5.value());

        moves = Reduce::Sum(np, pincrement_moves);
        reDist = Reduce::Sum(np, pincrement_reDist);
        maxspeed_proc = amrex::max(maxspeed_proc, Reduce::Max(np, pincrement_maxspeed));
        maxdist_proc  = amrex::max(maxdist_proc, Reduce::Max(np, pincrement_maxdist));
        //std::cout << "MAXDISTPROC: " << maxdist_proc << "\n";

        diffinst_proc += Reduce::Sum(np, pincrement_diffest);
    }

    // gather statistics
    ParallelDescriptor::ReduceIntSum(np_proc);
    ParallelDescriptor::ReduceRealSum(check);
    ParallelDescriptor::ReduceRealSum(moves_proc);
    ParallelDescriptor::ReduceIntSum(moves);
    ParallelDescriptor::ReduceRealMax(maxspeed_proc);
    ParallelDescriptor::ReduceRealMax(maxdist_proc);
    ParallelDescriptor::ReduceRealSum(diffinst_proc);
    ParallelDescriptor::ReduceIntSum(reDist);

    // write out global diagnostics
    if (ParallelDescriptor::IOProcessor()) {
        Print() << "I see " << np_proc << " particles\n";
        if (move_tog == 2 && all_dry != 1) {
            Print() << "Fraction of midpoint moves rejected: " << check/moves_proc << "\n";
        }
        Print() << reDist << " particles to be redistributed.\n";
        Print() << "Maximum observed speed: " << sqrt(maxspeed_proc) << "\n";
        Print() << "Maximum observed displacement (fraction of radius): " << maxdist_proc << "\n";
        //Print() << "Average diffusion coefficient: " << diffinst_proc/np_proc << "\n";
    }
    //    if(reDist != 0)
    //    {
        Redistribute();
        doRedist = 1;
    //    }
}

void FhdParticleContainer::SpreadIonsGPU(const Real* dxFluid, const Real* dxE, const Geometry geomF,
                                      const std::array<MultiFab, AMREX_SPACEDIM>& umac,
                                      const std::array<MultiFab, AMREX_SPACEDIM>& coords,
                                      std::array<MultiFab, AMREX_SPACEDIM>& efield,
                                      std::array<MultiFab, AMREX_SPACEDIM>& source,
                                      std::array<MultiFab, AMREX_SPACEDIM>& sourceTemp)
{
    BL_PROFILE_VAR("SpreadIons()",SpreadIons);

    const int lev = 0;
    const Real* dx = Geom(lev).CellSize();
    const Real* plo = Geom(lev).ProbLo();
    const Real* phi = Geom(lev).ProbHi();

    double potential = 0;

#ifdef _OPENMP
#pragma omp parallel
#endif

    for (FhdParIter pti(*this, lev); pti.isValid(); ++pti)
    {
        const int grid_id = pti.index();
        const int tile_id = pti.LocalTileIndex();
        const Box& tile_box  = pti.tilebox();

        auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];
        auto& particles = particle_tile.GetArrayOfStructs();
        const int np = particles.numParticles();

        emf_gpu(particles,
                      efield[0][pti], efield[1][pti], efield[2][pti],
                      AMREX_ZFILL(plo), AMREX_ZFILL(dxE));
    }

    if(fluid_tog != 0)
    {
        //spread_ions_fhd_gpu(particles,
        //                 sourceTemp[0][pti], sourceTemp[1][pti], sourceTemp[2][pti],
        //                 AMREX_ZFILL(plo),
        //                 AMREX_ZFILL(dxFluid));

        SpreadMarkersGpu(lev, sourceTemp, coords, dxFluid, 1);
    }

    if(fluid_tog != 0)
    {
        for (int i=0; i<AMREX_SPACEDIM; ++i) {
            MultiFabPhysBCDomainStress(sourceTemp[i], geomF, i);
            MultiFabPhysBCMacStress(sourceTemp[i], geomF, i);
        }

        sourceTemp[0].SumBoundary(geomF.periodicity());
        sourceTemp[1].SumBoundary(geomF.periodicity());
#if (AMREX_SPACEDIM == 3)
        sourceTemp[2].SumBoundary(geomF.periodicity());
#endif

        MultiFab::Add(source[0],sourceTemp[0],0,0,source[0].nComp(),source[0].nGrow());
        MultiFab::Add(source[1],sourceTemp[1],0,0,source[1].nComp(),source[1].nGrow());
#if (AMREX_SPACEDIM == 3)
        MultiFab::Add(source[2],sourceTemp[2],0,0,source[2].nComp(),source[2].nGrow());
#endif

        source[0].FillBoundary(geomF.periodicity());
        source[1].FillBoundary(geomF.periodicity());
#if (AMREX_SPACEDIM == 3)
        source[2].FillBoundary(geomF.periodicity());
#endif
    }
}

void FhdParticleContainer::SpreadIonsGPU(const Real* dxFluid, const Geometry geomF,
                                      const std::array<MultiFab, AMREX_SPACEDIM>& umac,
                                      const std::array<MultiFab, AMREX_SPACEDIM>& coords,
                                      std::array<MultiFab, AMREX_SPACEDIM>& source,
                                      std::array<MultiFab, AMREX_SPACEDIM>& sourceTemp)
{
    BL_PROFILE_VAR("SpreadIons()",SpreadIons);

    const int lev = 0;
    const Real* dx = Geom(lev).CellSize();
    const Real* plo = Geom(lev).ProbLo();
    const Real* phi = Geom(lev).ProbHi();

#ifdef _OPENMP
#pragma omp parallel
#endif

    //for (FhdParIter pti(* this, lev); pti.isValid(); ++pti)
    //{
    //    const int grid_id = pti.index();
    //    const int tile_id = pti.LocalTileIndex();
    //    const Box& tile_box  = pti.tilebox();
    //
    //    auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];
    //    auto& particles = particle_tile.GetArrayOfStructs();
    //    const int np = particles.numParticles();

        if(fluid_tog != 0)
        {
            //spread_ions_fhd_gpu(particles,
            //                 sourceTemp[0][pti], sourceTemp[1][pti], sourceTemp[2][pti],
            //                 AMREX_ZFILL(plo),
            //                 AMREX_ZFILL(dxFluid));
            SpreadMarkersGpu(lev, sourceTemp, coords, dxFluid, 1);
        }
    //}

    if(fluid_tog != 0)
    {
        for (int i=0; i<AMREX_SPACEDIM; ++i) {
            MultiFabPhysBCDomainStress(sourceTemp[i], geomF, i);
            MultiFabPhysBCMacStress(sourceTemp[i], geomF, i);
        }

        sourceTemp[0].SumBoundary(geomF.periodicity());
        sourceTemp[1].SumBoundary(geomF.periodicity());
#if (AMREX_SPACEDIM == 3)
        sourceTemp[2].SumBoundary(geomF.periodicity());
#endif

        MultiFab::Add(source[0],sourceTemp[0],0,0,source[0].nComp(),source[0].nGrow());
        MultiFab::Add(source[1],sourceTemp[1],0,0,source[1].nComp(),source[1].nGrow());
#if (AMREX_SPACEDIM == 3)
        MultiFab::Add(source[2],sourceTemp[2],0,0,source[2].nComp(),source[2].nGrow());
#endif

        source[0].FillBoundary(geomF.periodicity());
        source[1].FillBoundary(geomF.periodicity());
#if (AMREX_SPACEDIM == 3)
        source[2].FillBoundary(geomF.periodicity());
#endif
    }
}

void FhdParticleContainer::RadialDistribution(long totalParticles, const int step, const species* particleInfo)
{
    const int lev = 0;
    int bin;
    double domx, domy, domz, totalDist, temp;

    domx = (prob_hi[0] - prob_lo[0]);
    domy = (prob_hi[1] - prob_lo[1]);
    domz = (prob_hi[2] - prob_lo[2]);

    Real posx[totalParticles];
    Real posy[totalParticles];
    Real posz[totalParticles];
    int  species[totalParticles];

    Real charge[totalParticles];

    Print() << "Calculating radial distribution\n";

    // collect particle positions onto one processor
    PullDown(0, posx, -1, totalParticles);
    Print() << "HERE1\n";
    PullDown(0, posy, -2, totalParticles);
    Print() << "HERE2\n";
    PullDown(0, posz, -3, totalParticles);
    Print() << "HERE3\n";
    PullDown(0, charge, 27, totalParticles);
    Print() << "HERE4\n";
    PullDownInt(0, species, 4, totalParticles);
    Print() << "HERE5\n";

    // outer radial extent
    totalDist = totalBins*binSize;

    // this is the bin "hit count"
    RealVector radDist   (totalBins, 0.);
    RealVector radDist_pp(totalBins, 0.);
    RealVector radDist_pm(totalBins, 0.);
    RealVector radDist_mm(totalBins, 0.);

    double nearest[totalParticles*(nspecies+1)];
    double nn[nspecies*nspecies+1];

    for(int k = 0; k < (nspecies+1)*totalParticles; k++)
    {
        nearest[k] = 0;
    }

    for(int k = 0; k < nspecies*nspecies+1; k++)
    {
        nn[k] = 0;
    }

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (FhdParIter pti(*this, lev); pti.isValid(); ++pti) {

        const int grid_id = pti.index();
        const int tile_id = pti.LocalTileIndex();

        auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];
        auto& particles = particle_tile.GetArrayOfStructs();
        const int np = particles.numParticles();

        // loop over particles
        for (int i = 0; i < np; ++i) {

            ParticleType & part = particles[i];

            int iilo = (part.pos(0)-searchDist <= prob_lo[0]) ? -1 : 0;
            int iihi = (part.pos(0)+searchDist >= prob_hi[0]) ?  1 : 0;

            int jjlo = (part.pos(1)-searchDist <= prob_lo[1]) ? -1 : 0;
            int jjhi = (part.pos(1)+searchDist >= prob_hi[1]) ?  1 : 0;

            int kklo = (part.pos(2)-searchDist <= prob_lo[2]) ? -1 : 0;
            int kkhi = (part.pos(2)+searchDist >= prob_hi[2]) ?  1 : 0;

            double rad, dx, dy, dz;

            int id = part.id() - 1;

            // loop over other particles
            for(int j = 0; j < totalParticles; j++) {

                // assume triply periodic, check the domain and the 8 periodic images
                for(int ii = iilo; ii <= iihi; ii++) {
                for(int jj = jjlo; jj <= jjhi; jj++) {
                for(int kk = kklo; kk <= kkhi; kk++) {

                    // get distance between particles
                    dx = part.pos(0)-posx[j] - ii*domx;
                    dy = part.pos(1)-posy[j] - jj*domy;
                    dz = part.pos(2)-posz[j] - kk*domz;

                    int jSpec = species[j]-1;
                    rad = sqrt(dx*dx + dy*dy + dz*dz);

                    if((nearest[id*(nspecies+1) + jSpec] == 0 || nearest[id*(nspecies+1) + jSpec] > rad) && rad != 0)
                    {
                        nearest[id*(nspecies +1)+ jSpec] = rad;
                    }

                    if((nearest[id*(nspecies+1) + nspecies] == 0 || nearest[id*(nspecies+1) + nspecies] > rad) && rad != 0)
                    {
                        nearest[id*(nspecies +1) + nspecies] = rad;
//                        Print() << "Particle " << i << " species " << jSpec << ", " << nearest[i*nspecies + jSpec] << "\n";
                    }

                    // if particles are close enough, increment the bin
                    if(rad < totalDist && rad > 0.) {

                        bin = (int)amrex::Math::floor(rad/binSize);
                        radDist[bin]++;

                        if (part.rdata(FHD_realData::q) > 0) {
                            if (charge[j] > 0) {
                                radDist_pp[bin]++;
                            }
                            else if (charge[j] < 0) {
                                radDist_pm[bin]++;
                            }
                        }
                        else if (part.rdata(FHD_realData::q) < 0) {
                            if (charge[j] > 0) {
                                radDist_pm[bin]++;
                            }
                            else if (charge[j] < 0) {
                                radDist_mm[bin]++;
                            }
                        }
                    }
                }
                }
                }
            }  // loop over j (total particles)
        } // loop over i (np; local particles)
    }

    // compute total number density
    double n0_total = 0.;
    for (int i=0; i<nspecies; ++i) {
        n0_total += particleInfo[i].n0;
    }

    // collect the hit count
    ParallelDescriptor::ReduceRealSum(radDist   .dataPtr(),totalBins);
    ParallelDescriptor::ReduceRealSum(radDist_pp.dataPtr(),totalBins);
    ParallelDescriptor::ReduceRealSum(radDist_pm.dataPtr(),totalBins);
    ParallelDescriptor::ReduceRealSum(radDist_mm.dataPtr(),totalBins);
    ParallelDescriptor::ReduceRealSum(nearest,(nspecies+1)*totalParticles);

    int specCount[nspecies];

    for(int i=0;i<nspecies;i++)
    {
        specCount[i] = 0;
    }

    for(int i=0;i<totalParticles;i++)
    {
        int iSpec = species[i]-1;

        specCount[iSpec]++;

        for(int j=0;j<nspecies;j++)
        {
            nn[(iSpec)*nspecies + j] = nn[(iSpec)*nspecies + j] + nearest[i*(nspecies+1)+j];

            //std::cout << (iSpec)*nspecies + j  << ", " << nearest[i*(nspecies+1)+j] << "\n";

        }

        nn[nspecies*nspecies] = nn[nspecies*nspecies] + nearest[i*(nspecies+1)+nspecies];
    }

    for(int i=0;i<nspecies;i++)
    {
        for(int j=0;j<nspecies;j++)
        {
            if(i==j)
            {
                nn[i*nspecies+j] = nn[i*nspecies+j]/(specCount[i]-1);
            }else
            {
                nn[i*nspecies+j] = nn[i*nspecies+j]/specCount[i];
            }

            //Print() << nn[i*nspecies+j] << "\n";

        }
    }

    nn[nspecies*nspecies] = nn[nspecies*nspecies]/(totalParticles - 1);

    //Print() << nn[nspecies*nspecies+1] << "\n";


    // normalize by 1 / (number density * bin volume * total particle count)
    for(int i=0;i<totalBins;i++) {
        radDist   [i] *= 1./(n0_total*binVolRadial[i]*(double)totalParticles);
        radDist_pp[i] *= 1./(n0_total*binVolRadial[i]*(double)totalParticles);
        radDist_pm[i] *= 1./(n0_total*binVolRadial[i]*(double)totalParticles);
        radDist_mm[i] *= 1./(n0_total*binVolRadial[i]*(double)totalParticles);
    }

    // increment number of snapshots
    radialStatsCount++;
    int stepsminusone = radialStatsCount - 1;
    double stepsinv = 1.0/(double)radialStatsCount;

    // update the mean radial distribution
    for(int i=0;i<totalBins;i++) {
        meanRadialDistribution   [i] = (meanRadialDistribution   [i]*stepsminusone + radDist   [i])*stepsinv;
        meanRadialDistribution_pp[i] = (meanRadialDistribution_pp[i]*stepsminusone + radDist_pp[i])*stepsinv;
        meanRadialDistribution_pm[i] = (meanRadialDistribution_pm[i]*stepsminusone + radDist_pm[i])*stepsinv;
        meanRadialDistribution_mm[i] = (meanRadialDistribution_mm[i]*stepsminusone + radDist_mm[i])*stepsinv;
    }

    for(int i=0;i<(nspecies*nspecies + 1);i++) {

        //Print() << i << ", " << nn[i]<< "\n";
        nearestN[i] = (nearestN[i]*stepsminusone + nn[i])*stepsinv;
    }

    // output mean radial distribution g(r) based on plot_int
    if (plot_int > 0 && step%plot_int == 0) {

        if(ParallelDescriptor::MyProc() == 0) {

            std::string filename = Concatenate("radialDistribution",step,9);;
            std::ofstream ofs(filename, std::ofstream::out);

            // normalize by
            for(int i=0;i<totalBins;i++) {
                ofs << (i+0.5)*binSize << " "
                    << meanRadialDistribution   [i] << " "
                    << meanRadialDistribution_pp[i] << " "
                    << meanRadialDistribution_pm[i] << " "
                    << meanRadialDistribution_mm[i] << std::endl;
            }
            ofs.close();
        }
    }

    if (plot_int > 0 && step%plot_int == 0) {

        if(ParallelDescriptor::MyProc() == 0) {

            std::string filename = "nearestNeighbour";
            std::ofstream ofs(filename, std::ofstream::app);

            for(int i=0;i<nspecies*nspecies+1;i++) {
                ofs << nearestN[i] << " ";
            }
            ofs << std::endl;
            ofs.close();
        }
    }

    // reset the radial distribution at n_steps_skip (if n_steps_skip > 0)
    // OR
    // reset the radial distribution every |n_steps_skip| (if n_steps_skip < 0)
    if ((n_steps_skip > 0 && step == n_steps_skip) ||
        (n_steps_skip < 0 && step%n_steps_skip == 0) ) {

        Print() << "Resetting radial distribution collection.\n";

        radialStatsCount = 0;

        for(int i=0;i<totalBins;i++) {
            meanRadialDistribution   [i] = 0;
            meanRadialDistribution_pp[i] = 0;
            meanRadialDistribution_pm[i] = 0;
            meanRadialDistribution_mm[i] = 0;
        }

        for(int i=0;i<nspecies*nspecies;i++) {
            nearestN[i] = 0;
        }
    }
}

void FhdParticleContainer::potentialDistribution(long totalParticles, const int step, const species* particleInfo)
{
    const int lev = 0;
    int bin;
    Real totalDist;

    Print() << "Calculating potential distribution\n";

    // outer radial extent
    totalDist = totalBins*binSize;

    // this is the bin "hit count"
    RealVector radDist   (totalBins, 0.);

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (FhdParIter pti(*this, lev); pti.isValid(); ++pti) {

        const int grid_id = pti.index();
        const int tile_id = pti.LocalTileIndex();

        auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];
        auto& particles = particle_tile.GetArrayOfStructs();
        const int np = particles.numParticles();

        // loop over particles
        for (int i = 0; i < np; ++i) {

            ParticleType & part = particles[i];

            Real potential = part.rdata(FHD_realData::potential);

            int bin = (int)floor(potential/binSize);
            if(bin < totalBins)
            {
                radDist[bin]++;
            }
        } // loop over i (np; local particles)
    }

    // collect the hit count
    ParallelDescriptor::ReduceRealSum(radDist   .dataPtr(),totalBins);

    // increment number of snapshots
    radialStatsCount++;
    int stepsminusone = radialStatsCount - 1;
    double stepsinv = 1.0/(double)radialStatsCount;

    // update the mean radial distribution
    for(int i=0;i<totalBins;i++) {
        meanRadialDistribution   [i] = (meanRadialDistribution   [i]*stepsminusone + (radDist[i])/(totalParticles*binSize))*stepsinv;
    }

    // output mean radial distribution g(r) based on plot_int
    if (plot_int > 0 && step%plot_int == 0) {

        if(ParallelDescriptor::MyProc() == 0) {

            std::string filename = Concatenate("potentialDistribution",step,9);;
            std::ofstream ofs(filename, std::ofstream::out);

            // normalize by
            for(int i=0;i<totalBins;i++) {
                ofs << (i+0.5)*binSize << " "
                    << meanRadialDistribution[i] << std::endl;
            }
            ofs.close();
        }
    }

    // reset the radial distribution at n_steps_skip (if n_steps_skip > 0)
    // OR
    // reset the radial distribution every |n_steps_skip| (if n_steps_skip < 0)
    if ((n_steps_skip > 0 && step == n_steps_skip) ||
        (n_steps_skip < 0 && step%n_steps_skip == 0) ) {

        Print() << "Resetting potential distribution collection.\n";

        radialStatsCount = 0;

        for(int i=0;i<totalBins;i++) {
            meanRadialDistribution   [i] = 0;
        }
    }
}

void FhdParticleContainer::CartesianDistribution(long totalParticles, const int step, const species* particleInfo)
{
    BL_PROFILE_VAR("CartesianDistribution()",CartesianDistribution);

    const int lev = 0;
    int bin;
    double domx, domy, domz, totalDist, temp;

    domx = (prob_hi[0] - prob_lo[0]);
    domy = (prob_hi[1] - prob_lo[1]);
    domz = (prob_hi[2] - prob_lo[2]);

    Real posx[totalParticles];
    Real posy[totalParticles];
    Real posz[totalParticles];

    Real charge[totalParticles];

    Print() << "Calculating Cartesian distribution\n";

    // collect particle positions onto one processor
    PullDown(0, posx, -1, totalParticles);
    PullDown(0, posy, -2, totalParticles);
    PullDown(0, posz, -3, totalParticles);
    PullDown(0, charge, 27, totalParticles);

    // outer extent
    totalDist = totalBins*binSize;

    // this is the bin "hit count"
    RealVector XDist   (totalBins, 0.);
    RealVector XDist_pp(totalBins, 0.);
    RealVector XDist_pm(totalBins, 0.);
    RealVector XDist_mm(totalBins, 0.);
    RealVector YDist   (totalBins, 0.);
    RealVector YDist_pp(totalBins, 0.);
    RealVector YDist_pm(totalBins, 0.);
    RealVector YDist_mm(totalBins, 0.);
    RealVector ZDist   (totalBins, 0.);
    RealVector ZDist_pp(totalBins, 0.);
    RealVector ZDist_pm(totalBins, 0.);
    RealVector ZDist_mm(totalBins, 0.);

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (FhdParIter pti(*this, lev); pti.isValid(); ++pti) {

        const int grid_id = pti.index();
        const int tile_id = pti.LocalTileIndex();

        auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];
        auto& particles = particle_tile.GetArrayOfStructs();
        const int np = particles.numParticles();

        // loop over particles
        for (int i = 0; i < np; ++i) {

            ParticleType & part = particles[i];

            int iilo = (part.pos(0)-searchDist <= prob_lo[0]) ? -1 : 0;
            int iihi = (part.pos(0)+searchDist >= prob_hi[0]) ?  1 : 0;

            int jjlo = (part.pos(1)-searchDist <= prob_lo[1]) ? -1 : 0;
            int jjhi = (part.pos(1)+searchDist >= prob_hi[1]) ?  1 : 0;

            int kklo = (part.pos(2)-searchDist <= prob_lo[2]) ? -1 : 0;
            int kkhi = (part.pos(2)+searchDist >= prob_hi[2]) ?  1 : 0;

            double dist, dx, dy, dz;
            // loop over other particles
            for(int j = 0; j < totalParticles; j++)
            {
                // assume triply periodic, check the domain and the 8 periodic images
                for(int ii = iilo; ii <= iihi; ii++)
                {
                for(int jj = jjlo; jj <= jjhi; jj++)
                {
                for(int kk = kklo; kk <= kkhi; kk++)
                {
                    // get distance between particles
                    dx = amrex::Math::abs(part.pos(0)-posx[j] - ii*domx);
                    dy = amrex::Math::abs(part.pos(1)-posy[j] - jj*domy);
                    dz = amrex::Math::abs(part.pos(2)-posz[j] - kk*domz);

                    dist = sqrt(dx*dx + dy*dy + dz*dz);

                    // if particles are close enough, increment the bin
                    if (dist > 0.) {
                        if(dx < totalDist && dy < searchDist && dz < searchDist) {

                            bin = (int)amrex::Math::floor(dx/binSize);
                            XDist[bin]++;

                            if (part.rdata(FHD_realData::q) > 0) {
                                if (charge[j] > 0) {
                                    XDist_pp[bin]++;
                                }
                                else if (charge[j] < 0) {
                                    XDist_pm[bin]++;
                                }
                            }
                            else if (part.rdata(FHD_realData::q) < 0) {
                                if (charge[j] > 0) {
                                    XDist_pm[bin]++;
                                }
                                else if (charge[j] < 0) {
                                    XDist_mm[bin]++;
                                }
                            }
                        }
                        if(dy < totalDist && dx < searchDist && dz < searchDist) {

                            bin = (int)amrex::Math::floor(dy/binSize);
                            YDist[bin]++;

                            if (part.rdata(FHD_realData::q) > 0) {
                                if (charge[j] > 0) {
                                    YDist_pp[bin]++;
                                }
                                else if (charge[j] < 0) {
                                    YDist_pm[bin]++;
                                }
                            }
                            else if (part.rdata(FHD_realData::q) < 0) {
                                if (charge[j] > 0) {
                                    YDist_pm[bin]++;
                                }
                                else if (charge[j] < 0) {
                                    YDist_mm[bin]++;
                                }
                            }
                        }
                        if(dz < totalDist && dx < searchDist && dy < searchDist) {

                            bin = (int)amrex::Math::floor(dz/binSize);
                            ZDist[bin]++;

                            if (part.rdata(FHD_realData::q) > 0) {
                                if (charge[j] > 0) {
                                    ZDist_pp[bin]++;
                                }
                                else if (charge[j] < 0) {
                                    ZDist_pm[bin]++;
                                }
                            }
                            else if (part.rdata(FHD_realData::q) < 0) {
                                if (charge[j] > 0) {
                                    ZDist_pm[bin]++;
                                }
                                else if (charge[j] < 0) {
                                    ZDist_mm[bin]++;
                                }
                            }
                        }
                    }
                }
                }
                }
            }
        }
    }

    // compute total number density
    double n0_total = 0.;
    for (int i=0; i<nspecies; ++i) {
        n0_total += particleInfo[i].n0;
    }

    // collect the hit count
    ParallelDescriptor::ReduceRealSum(XDist   .dataPtr(),totalBins);
    ParallelDescriptor::ReduceRealSum(XDist_pp.dataPtr(),totalBins);
    ParallelDescriptor::ReduceRealSum(XDist_pm.dataPtr(),totalBins);
    ParallelDescriptor::ReduceRealSum(XDist_mm.dataPtr(),totalBins);
    ParallelDescriptor::ReduceRealSum(YDist   .dataPtr(),totalBins);
    ParallelDescriptor::ReduceRealSum(YDist_pp.dataPtr(),totalBins);
    ParallelDescriptor::ReduceRealSum(YDist_pm.dataPtr(),totalBins);
    ParallelDescriptor::ReduceRealSum(YDist_mm.dataPtr(),totalBins);
    ParallelDescriptor::ReduceRealSum(ZDist   .dataPtr(),totalBins);
    ParallelDescriptor::ReduceRealSum(ZDist_pp.dataPtr(),totalBins);
    ParallelDescriptor::ReduceRealSum(ZDist_pm.dataPtr(),totalBins);
    ParallelDescriptor::ReduceRealSum(ZDist_mm.dataPtr(),totalBins);

    // normalize by 1 / (number density * bin volume * total particle count)
    for(int i=0;i<totalBins;i++) {
        XDist   [i] *= 1./(n0_total*binVolCartesian*(double)totalParticles);
        XDist_pp[i] *= 1./(n0_total*binVolCartesian*(double)totalParticles);
        XDist_pm[i] *= 1./(n0_total*binVolCartesian*(double)totalParticles);
        XDist_mm[i] *= 1./(n0_total*binVolCartesian*(double)totalParticles);
        YDist   [i] *= 1./(n0_total*binVolCartesian*(double)totalParticles);
        YDist_pp[i] *= 1./(n0_total*binVolCartesian*(double)totalParticles);
        YDist_pm[i] *= 1./(n0_total*binVolCartesian*(double)totalParticles);
        YDist_mm[i] *= 1./(n0_total*binVolCartesian*(double)totalParticles);
        ZDist   [i] *= 1./(n0_total*binVolCartesian*(double)totalParticles);
        ZDist_pp[i] *= 1./(n0_total*binVolCartesian*(double)totalParticles);
        ZDist_pm[i] *= 1./(n0_total*binVolCartesian*(double)totalParticles);
        ZDist_mm[i] *= 1./(n0_total*binVolCartesian*(double)totalParticles);
    }

    // increment number of snapshots
    cartesianStatsCount++;
    int stepsminusone = cartesianStatsCount - 1;
    double stepsinv = 1.0/(double)cartesianStatsCount;

    // update the mean Cartesian distribution
    for(int i=0;i<totalBins;i++) {
        meanXDistribution   [i] = (meanXDistribution   [i]*stepsminusone + XDist   [i])*stepsinv;
        meanXDistribution_pp[i] = (meanXDistribution_pp[i]*stepsminusone + XDist_pp[i])*stepsinv;
        meanXDistribution_pm[i] = (meanXDistribution_pm[i]*stepsminusone + XDist_pm[i])*stepsinv;
        meanXDistribution_mm[i] = (meanXDistribution_mm[i]*stepsminusone + XDist_mm[i])*stepsinv;
        meanYDistribution   [i] = (meanYDistribution   [i]*stepsminusone + YDist   [i])*stepsinv;
        meanYDistribution_pp[i] = (meanYDistribution_pp[i]*stepsminusone + YDist_pp[i])*stepsinv;
        meanYDistribution_pm[i] = (meanYDistribution_pm[i]*stepsminusone + YDist_pm[i])*stepsinv;
        meanYDistribution_mm[i] = (meanYDistribution_mm[i]*stepsminusone + YDist_mm[i])*stepsinv;
        meanZDistribution   [i] = (meanZDistribution   [i]*stepsminusone + ZDist   [i])*stepsinv;
        meanZDistribution_pp[i] = (meanZDistribution_pp[i]*stepsminusone + ZDist_pp[i])*stepsinv;
        meanZDistribution_pm[i] = (meanZDistribution_pm[i]*stepsminusone + ZDist_pm[i])*stepsinv;
        meanZDistribution_mm[i] = (meanZDistribution_mm[i]*stepsminusone + ZDist_mm[i])*stepsinv;
    }

    // output mean Cartesian distribution g(x), g(y), g(z) based on plot_int
    if (plot_int > 0 && step%plot_int == 0) {

        if(ParallelDescriptor::MyProc() == 0) {

            std::string filename = Concatenate("cartesianDistribution",step,9);;
            std::ofstream ofs(filename, std::ofstream::out);

            for(int i=0;i<totalBins;i++) {
                ofs << (i+0.5)*binSize << " "
                    << meanXDistribution   [i] << " "
                    << meanXDistribution_pp[i] << " "
                    << meanXDistribution_pm[i] << " "
                    << meanXDistribution_mm[i] << " "
                    << meanYDistribution   [i] << " "
                    << meanYDistribution_pp[i] << " "
                    << meanYDistribution_pm[i] << " "
                    << meanYDistribution_mm[i] << " "
                    << meanZDistribution   [i] << " "
                    << meanZDistribution_pp[i] << " "
                    << meanZDistribution_pm[i] << " "
                    << meanZDistribution_mm[i] << " " << std::endl;
            }
            ofs.close();
        }
    }

    // reset the Cartesian distribution at n_steps_skip (if n_steps_skip > 0)
    // OR
    // reset the Cartesian distribution every |n_steps_skip| (if n_steps_skip < 0)
    if ((n_steps_skip > 0 && step == n_steps_skip) ||
        (n_steps_skip < 0 && step%n_steps_skip == 0) ) {

        Print() << "Resetting Cartesian distribution collection.\n";

        cartesianStatsCount = 0;

        for(int i=0;i<totalBins;i++) {
            meanXDistribution   [i] = 0;
            meanXDistribution_pp[i] = 0;
            meanXDistribution_pm[i] = 0;
            meanXDistribution_mm[i] = 0;
            meanYDistribution   [i] = 0;
            meanYDistribution_pp[i] = 0;
            meanYDistribution_pm[i] = 0;
            meanYDistribution_mm[i] = 0;
            meanZDistribution   [i] = 0;
            meanZDistribution_pp[i] = 0;
            meanZDistribution_pm[i] = 0;
            meanZDistribution_mm[i] = 0;
        }
    }
}

void FhdParticleContainer::collectFieldsGPU(const Real dt, const Real* dxPotential,
                                         const MultiFab& RealCenterCoords, const Geometry geomP, MultiFab& charge, MultiFab& chargeTemp,
                                         MultiFab& mass, MultiFab& massTemp)
{
    BL_PROFILE_VAR("collectFields()",collectFields);

    const int lev = 0;
    const Real* dx = Geom(lev).CellSize();
    const Real* plo = Geom(lev).ProbLo();
    const Real* phi = Geom(lev).ProbHi();

#ifdef _OPENMP
#pragma omp parallel
#endif

    charge.setVal(0.0);
    chargeTemp.setVal(0.0);

    //mass.setVal(0.0);
    //massTemp.setVal(0.0);

    for (FhdParIter pti(*this, lev); pti.isValid(); ++pti)
    {
        const int grid_id = pti.index();
        const int tile_id = pti.LocalTileIndex();
        const Box& tile_box  = pti.tilebox();

        auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];
        auto& particles = particle_tile.GetArrayOfStructs();
        const int np = particles.numParticles();

        collect_charge_gpu(particles, chargeTemp[pti], AMREX_ZFILL(geomP.ProbLo()), AMREX_ZFILL(dxPotential));
    }

    MultiFabPhysBCCharge(chargeTemp, geomP);

    chargeTemp.SumBoundary(geomP.periodicity());
    //massTemp.SumBoundary(geomP.periodicity());

    MultiFab::Add(charge,chargeTemp,0,0,charge.nComp(),charge.nGrow());
    //MultiFab::Add(mass,massTemp,0,0,mass.nComp(),mass.nGrow());

    charge.FillBoundary(geomP.periodicity());
}

void FhdParticleContainer::EvaluateStats(MultiFab& particleInstant,
                                         MultiFab& particleMeans, species particleInfo, const Real delt, int steps)
{
    BL_PROFILE_VAR("EvaluateStats()",EvaluateStats);

    const int lev = 0;
    const double Neff = particleInfo.Neff;
    const double n0 = particleInfo.n0;
    const double T0 = particleInfo.T;

    double tp = 0;
    double te = 0;
    double tm = 0;

    double totalMass;

    Gpu::DeviceVector<Real> avcurrent_tile (3);

//    GpuArray<int, 3> avcurrent_tile = {0,0,0};
    RealVector avcurrent_proc(3);

    std::fill(avcurrent_proc.begin(), avcurrent_proc.end(), 0.);
    std::fill(avcurrent_tile.begin(), avcurrent_tile.end(), 0.);

    Real* avc_tile = avcurrent_tile.dataPtr();

    BoxArray ba = particleMeans.boxArray();
    long cellcount = ba.numPts();

    const Real* dx = Geom(lev).CellSize();
    const Real dxInv = 1.0/dx[0];
    const Real cellVolInv = 1.0/(dx[0]*dx[0]*dx[0]);

    const Real stepsInv = 1.0/steps;
    const int stepsMinusOne = steps-1;

    // zero instantaneous values
    particleInstant.setVal(0.);

    for (FhdParIter pti(*this, lev); pti.isValid(); ++pti)
    {

        PairIndex index(pti.index(), pti.LocalTileIndex());
        const int np = this->GetParticles(lev)[index].numRealParticles();
        auto& plev = this->GetParticles(lev);
        auto& ptile = plev[index];
        auto& aos   = ptile.GetArrayOfStructs();
        const Box& tile_box  = pti.tilebox();
        ParticleType* particles = aos().dataPtr();

        tp = tp + np;

        GpuArray<int, 3> bx_lo = {tile_box.loVect()[0], tile_box.loVect()[1], tile_box.loVect()[2]};
        GpuArray<int, 3> bx_hi = {tile_box.hiVect()[0], tile_box.hiVect()[1], tile_box.hiVect()[2]};

//        Array4<Real> part_inst = (&particleInstant[pti])->array();

        Array4<Real> part_inst = particleInstant[pti].array();
        Array4<Real> part_mean = particleMeans[pti].array();

        //const Array4<Real>& data = charge.array(mfi);
        // Updates multfab
        AMREX_FOR_1D( np, ni,
        {
            ParticleType & part = particles[ni];

            int i = floor(part.pos(0)*dxInv);
            int j = floor(part.pos(1)*dxInv);
            int k = floor(part.pos(2)*dxInv);

            amrex::Gpu::Atomic::Add(&part_inst(i,j,k,0), 1.0);
            amrex::Gpu::Atomic::Add(&part_inst(i,j,k,1), part.rdata(FHD_realData::mass));

            amrex::Gpu::Atomic::Add(&part_inst(i,j,k,2), part.rdata(FHD_realData::velx));
            amrex::Gpu::Atomic::Add(&part_inst(i,j,k,3), part.rdata(FHD_realData::vely));
            amrex::Gpu::Atomic::Add(&part_inst(i,j,k,4), part.rdata(FHD_realData::velz));

            amrex::Gpu::Atomic::Add(&part_inst(i,j,k,5), part.rdata(FHD_realData::velx)*part.rdata(FHD_realData::q));
            amrex::Gpu::Atomic::Add(&part_inst(i,j,k,6), part.rdata(FHD_realData::vely)*part.rdata(FHD_realData::q));
            amrex::Gpu::Atomic::Add(&part_inst(i,j,k,7), part.rdata(FHD_realData::velz)*part.rdata(FHD_realData::q));

            amrex::Gpu::Atomic::Add(&part_inst(i,j,k,7+part.idata(FHD_intData::species)), part.rdata(FHD_realData::q));
        });

        amrex::ParallelFor(tile_box,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real membersInv = 1.0/part_inst(i,j,k,0);

            part_inst(i,j,k,1) = part_inst(i,j,k,1)*cellVolInv;

            part_inst(i,j,k,2) = part_inst(i,j,k,2)*membersInv;
            part_inst(i,j,k,3) = part_inst(i,j,k,3)*membersInv;
            part_inst(i,j,k,4) = part_inst(i,j,k,4)*membersInv;

            part_inst(i,j,k,5) = part_inst(i,j,k,5)*cellVolInv;
            part_inst(i,j,k,6) = part_inst(i,j,k,6)*cellVolInv;
            part_inst(i,j,k,7) = part_inst(i,j,k,7)*cellVolInv;

            for(int l=0;l<nspecies;l++)
            {
                part_inst(i,j,k,8 + l) = part_inst(i,j,k,8 + l)*cellVolInv;
            }
        });

        amrex::ParallelFor(tile_box,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            part_mean(i,j,k,1)  = (part_mean(i,j,k,1)*stepsMinusOne + part_inst(i,j,k,1))*stepsInv;
            part_mean(i,j,k,2)  = (part_mean(i,j,k,2)*stepsMinusOne + part_inst(i,j,k,2))*stepsInv;
            part_mean(i,j,k,3)  = (part_mean(i,j,k,3)*stepsMinusOne + part_inst(i,j,k,3))*stepsInv;
            part_mean(i,j,k,4)  = (part_mean(i,j,k,4)*stepsMinusOne + part_inst(i,j,k,4))*stepsInv;
            part_mean(i,j,k,5)  = (part_mean(i,j,k,5)*stepsMinusOne + part_inst(i,j,k,5))*stepsInv;
            part_mean(i,j,k,6)  = (part_mean(i,j,k,6)*stepsMinusOne + part_inst(i,j,k,6))*stepsInv;
            part_mean(i,j,k,7)  = (part_mean(i,j,k,7)*stepsMinusOne + part_inst(i,j,k,7))*stepsInv;

            for(int l=0;l<nspecies;l++)
            {
                part_mean(i,j,k,8 + l)  = (part_mean(i,j,k,8 + l)*stepsMinusOne + part_inst(i,j,k,8 + l))*stepsInv;
            }

            avc_tile[0] = avc_tile[0] + part_mean(i,j,k,5);
            avc_tile[1] = avc_tile[1] + part_mean(i,j,k,6);
            avc_tile[2] = avc_tile[2] + part_mean(i,j,k,7);
        });

        for (int i=0; i<3; ++i) {
            avcurrent_proc[i] += avcurrent_tile[i];
        }
    }

    // gather statistics
    ParallelDescriptor::ReduceRealSum(avcurrent_proc.dataPtr(),3);

    // write out current mean and variance to file
    if(ParallelDescriptor::MyProc() == 0) {
        std::string filename = "currentEst";
        std::ofstream ofs(filename, std::ofstream::app);
        ofs << avcurrent_proc[0]/cellcount << "  "
            << avcurrent_proc[1]/cellcount << "  "
            << avcurrent_proc[2]/cellcount <<  "\n";
        ofs.close();
    }
}

void FhdParticleContainer::WriteParticlesAscii(std::string asciiName)
{
    BL_PROFILE_VAR("WriteParticlesAscii()",WriteParticlesAscii);

    WriteAsciiFile(asciiName);
}

void
FhdParticleContainer::PrintParticles()
{
    BL_PROFILE_VAR("PrintParticles()",PrintParticles);

    int lev =0;

    for(FhdParIter pti(* this, lev); pti.isValid(); ++pti)
    {
        PairIndex index(pti.index(), pti.LocalTileIndex());

        AoS & particles = this->GetParticles(lev).at(index).GetArrayOfStructs();
        long np = this->GetParticles(lev).at(index).numParticles();

        for(int i=0; i<np; ++i)
        {
            ParticleType & part = particles[i];

            if(part.idata(FHD_intData::pinned) == 0)
            {
                double bigM  = part.rdata(FHD_realData::totalDiff)/(T_init[0]*k_B);
                double absForce = sqrt(part.rdata(FHD_realData::forcex)*part.rdata(FHD_realData::forcex) + part.rdata(FHD_realData::forcey)*part.rdata(FHD_realData::forcey) + part.rdata(FHD_realData::forcez)*part.rdata(FHD_realData::forcez));

                std::cout << scientific << setprecision(15) << "Particle " << ParallelDescriptor::MyProc() << ", " << part.id() << ", force: " << part.rdata(FHD_realData::forcex) << ", " << part.rdata(FHD_realData::forcey) << ", " << part.rdata(FHD_realData::forcez) << ", " << absForce << std::endl;
                std::cout << scientific << setprecision(15) << "Particle " << ParallelDescriptor::MyProc() << ", " << part.id() << ", position/q: " << part.pos(0) << ", " << part.pos(1) << ", " << part.pos(2) << ", " << part.rdata(FHD_realData::q) << std::endl;
                std::cout << scientific << setprecision(15) << "Particle " << ParallelDescriptor::MyProc() << ", " << part.id() << ", vel: " << part.rdata(FHD_realData::velx) << ", " << part.rdata(FHD_realData::vely) << ", " << part.rdata(FHD_realData::velz) << std::endl;
                std::cout << scientific << setprecision(15) << "Particle " << ParallelDescriptor::MyProc() << ", " << part.id() << ", normalised mobility: " << (part.rdata(FHD_realData::velx)/part.rdata(FHD_realData::forcex))/bigM << ", " << (part.rdata(FHD_realData::vely)/part.rdata(FHD_realData::forcey))/bigM << ", " << (part.rdata(FHD_realData::velz)/part.rdata(FHD_realData::forcez))/bigM << std::endl;

                if(part.id() == 1)
                {
                    std::string filename = "threepmForce";
                    std::ofstream ofs(filename, std::ofstream::app);
                    ofs << part.rdata(FHD_realData::forcey) << ", " << absForce << std::endl;
                }
            }
        }
    }
}

void
FhdParticleContainer::TwoParticleCorrelation()
{
    BL_PROFILE_VAR("TwoParticleCorrelation()",PrintParticles);

    int lev =0;

    for(FhdParIter pti(* this, lev); pti.isValid(); ++pti)
    {
        PairIndex index(pti.index(), pti.LocalTileIndex());

        AoS & particles = this->GetParticles(lev).at(index).GetArrayOfStructs();
        long np = this->GetParticles(lev).at(index).numParticles();

        for(int i=0;i < np;i++)
        {
            ParticleType & part = particles[i];

            string fileX = Concatenate("partX", part.id());
            string fileY = Concatenate("partY", part.id());
            string fileZ = Concatenate("partZ", part.id());

            std::ofstream ofs1(fileX, std::ofstream::app);
            ofs1 << setprecision(15) << part.rdata(FHD_realData::velx) << std::endl;

            std::ofstream ofs2(fileY, std::ofstream::app);
            ofs2 << setprecision(15) << part.rdata(FHD_realData::vely) << std::endl;

            std::ofstream ofs3(fileZ, std::ofstream::app);
            ofs3 << setprecision(15) << part.rdata(FHD_realData::velz) << std::endl;

            ofs1.close();
            ofs2.close();
            ofs3.close();
        }
    }
}

void
FhdParticleContainer::SetPosition(int id, Real x, Real y, Real z)
{
    BL_PROFILE_VAR("SetPosition()",SetPosition);

    int lev =0;

    for(FhdParIter pti(* this, lev); pti.isValid(); ++pti)
    {
        PairIndex index(pti.index(), pti.LocalTileIndex());

        AoS & particles = this->GetParticles(lev).at(index).GetArrayOfStructs();
        long np = this->GetParticles(lev).at(index).numParticles();

        for(int i=0;i < np;i++)
        {
            ParticleType & part = particles[i];

            if(part.id() == id)
            {
                part.pos(0) = x;
                part.pos(1) = y;
                part.pos(2) = z;

                //std::cout << "Rank " << ParallelDescriptor::MyProc() << " particle " << id << " moving to " << x << ", " << y << ", " << z << std::endl;
            }
        }
    }

    clearNeighbors();
    Redistribute();
    fillNeighbors();
}

void
FhdParticleContainer::SetVel(int id, Real x, Real y, Real z)
{
    BL_PROFILE_VAR("SetVel()",SetVel);

    int lev =0;

    for(FhdParIter pti(* this, lev); pti.isValid(); ++pti)
    {
        PairIndex index(pti.index(), pti.LocalTileIndex());

        AoS & particles = this->GetParticles(lev).at(index).GetArrayOfStructs();
        long np = this->GetParticles(lev).at(index).numParticles();

        for(int i=0;i < np;i++)
        {
            ParticleType & part = particles[i];

            if(part.id() == id)
            {
                part.rdata(FHD_realData::velx) = x;
                part.rdata(FHD_realData::vely) = y;
                part.rdata(FHD_realData::velz) = z;
            }
            std::cout << "Rank " << ParallelDescriptor::MyProc() << " particle " << id << " moving to " << x << ", " << y << ", " << z << std::endl;
        }
    }
}

void
FhdParticleContainer::SetForce(int id, Real x, Real y, Real z)
{
    BL_PROFILE_VAR("SetForce()",SetVel);

    int lev =0;

    for(FhdParIter pti(* this, lev); pti.isValid(); ++pti)
    {
        PairIndex index(pti.index(), pti.LocalTileIndex());

        AoS & particles = this->GetParticles(lev).at(index).GetArrayOfStructs();
        long np = this->GetParticles(lev).at(index).numParticles();

        for(int i=0;i < np;i++)
        {
            ParticleType & part = particles[i];

            if(part.id() == id)
            {
                part.rdata(FHD_realData::forcex) = x;
                part.rdata(FHD_realData::forcey) = y;
                part.rdata(FHD_realData::forcez) = z;
                std::cout << "Rank " << ParallelDescriptor::MyProc() << " particle " << id << " forcing with " << x << ", " << y << ", " << z << std::endl;
            }else
            {
                part.rdata(FHD_realData::forcex) = 0;
                part.rdata(FHD_realData::forcey) = 0;
                part.rdata(FHD_realData::forcez) = 0;
            }
        }
    }
}

void
FhdParticleContainer::clearMobilityMatrix()
{
    for(int i=0;i<(3*totalPinnedMarkers);i++)
    {
        for(int j=0;j<(3*totalPinnedMarkers);j++)
        {
            pinMatrix[i*3*totalPinnedMarkers + j] = 0;
            //Print() << i << ", " << j << ", " << i*3*totalPinnedMarkers + j << ", " << pinMatrix[i*totalPinnedMarkers + j] << std::endl;
        }
    }
}

void
FhdParticleContainer::fillMobilityMatrix(int id, int comp)
{
    int lev = 0;

    Real velx[totalMarkers];
    Real vely[totalMarkers];
    Real velz[totalMarkers];
    int  pinned[totalMarkers];
    Real  forcex[totalMarkers];
    Real  forcey[totalMarkers];
    Real  forcez[totalMarkers];

    int idMap[totalMarkers];

    PullDown(0, velx, FHD_realData::velx, totalMarkers);
    PullDown(0, vely, FHD_realData::vely, totalMarkers);
    PullDown(0, velz, FHD_realData::velz, totalMarkers);
    PullDownInt(0, pinned, FHD_intData::pinned, totalMarkers);

    for(FhdParIter pti(* this, lev); pti.isValid(); ++pti)
    {
        PairIndex index(pti.index(), pti.LocalTileIndex());

        AoS & particles = this->GetParticles(lev).at(index).GetArrayOfStructs();
        long np = this->GetParticles(lev).at(index).numParticles();

//        AMREX_FOR_1D( np, i,
        for(int i=0;i<np;i++)
        {
            ParticleType & part = particles[i];

            velx[part.id()-1] = part.rdata(FHD_realData::velx);
            vely[part.id()-1] = part.rdata(FHD_realData::vely);
            velz[part.id()-1] = part.rdata(FHD_realData::velz);
        }
    }

    Real velpin[totalPinnedMarkers*3];
    int k = 0;

    for(int i=0;i<totalMarkers;i++)
    {
        if(pinned[i] == 1)
        {
            velpin[k] = velx[i];
            velpin[k+1] = vely[i];
            velpin[k+2] = velz[i];

            k = k+3;
        }
    }

    k=1;
    for(int i=0;i<totalMarkers;i++)
    {
        if(pinned[i] == 1)
        {
            idMap[i] = k;
            k++;
        }
    }

    int realID = idMap[id-1];

    int matrixSize = 3*totalPinnedMarkers;

    for(int i=0;i<matrixSize;i=i+3)
    {
        int j = 3*(realID-1) + comp;

        pinMatrix[i*matrixSize +j] += velpin[i];
        pinMatrix[(i+1)*matrixSize +j] += velpin[i+1];
        pinMatrix[(i+2)*matrixSize +j] += velpin[i+2];

//        Print() << "MAT: " << i*matrixSize +j << ", " << pinMatrix[i*matrixSize+j] << std::endl;
//        Print() << "MAT: " << (i+1)*matrixSize +j << ", " << pinMatrix[(i+1)*matrixSize+j] << std::endl;
//        Print() << "MAT: " << (i+2)*matrixSize +j << ", " << pinMatrix[(i+2)*matrixSize+j] << std::endl;
    }

//    Print() << "Real ID: " << realID << std::endl;
//    Print() << "MAT: " << pinMatrix[0] << std::endl;
//    Print() << "MAT: " << pinMatrix[2*matrixSize] << std::endl;


}

void
FhdParticleContainer::writeVel(int id)
{
    int lev =0;

    Real x = 0;
    Real y = 0;
    Real z = 0;

    for(FhdParIter pti(* this, lev); pti.isValid(); ++pti)
    {
        PairIndex index(pti.index(), pti.LocalTileIndex());

        AoS & particles = this->GetParticles(lev).at(index).GetArrayOfStructs();
        long np = this->GetParticles(lev).at(index).numParticles();

        for(int i=0;i < np;i++)
        {
            ParticleType & part = particles[i];

            if(part.id() == id)
            {
                x = part.rdata(FHD_realData::velx);
                y = part.rdata(FHD_realData::vely);
                z = part.rdata(FHD_realData::velz);
            }
        }
    }

    ParallelDescriptor::ReduceRealSum(x);
    ParallelDescriptor::ReduceRealSum(y);
    ParallelDescriptor::ReduceRealSum(z);

    if(ParallelDescriptor::MyProc() == 0) {
        std::string filename = "velOut";
        std::ofstream ofs(filename, std::ofstream::app);

        ofs << x << std::endl;
        ofs << y << std::endl;
        ofs << z << std::endl;

        ofs.close();
    }
}

void
FhdParticleContainer::invertMatrix()
{
    // timer for profiling
    BL_PROFILE_VAR("invertMatrix()",initRankLists);

    if(ParallelDescriptor::MyProc() == 0) {
        int N = 3*totalPinnedMarkers;

        Real* inv = new Real[N*N];
        Real* A = new Real[N*N];

        Real time1 = ParallelDescriptor::second();

        for(int i=0;i<N;i++)
        {
            for(int j=0;j<N;j++)
            {
                A[i*N + j] = pinMatrix[i*N + j];
            }
        }

        for(int i = 0; i < N; i++){
            for(int j = 0; j < N; j++){
                if(i==j)
                    inv[i*N + j] = 1.0;
                else
                    inv[i*N + j] = 0.0;
            }
        }

        for(int i = 0; i < N; i++){
            for(int j = 0; j < N; j++){
                if(i!=j){
                    Real ratio = A[j*N + i]/A[i*N + i];
                    for(int k = 0; k < N; k++){
                        A[j*N + k] -= ratio * A[i*N + k];
                        inv[j*N + k] -= ratio * inv[i*N + k];
                    }
                }
            }
        }

        for(int i = 0; i < N; i++){
            Real a = A[i*N + i];
            for(int j = 0; j < N; j++){
                A[i*N + j] /= a;
                inv[i*N + j] /= a;
            }
        }


//        std::string filename1 = "matOut";
//        std::ofstream ofs1(filename1, std::ofstream::app);

//        int matrixSize = 3*totalPinnedMarkers;
//
//        for(int i=0;i<matrixSize;i++)
//        {
//            for(int j=0;j<matrixSize;j++)
//            {
//                ofs1 << setprecision(15) << pinMatrix[i*matrixSize +j] << std::endl;

//            }
//        }
//
//        ofs1.close();

        Real time2 = ParallelDescriptor::second() - time1;

        Print() << "Inverse matrix calculated in " << time2 << " seconds.\n";

//        std::string filename1 = "permOut";
//        remove("permOut");
//        std::ofstream ofs1(filename1, std::ofstream::app);

//
//        for(int i=0;i<N;i++)
//        {
//            for(int j=0;j<N;j++)
//            {
//                ofs1 << setprecision(15) << A[i*N + j] << std::endl;

//            }
//        }
//
//        ofs1.close();


        std::string filename = "invOut";
        remove("invOut");
        //std::ofstream ofs2(filename2, std::ofstream::app);
        ofstream ofs( filename, ios::binary );

        for(int i=0;i<N;i++)
        {
            for(int j=0;j<N;j++)
            {
                //ofs2 << setprecision(15) << inv[i*N + j] << std::endl;
                Real element = inv[i*N + j];
                ofs.write( reinterpret_cast<char*>( &element ), sizeof element );
            }
        }

//        ofs2.close();
        ofs.close();
        Print() << "Inverse matrix written\n";

        delete inv;
        delete A;
    }
}

void
FhdParticleContainer::writeMat()
{
    if(ParallelDescriptor::MyProc() == 0) {
        std::string filename = "matOut";

        remove("matOut");
        //std::ofstream ofs(filename, std::ofstream::app);
        ofstream ofs( filename, ios::binary );

        int matrixSize = 3*totalPinnedMarkers;

        for(int i=0;i<matrixSize;i++)
        {
            for(int j=0;j<matrixSize;j++)
            {
                //ofs << setprecision(15) << pinMatrix[i*matrixSize +j] << std::endl;

                Real element = pinMatrix[i*matrixSize +j];
                ofs.write( reinterpret_cast<char*>( &element ), sizeof element );
            }
        }

        ofs.close();
        Print() << "WRITTEN\n";
    }
}

void
FhdParticleContainer::MeanSqrCalc(int lev, int step) {

    BL_PROFILE_VAR("MeanSqrCalc()",MeanSqrCalc);

    Real diffTotal = 0;
    Real tt = 0.; // set to zero to protect against grids with no particles
    long nTotal = 0;
    Real sumPosQ[3] = {0,0,0};
    Real sqrDispX[nspecies];
    Real sqrDispY[nspecies];
    Real sqrDispZ[nspecies];
    Real sqrDisp[nspecies];
    Real travelTime[nspecies];
    int specCount[nspecies];

    int stepstat[nspecies];

    for(int i=0;i<nspecies;i++)
    {
        stepstat[i] = fmod(step-1,msd_int[i]);
        cout << "remainder: " << stepstat[i] << endl;
    }

    for(int i=0;i<nspecies;i++)
    {
        sqrDispX[i]=0;
        sqrDispY[i]=0;
        sqrDispZ[i]=0;
        sqrDisp[i]=0;
        specCount[i]=0;
    }

    for (MyIBMarIter pti(* this, lev); pti.isValid(); ++pti) {

        TileIndex index(pti.index(), pti.LocalTileIndex());

        AoS & particles = this->GetParticles(lev).at(index).GetArrayOfStructs();
        long np = this->GetParticles(lev).at(index).numParticles();
        nTotal += np;

        for (int i=0; i<np; ++i) {
            ParticleType & part = particles[i];
            int spec = part.idata(FHD_intData::species)-1;

            if(stepstat[spec] == 0)
            {
                for (int d=0; d<AMREX_SPACEDIM; ++d){
                    part.rdata(FHD_realData::ox + d) = part.rdata(FHD_realData::ax + d);
                }
                part.rdata(FHD_realData::travelTime) = 0;
            }
        }

        for (int i=0; i<np; ++i) {
            ParticleType & part = particles[i];

            int spec = part.idata(FHD_intData::species)-1;

            Real dispX = pow(part.rdata(FHD_realData::ax)-part.rdata(FHD_realData::ox),2);
            Real dispY = pow(part.rdata(FHD_realData::ay)-part.rdata(FHD_realData::oy),2);
            Real dispZ = pow(part.rdata(FHD_realData::az)-part.rdata(FHD_realData::oz),2);

            sqrDispX[spec] += dispX;
            sqrDispY[spec] += dispY;
            sqrDispZ[spec] += dispZ;
            sqrDisp[spec] += dispX + dispY + dispZ;

            travelTime[spec] = part.rdata(FHD_realData::travelTime);
            specCount[spec]++;
        }
    }

    for(int i=0;i<nspecies;i++)
    {
        Real temp = sqrDispX[i];
        ParallelDescriptor::ReduceRealSum(temp);
        sqrDispX[i] = temp;

        temp = sqrDispY[i];
        ParallelDescriptor::ReduceRealSum(temp);
        sqrDispY[i] = temp;

        temp = sqrDispZ[i];
        ParallelDescriptor::ReduceRealSum(temp);
        sqrDispZ[i] = temp;

        temp = sqrDisp[i];
        ParallelDescriptor::ReduceRealSum(temp);
        sqrDisp[i] = temp;

        temp = travelTime[i];
        ParallelDescriptor::ReduceRealMax(temp);
        travelTime[i] = temp;

        int itemp = specCount[i];
        ParallelDescriptor::ReduceIntSum(itemp);
        specCount[i] = itemp;
    }

    for(int i=0;i<nspecies;i++)
    {
        if(specCount[i] != 0)
        {
            sqrDispX[i] /= specCount[i];
            sqrDispY[i] /= specCount[i];
            sqrDispZ[i] /= specCount[i];
            sqrDisp[i] /= specCount[i];
        }
    }

    if(ParallelDescriptor::MyProc() == 0) {
        for(int i=0;i<nspecies;i++)
        {
            if(msd_int[i] > 0)
            {
                std::string specname = Concatenate("msdEst_",i+1);
                std::ofstream ofs(specname, std::ofstream::app);

                if(stepstat[i]==0)
                {
                    ofs << std::endl;
                }else if(stepstat[i]<(msd_len[i]+1))
                {
                    ofs << travelTime[i] << "  " << sqrDispX[i] << "  " << sqrDispY[i] << "  "<< sqrDispZ[i] << "  "<< sqrDisp[i] << std::endl;
                }

                ofs.close();
            }
        }
    }
}

void
FhdParticleContainer::GetAllParticlePositions(Real* posx, Real* posy, Real* posz, int totalParticles) {
    // collect particle positions onto one processor
    PullDown(0, posx, -1, totalParticles);
    PullDown(0, posy, -2, totalParticles);
    PullDown(0, posz, -3, totalParticles);
}

void
FhdParticleContainer::BuildCorrectionTable(const Real* dx, int setMeasureFinal) {

    BL_PROFILE_VAR("BuildCorrectionTable()",BuildCorrectionTable);

    int lev = 0;

    Real x0,y0,z0,x1,y1,z1, costheta, sintheta, cosphi, sinphi;

    Real dr = threepmCurrentBin*(threepmRange/threepmBins)*dx[0];

    if(threepmCurrentBin == threepmBins)
    {
        setMeasureFinal = 2;
    }
    if(threepmCurrentBin > threepmBins)
    {
        setMeasureFinal = 3;
    }

    for (MyIBMarIter pti(* this, lev); pti.isValid(); ++pti) {

        TileIndex index(pti.index(), pti.LocalTileIndex());

        AoS & particles = this->GetParticles(lev).at(index).GetArrayOfStructs();
        long np = this->GetParticles(lev).at(index).numParticles();

        if(setMeasureFinal == 0)
        {
            ParticleType & part0 = particles[0];
            ParticleType & part1 = particles[1];

            x0 = prob_lo[0] + 0.25*(prob_hi[0]-prob_lo[0]) + amrex::Random()*(prob_hi[0]-prob_lo[0])*0.5;
            y0 = prob_lo[1] + 0.25*(prob_hi[0]-prob_lo[0]) + amrex::Random()*(prob_hi[1]-prob_lo[1])*0.5;
            z0 = prob_lo[2] + 0.25*(prob_hi[0]-prob_lo[0]) + amrex::Random()*(prob_hi[2]-prob_lo[2])*0.5;

            SetPosition(1, x0, y0, z0);

            costheta = 2.*amrex::Random() - 1.;
            sintheta = 1. - costheta*costheta;

            Real phi = amrex::Random() * 2. * 3.14159265358979;
            cosphi = std::cos(phi);
            sinphi = std::sin(phi);

            x1 = x0 + dr*sintheta*cosphi;
            y1 = y0 + dr*sintheta*sinphi;
            z1 = z0 + dr*costheta;

            SetPosition(2, x1, y1, z1);
        }

        if(setMeasureFinal == 1)
        {
            ParticleType & part0 = particles[0];
            ParticleType & part1 = particles[1];

            Real ee = (permittivity*4*3.14159265359);

            Real re = threepmCurrentBin*(threepmRange/threepmBins)/(1);

            Real forceNorm = -(dx[0]*dx[0])*ee/(part0.rdata(FHD_realData::q)*part1.rdata(FHD_realData::q));

            //Print() << "Force: " << part0.rdata(FHD_realData::forcex) << ", " << part0.rdata(FHD_realData::forcey) << ", " << part0.rdata(FHD_realData::forcez) << std::endl;

            Real forceMag = sqrt(pow(part0.rdata(FHD_realData::forcex),2) + pow(part0.rdata(FHD_realData::forcey),2) + pow(part0.rdata(FHD_realData::forcez),2))*forceNorm;

           // Print() << "ForceMag: " << forceNorm << std::endl;

           // Print() << "currentPre: " << threepmVals[threepmCurrentBin] << std::endl;

           // Print() << "Bin: " << threepmCurrentBin << std::endl;

            if(threepmCurrentSample == 1)
            {
                threepmVals[threepmCurrentBin] = 0;
            }

            threepmVals[threepmCurrentBin] = (threepmVals[threepmCurrentBin]*(threepmCurrentSample - 1) + forceMag)/threepmCurrentSample;

          //  Print() << "currentPost: " << threepmVals[threepmCurrentBin] << std::endl;

         //   Print() << "Sample: " << threepmCurrentSample << std::endl;

            if(forceMag < threepmMin[threepmCurrentBin])
            {
                threepmMin[threepmCurrentBin] = forceMag;
            }

            if(forceMag > threepmMax[threepmCurrentBin])
            {
                threepmMax[threepmCurrentBin] = forceMag;
            }

            Print() << "Bin " << threepmCurrentBin << " norm: " << threepmVals[threepmCurrentBin] << endl;

            threepmCurrentSample++;

            if(threepmCurrentSample > threepmSamples)
            {
                threepmCurrentSample = 1;
                threepmCurrentBin++;
            }
        }
    }

    if(setMeasureFinal == 2)
    {
        Print() << "Outputting correction data\n";
        std::string filename = "threepmPoints";
        std::ofstream ofs0(filename, std::ofstream::out);

        // normalize by
        for(int i=0;i<threepmBins;i++) {
            ofs0 << threepmPoints[i] << ", ";
        }
        ofs0 << std::endl;
        ofs0.close();

        filename = "threepmMax";
        std::ofstream ofs1(filename, std::ofstream::out);

        // normalize by
        for(int i=0;i<threepmBins;i++) {
            ofs1 << threepmMax[i] << ", ";
        }
        ofs1 << std::endl;
        ofs1.close();

        filename = "threepmMin";
        std::ofstream ofs2(filename, std::ofstream::out);

        // normalize by
        for(int i=0;i<threepmBins;i++) {
            ofs2 << threepmMin[i] << ", ";
        }
        ofs2 << std::endl;
        ofs2.close();

        filename = "threepmVals";
        std::ofstream ofs3(filename, std::ofstream::out);

        // normalize by
        for(int i=0;i<threepmBins;i++) {
            ofs3 << threepmVals[i] << ", ";
        }
        ofs3 << std::endl;
        ofs3.close();

        threepmCurrentBin++;
    }
}

//void
//FhdParticleContainer::correctCellVectors(int old_index, int new_index,
//                                         int grid, const ParticleType& p)
//{
//    if (not p.idata(FHD_intData::sorted)) return;
//    IntVect iv(p.idata(FHD_intData::i), p.idata(FHD_intData::j), p.idata(FHD_intData::k));
//
//    auto& cell_vector = m_cell_vectors[grid](iv);
//    for (int i = 0; i < static_cast<int>(cell_vector.size()); ++i) {
//        if (cell_vector[i] == old_index + 1) {
//            cell_vector[i] = new_index + 1;
//            return;
//        }
//    }
//}

int
FhdParticleContainer::numWrongCell()
{
    BL_PROFILE_VAR("numWrongCell()",numWrongCell);

    const int lev = 0;
    int num_wrong = 0;

#ifdef _OPENMP
#pragma omp parallel reduction(+:num_wrong)
#endif
    for (FhdParIter pti(*this, lev); pti.isValid(); ++pti)
    {
        auto& particles = pti.GetArrayOfStructs();
        const int np    = pti.numParticles();
        for(int pindex = 0; pindex < np; ++pindex) {
            const ParticleType& p = particles[pindex];
            const IntVect& iv = this->Index(p, lev);
            if ((iv[0] != p.idata(FHD_intData::i)) or (iv[1] != p.idata(FHD_intData::j)) or (iv[2] != p.idata(FHD_intData::k))) {
                num_wrong += 1;
            }
        }
    }

    ParallelDescriptor::ReduceIntSum(num_wrong, ParallelDescriptor::IOProcessorNumber());
    return num_wrong;
}

void FhdParticleContainer::PostRestart()
{
    BL_PROFILE_VAR("PostRestart()",PostRestart);

}
