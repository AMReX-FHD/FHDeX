#include "INS_functions.H"
#include "common_functions.H"
#include "FhdParticleContainer.H"

void getCellVols(MultiFab & vols, const Geometry & Geom, int samples)
{

    const Real* dx = Geom.CellSize();
    const Real* plo = Geom.ProbLo();

    for ( MFIter mfi(vols); mfi.isValid(); ++mfi ) {

        const Box& bx = mfi.validbox();
        get_cell_vols(BL_TO_FORTRAN_3D(vols[mfi]), ZFILL(dx), &samples, plo);

    }
}



//Anthony Stuff

void FhdParticleContainer::InitParticlesDSMC(species* particleInfo, int pL, int pR, Real tL, Real tR) {
    //distribute particles in two boxes, assuming domain is split in center.
    //Place specified number of particles at specified temp in each box.

    const int lev = 0;
    const Geometry& geom = Geom(lev);
    const Real* dx = geom.CellSize();
    const Real* plo = geom.ProbLo();
    const Real* phi = geom.ProbHi();

    //constants used to randomly distribute particles in each box
    const Real xMinL = plo[0];
    const Real xMinR = (plo[0]+phi[0])/2.0;
    const Real xRange = (phi[0]-plo[0])/2.0;
    Real xMin, temp;

    //printf("Procs: %d\n", ParallelDescriptor::NProcs());

    int pcount = 0;

    for (MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        const Box& tile_box  = mfi.tilebox();
        const RealBox tile_realbox{tile_box, geom.CellSize(), geom.ProbLo()};
        const int grid_id = mfi.index();
        const int tile_id = mfi.LocalTileIndex();
        auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];

        //Assuming tile=box for now, i.e. no tiling.
        IntVect smallEnd = tile_box.smallEnd();
        IntVect bigEnd = tile_box.bigEnd();


        for(int i_spec=0; i_spec < nspecies; i_spec++)
        {
            for (int i_part=0; i_part<particleInfo[i_spec].total;i_part++)
            {

                //set coefficients based on if this particle is going in left or right box
                if (i_part < pL) {
                    xMin = xMinL; temp = tL;
                }
                else{
                    xMin = xMinR; temp = tR;
                }

                ParticleType p;
                p.id()  = ParticleType::NextID();
                p.cpu() = ParallelDescriptor::MyProc();
                p.idata(IntData::sorted) = 0;

                p.pos(0) = xMin + get_uniform_func()*xRange;
                printf("%d %f\n", i_part,p.pos(0));
                p.pos(1) = plo[1] + get_uniform_func()*(phi[1]-plo[1]);
#if (BL_SPACEDIM == 3)
                p.pos(2) = plo[2] + get_uniform_func()*(phi[2]-plo[2]);
#endif


                p.rdata(RealData::q) = 0;

                //original position stored for MSD calculations
                p.rdata(RealData::ox) = p.pos(0);
                p.rdata(RealData::oy) = p.pos(1);
#if (BL_SPACEDIM == 3)
                p.rdata(RealData::oz) = p.pos(2);
#endif

                p.rdata(RealData::vx) = sqrt(particleInfo[i_spec].R*temp)*get_particle_normal_func();
                p.rdata(RealData::vy) = sqrt(particleInfo[i_spec].R*temp)*get_particle_normal_func();
                p.rdata(RealData::vz) = sqrt(particleInfo[i_spec].R*temp)*get_particle_normal_func();

                p.rdata(RealData::fx) = 0;
                p.rdata(RealData::fy) = 0;
                p.rdata(RealData::fz) = 0;

                p.rdata(RealData::ux) = 0;
                p.rdata(RealData::uy) = 0;
                p.rdata(RealData::uz) = 0;

                //Print() << "Pos: " << p.pos(0) << ", " << p.pos(1) << ", " << p.pos(2) << p.rdata(RealData::vx) << ", " << p.rdata(RealData::vy) << ", " << p.rdata(RealData::vz) << "\n" ;

                p.rdata(RealData::mass) = particleInfo[i_spec].m; //mass
                p.rdata(RealData::R) = particleInfo[i_spec].R; //R
                p.rdata(RealData::radius) = particleInfo[i_spec].d/2.0; //radius

                p.idata(IntData::species) = i_spec +1;

                p.rdata(RealData::ox) = 0;
                p.rdata(RealData::oy) = 0;
                p.rdata(RealData::oz) = 0;

                p.rdata(RealData::ax) = 0;
                p.rdata(RealData::ay) = 0;
                p.rdata(RealData::az) = 0;

                p.rdata(RealData::accelFactor) = 0;
                p.rdata(RealData::dragFactor) = 0;
                p.rdata(RealData::travelTime) = 0;
                p.rdata(RealData::diffAv) = 0;
                p.rdata(RealData::stepCount) = 0;
                p.rdata(RealData::multi) = 0;
                p.rdata(RealData::dryDiff) = 0;
                p.rdata(RealData::wetDiff) = 0;
                p.rdata(RealData::totalDiff) = 0;
                p.rdata(RealData::sigma) = 0;
                p.rdata(RealData::eepsilon) = 0;
                p.rdata(RealData::potential) = 0;

                particle_tile.push_back(p);

                pcount++;
            }
//
        }
    }

    Redistribute();
    UpdateCellVectors();
    ReBin();
}

void FhdParticleContainer::getParticleMapping(std::map<int,int>& particleMap, int xCells, int num_boxes, int pL, int pR) {
    //get mapping from tile number to number of particles in that tile

    //get geometry to iterate over
    const int lev = 0;
    const Geometry& geom = Geom(lev);

    //variables for particle distribution
    int box_per_side = num_boxes / 2; //assumes evenly split, syymetric domain
    int ppbL = pL / box_per_side; int ppbR = pR / box_per_side; //particles per box
    int remL = pL % box_per_side; int remR = pR % box_per_side; //remainder per side
    int countL = 0; int countR = 0; //count how many boxes considered on each side

    //printf("%d %d\n", pL, pR);
    //printf("%d %d %d %d\n", ppbL, ppbR, remL, remR);

    //deal with random distribution of the remainder - overcounts for many processes
    std::vector<int> boxL; std::vector<int> boxR;
    for (int i = 0; i < box_per_side; i++) {
        boxL.push_back(i); boxR.push_back(i);
    }
    std::random_shuffle( boxL.begin(), boxL.end() );
    std::random_shuffle( boxR.begin(), boxR.end() );
    int addLeft[box_per_side]; int addRight[box_per_side];
    for (int i = 0; i < box_per_side; i++) {
        addLeft[i] = 0; addRight[i] = 0;
    }
    for (int i = 0; i < remL; i++) addLeft[boxL[i]] = 1;
    for (int i = 0; i < remR; i++) addRight[boxR[i]] = 1;


    //loop over boxes distributing particles
    for (MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        const Box& tile_box  = mfi.tilebox();
        const RealBox tile_realbox{tile_box, geom.CellSize(), geom.ProbLo()};
        const int grid_id = mfi.index();
        const int tile_id = mfi.LocalTileIndex();

        //Assuming tile=box for now, i.e. no tiling.
        IntVect smallEnd = tile_box.smallEnd();
        IntVect bigEnd = tile_box.bigEnd();

        //printf("s, b, id = %d %d %d\n", smallEnd[0], bigEnd[0], grid_id);

        //check which side, place particles
        int add;
        if (smallEnd[0] < xCells / 2) {
            add = ppbL+addLeft[countL];
            particleMap.insert(std::pair<int,int>(grid_id, add));
            countL++;
        }
        else {
            add = ppbR+addRight[countR];
            particleMap.insert(std::pair<int,int>(grid_id, add));
            countR++;
        }

    }

    //for (int i = 0; i < 8; i++) printf("%d %d\n", i, particleMap[i]);
}

void FhdParticleContainer::InitParticlesDSMCtest(species* particleInfo, int num_boxes, int pL, int pR, Real tL, Real tR) {
    //distribute particles in two boxes, assuming domain is split in center.
    //Place specified number of particles at specified temp in each box.

    const int lev = 0;
    const Geometry& geom = Geom(lev);
    const Real* dx = geom.CellSize();
    const Real* plo = geom.ProbLo();
    const Real* phi = geom.ProbHi();

    //printf("%f \n", phi[0]);

    const int xCells = (int)(phi[0]/dx[0]);

    //constants used to randomly distribute particles in each box
    Real temp;

    //printf("Procs: %d\n", ParallelDescriptor::NProcs());

    int pcount = 0;
    std::map<int,int> particleMap;
    getParticleMapping(particleMap, xCells, num_boxes, pL, pR);
    //MPI_Barrier(MPI_COMM_WORLD);
    ParallelDescriptor::Barrier();

    for (MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        const Box& tile_box  = mfi.tilebox();
        const RealBox tile_realbox{tile_box, geom.CellSize(), geom.ProbLo()};
        const int grid_id = mfi.index();
        const int tile_id = mfi.LocalTileIndex();
        auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];

        //Assuming tile=box for now, i.e. no tiling.
        IntVect smallEnd = tile_box.smallEnd();
        IntVect bigEnd = tile_box.bigEnd();

        //get bounds and ranges for particle generation
        const Real xRange = (bigEnd[0]-smallEnd[0]+1)*dx[0];
        const Real xMin = smallEnd[0]*dx[0];
        const Real yRange = (bigEnd[1]-smallEnd[1]+1)*dx[1];
        const Real yMin = smallEnd[1]*dx[1];
#if (BL_SPACEDIM == 3)
        const Real zRange = (bigEnd[2]-smallEnd[2]+1)*dx[2];
        const Real zMin = smallEnd[2]*dx[2];
#endif

        //printf("s, b, id = %d %d %d\n", smallEnd[0], bigEnd[0], grid_id);

        for(int i_spec=0; i_spec < nspecies; i_spec++)
        {
            for (int i_part=0; i_part<particleMap[grid_id];i_part++)
            {
                //set coefficients based on if this particle is going in left or right box
                if (smallEnd[0] < xCells / 2) {
                    temp = tL;
                }
                else{
                    temp = tR;
                }

                ParticleType p;
                p.id()  = ParticleType::NextID();
                p.cpu() = ParallelDescriptor::MyProc();
                p.idata(IntData::sorted) = 0;

                p.pos(0) = xMin + get_uniform_func()*xRange;
                p.pos(1) = yMin + get_uniform_func()*yRange;
#if (BL_SPACEDIM == 3)
                p.pos(2) = zMin + get_uniform_func()*zRange;
#endif
                printf("%d %f %f\n", i_part,p.pos(0), p.pos(1));

                p.rdata(RealData::q) = 0;

                //original position stored for MSD calculations
                p.rdata(RealData::ox) = p.pos(0);
                p.rdata(RealData::oy) = p.pos(1);
#if (BL_SPACEDIM == 3)
                p.rdata(RealData::oz) = p.pos(2);
#endif

                p.rdata(RealData::vx) = sqrt(particleInfo[i_spec].R*temp)*get_particle_normal_func();
                p.rdata(RealData::vy) = sqrt(particleInfo[i_spec].R*temp)*get_particle_normal_func();
                p.rdata(RealData::vz) = sqrt(particleInfo[i_spec].R*temp)*get_particle_normal_func();

                p.rdata(RealData::fx) = 0;
                p.rdata(RealData::fy) = 0;
                p.rdata(RealData::fz) = 0;

                p.rdata(RealData::ux) = 0;
                p.rdata(RealData::uy) = 0;
                p.rdata(RealData::uz) = 0;

                //Print() << "Pos: " << p.pos(0) << ", " << p.pos(1) << ", " << p.pos(2) << p.rdata(RealData::vx) << ", " << p.rdata(RealData::vy) << ", " << p.rdata(RealData::vz) << "\n" ;

                p.rdata(RealData::mass) = particleInfo[i_spec].m; //mass
                p.rdata(RealData::R) = particleInfo[i_spec].R; //R
                p.rdata(RealData::radius) = particleInfo[i_spec].d/2.0; //radius

                p.idata(IntData::species) = i_spec +1;

                p.rdata(RealData::ox) = 0;
                p.rdata(RealData::oy) = 0;
                p.rdata(RealData::oz) = 0;

                p.rdata(RealData::ax) = 0;
                p.rdata(RealData::ay) = 0;
                p.rdata(RealData::az) = 0;

                p.rdata(RealData::accelFactor) = 0;
                p.rdata(RealData::dragFactor) = 0;
                p.rdata(RealData::travelTime) = 0;
                p.rdata(RealData::diffAv) = 0;
                p.rdata(RealData::stepCount) = 0;
                p.rdata(RealData::multi) = 0;
                p.rdata(RealData::dryDiff) = 0;
                p.rdata(RealData::wetDiff) = 0;
                p.rdata(RealData::totalDiff) = 0;
                p.rdata(RealData::sigma) = 0;
                p.rdata(RealData::eepsilon) = 0;
                p.rdata(RealData::potential) = 0;

                particle_tile.push_back(p);

                pcount++;
            }
//
        }
    }

    Redistribute();
    UpdateCellVectors();
    ReBin();
}

void FhdParticleContainer::ApplyThermostat(species* particleInfo, MultiFab& cellVols,
                                           paramPlane* paramplanes, int ns, Real tL, Real tR) {
    //compute the cell temperatures as variance of velocity. Re-scale to given values

    const int lev = 0;
    const Real Neff = particleInfo->Neff;
    Real vFix; //generate new velocities if 0 variance
    Real varTol = 1e-5; //tolerance for variance to be 0

    //declare storage for variables going to fortran
    int  pL = 0, pR = 0;
    Real varL = 0, varR = 0;
    Real vLx = 0, vRx = 0;
    Real vLy = 0, vRy = 0;
    Real vLz = 0, vRz = 0;
    Real meanLx = 0; Real meanRx = 0;
    Real meanLy = 0; Real meanRy = 0;
    Real meanLz = 0; Real meanRz = 0;

    //get the total particle number and velocity on each side
    for (FhdParIter pti(*this, lev); pti.isValid(); ++pti) {
        const int grid_id = pti.index();
        const int tile_id = pti.LocalTileIndex();
        const Box& tile_box  = pti.tilebox();

        auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];
        auto& parts = particle_tile.GetArrayOfStructs();
        const int Np = parts.numParticles();

        //printf("Get velocity, np = %d\n", Np);
        getVelocity(parts.data(),
                    ARLIM_3D(tile_box.loVect()),
                    ARLIM_3D(tile_box.hiVect()),
                    m_vector_ptrs[grid_id].dataPtr(),
                    m_vector_size[grid_id].dataPtr(),
                    ARLIM_3D(m_vector_ptrs[grid_id].loVect()),
                    ARLIM_3D(m_vector_ptrs[grid_id].hiVect()),
                    BL_TO_FORTRAN_3D(cellVols[pti]), &Neff, &Np,
                    paramplanes, &ns, &pL, &pR, &vLx, &vRx, &vLy, &vRy,
                    &vLz, &vRz);


    }

    //perform reductions
    ParallelDescriptor::ReduceRealSum(vLx);
    ParallelDescriptor::ReduceRealSum(vRx);
    ParallelDescriptor::ReduceRealSum(vLy);
    ParallelDescriptor::ReduceRealSum(vRy);
    ParallelDescriptor::ReduceRealSum(vLz);
    ParallelDescriptor::ReduceRealSum(vRz);
    ParallelDescriptor::ReduceIntSum(pL);
    ParallelDescriptor::ReduceIntSum(pR);


    //debug line
    //printf("Particle counts: %d %d\n", pL, pR);

    //get mean velocities per side
    if (pL > 0) {
        meanLx = vLx / pL;
        meanLy = vLy / pL;
        meanLz = vLz / pL;
    }
    if (pR > 0) {
        meanRx = vRx / pR;
        meanRy = vRy / pR;
        meanRz = vRz / pR;
    }


    //debug line
    //printf("%f %f %f %f %f %f\n", meanLx, meanRx, meanLy, meanRy, meanLz, meanRz);

    Real netV = amrex::Math::abs(meanLx) + amrex::Math::abs(meanRx) + amrex::Math::abs(meanLy) + amrex::Math::abs(meanRy) +
                amrex::Math::abs(meanLz) + amrex::Math::abs(meanRz);


    //if netV is 0, no thermostatting is necc
    if (netV < varTol) {
        return;
    }


    //get temperature on each side via velocity variance
    for (FhdParIter pti(*this, lev); pti.isValid(); ++pti) {
        const int grid_id = pti.index();
        const int tile_id = pti.LocalTileIndex();
        const Box& tile_box  = pti.tilebox();

        auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];
        auto& parts = particle_tile.GetArrayOfStructs();
        const int Np = parts.numParticles();

        getTemp(parts.data(),
                ARLIM_3D(tile_box.loVect()),
                ARLIM_3D(tile_box.hiVect()),
                m_vector_ptrs[grid_id].dataPtr(),
                m_vector_size[grid_id].dataPtr(),
                ARLIM_3D(m_vector_ptrs[grid_id].loVect()),
                ARLIM_3D(m_vector_ptrs[grid_id].hiVect()),
                BL_TO_FORTRAN_3D(cellVols[pti]), &Neff, &Np,
                paramplanes, &ns, &meanLx, &meanRx, &meanLy, &meanRy,
                &meanLz, &meanRz, &varL, &varR);

    }

    //perform reductions, get variance by scaling by population
    ParallelDescriptor::ReduceRealSum(varL); varL = varL/(3*pL);
    ParallelDescriptor::ReduceRealSum(varR); varR = varR/(3*pR);

    //make variances 0 in case of no particles
    if (pL == 0) varL = 0;
    if (pR == 0) varR = 0;

    //debug line
    //printf("Vars: %f %f\n", varL, varR);

    //if the variance is 0, regenerate those velocities and re-thermostat
    if (varL < varTol) {
        if (pL > 1) {
            for (FhdParIter pti(*this, lev); pti.isValid(); ++pti) {
                const int grid_id = pti.index();
                const int tile_id = pti.LocalTileIndex();
                const Box& tile_box  = pti.tilebox();

                auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];
                auto& parts = particle_tile.GetArrayOfStructs();
                const int Np = parts.numParticles();

                for (int i = 0; i < Np; i++) {
                    if (parts[i].pos(0) < 1) {
                        vFix = sqrt(particleInfo[0].R*tL)*get_particle_normal_func();
                        parts[i].rdata(RealData::vx) = vFix;
                        vFix = sqrt(particleInfo[0].R*tL)*get_particle_normal_func();
                        parts[i].rdata(RealData::vy) = vFix;
                        vFix = sqrt(particleInfo[0].R*tL)*get_particle_normal_func();
                        parts[i].rdata(RealData::vz) = vFix;
                    }
                }
            }

            ApplyThermostat(particleInfo, cellVols,paramplanes, ns, tL, tR);
            return;
        }
        else {//only 1 particle, set velocity by energy E=1/2mv^2=3/2kT
            for (FhdParIter pti(*this, lev); pti.isValid(); ++pti) {
                const int grid_id = pti.index();
                const int tile_id = pti.LocalTileIndex();
                const Box& tile_box  = pti.tilebox();

                auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];
                auto& parts = particle_tile.GetArrayOfStructs();
                const int Np = parts.numParticles();

                for (int i = 0; i < Np; i++) {
                    if (parts[i].pos(0) < 1) {
                        Real v2 = parts[i].rdata(RealData::vx)*parts[i].rdata(RealData::vx)+
                                  parts[i].rdata(RealData::vy)*parts[i].rdata(RealData::vy)+
                                  parts[i].rdata(RealData::vz)*parts[i].rdata(RealData::vz);
                        Real C = sqrt(3*particleInfo[0].R*tL/v2); //correction factor
                        parts[i].rdata(RealData::vx) *= C;
                        parts[i].rdata(RealData::vy) *= C;
                        parts[i].rdata(RealData::vz) *= C;
                    }
                }
                Real rC = sqrt(tR/varR);
                Real lC = 1.0;
                Real mean0 = 0;

                thermostat(parts.data(),
                        ARLIM_3D(tile_box.loVect()),
                        ARLIM_3D(tile_box.hiVect()),
                        m_vector_ptrs[grid_id].dataPtr(),
                        m_vector_size[grid_id].dataPtr(),
                        ARLIM_3D(m_vector_ptrs[grid_id].loVect()),
                        ARLIM_3D(m_vector_ptrs[grid_id].hiVect()),
                        BL_TO_FORTRAN_3D(cellVols[pti]), &Neff, &Np,
                        paramplanes, &ns, &mean0, &meanRx, &mean0, &meanRy,
                        &mean0, &meanRz, &lC, &rC);
            }
            return;
        }
    }
    if (varR < varTol) {
        if (pR > 1) {
            for (FhdParIter pti(*this, lev); pti.isValid(); ++pti) {
                const int grid_id = pti.index();
                const int tile_id = pti.LocalTileIndex();
                const Box& tile_box  = pti.tilebox();

                auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];
                auto& parts = particle_tile.GetArrayOfStructs();
                const int Np = parts.numParticles();

                for (int i = 0; i < Np; i++) {
                    if (parts[i].pos(0) > 1) {
                        vFix = sqrt(particleInfo[0].R*tR)*get_particle_normal_func();
                        parts[i].rdata(RealData::vx) = vFix;
                        vFix = sqrt(particleInfo[0].R*tR)*get_particle_normal_func();
                        parts[i].rdata(RealData::vy) = vFix;
                        vFix = sqrt(particleInfo[0].R*tR)*get_particle_normal_func();
                        parts[i].rdata(RealData::vz) = vFix;
                    }
                }
            }
            //abort();
            ApplyThermostat(particleInfo, cellVols,paramplanes, ns, tL, tR);
            return;
        }
        else {//only 1 particle in right, set velocity by energy E=1/2mv^2=3/2kT
            for (FhdParIter pti(*this, lev); pti.isValid(); ++pti) {
                const int grid_id = pti.index();
                const int tile_id = pti.LocalTileIndex();
                const Box& tile_box  = pti.tilebox();

                auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];
                auto& parts = particle_tile.GetArrayOfStructs();
                const int Np = parts.numParticles();

                for (int i = 0; i < Np; i++) {
                    if (parts[i].pos(0) > 1) {
                        Real v2 = parts[i].rdata(RealData::vx)*parts[i].rdata(RealData::vx)+
                                  parts[i].rdata(RealData::vy)*parts[i].rdata(RealData::vy)+
                                  parts[i].rdata(RealData::vz)*parts[i].rdata(RealData::vz);
                        Real C = sqrt(3*particleInfo[0].R*tR/v2); //correction factor

                        parts[i].rdata(RealData::vx) *= C;
                        parts[i].rdata(RealData::vy) *= C;
                        parts[i].rdata(RealData::vz) *= C;
                    }
                }

                Real lC = sqrt(tL/varL);
                Real rC = 1.0;
                Real mean0 = 0;

                thermostat(parts.data(),
                        ARLIM_3D(tile_box.loVect()),
                        ARLIM_3D(tile_box.hiVect()),
                        m_vector_ptrs[grid_id].dataPtr(),
                        m_vector_size[grid_id].dataPtr(),
                        ARLIM_3D(m_vector_ptrs[grid_id].loVect()),
                        ARLIM_3D(m_vector_ptrs[grid_id].hiVect()),
                        BL_TO_FORTRAN_3D(cellVols[pti]), &Neff, &Np,
                        paramplanes, &ns, &meanLx, &mean0, &meanLy, &mean0,
                        &meanLz, &mean0, &lC, &rC);
            }
            return;
        }

    }

    //only apply corrections if variance is nonzero
    if (varL > varTol && varR > varTol) {
        //compute correction factors

        //printf("Vars after conditional: %f %f\n", varL, varR);
        Real lC = sqrt(tL/varL);
        Real rC = sqrt(tR/varR);

        //printf("corrections: %f %f\n", lC, rC);

        //apply correction factors
        for (FhdParIter pti(*this, lev); pti.isValid(); ++pti) {
            const int grid_id = pti.index();
            const int tile_id = pti.LocalTileIndex();
            const Box& tile_box  = pti.tilebox();

            auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];
            auto& parts = particle_tile.GetArrayOfStructs();
            const int Np = parts.numParticles();


            thermostat(parts.data(),
                        ARLIM_3D(tile_box.loVect()),
                        ARLIM_3D(tile_box.hiVect()),
                        m_vector_ptrs[grid_id].dataPtr(),
                        m_vector_size[grid_id].dataPtr(),
                        ARLIM_3D(m_vector_ptrs[grid_id].loVect()),
                        ARLIM_3D(m_vector_ptrs[grid_id].hiVect()),
                        BL_TO_FORTRAN_3D(cellVols[pti]), &Neff, &Np,
                        paramplanes, &ns, &meanLx, &meanRx, &meanLy, &meanRy,
                        &meanLz, &meanRz, &lC, &rC);

        }
    }
}


void FhdParticleContainer::Resample(species* particleInfo, Real tL, Real tR) {
    //regenerate all particle positions and velocities, conserve box counts

    const int lev = 0;
    const Real Neff = particleInfo->Neff;

    const Geometry& geom = Geom(lev);
    const Real* dx = geom.CellSize();
    const Real* plo = geom.ProbLo();
    const Real* phi = geom.ProbHi();

    Real temp;

    //loop over boxes
    for (FhdParIter pti(*this, lev); pti.isValid(); ++pti) {
        const int grid_id = pti.index();
        const int tile_id = pti.LocalTileIndex();
        const Box& tile_box  = pti.tilebox();

        //Assuming tile=box for now, i.e. no tiling.
        IntVect smallEnd = tile_box.smallEnd();
        IntVect bigEnd = tile_box.bigEnd();

        //get bounds and ranges for particle generation
        const Real xRange = (bigEnd[0]-smallEnd[0]+1)*dx[0];
        const Real xMin = smallEnd[0]*dx[0];
        const Real yRange = (bigEnd[1]-smallEnd[1]+1)*dx[1];
        const Real yMin = smallEnd[1]*dx[1];
#if (BL_SPACEDIM == 3)
        const Real zRange = (bigEnd[2]-smallEnd[2]+1)*dx[2];
        const Real zMin = smallEnd[2]*dx[2];
#endif

        //get particle array
        auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];
        auto& parts = particle_tile.GetArrayOfStructs();
        const int Np = parts.numParticles();

        //loop over particles
        for (int i = 0; i < Np; i++) {
            if (parts[i].pos(0) < 1) {
                temp = tL;
            }
            else {
                temp = tR;
            }

            //regenerate positions
            parts[i].pos(0) = xMin + get_uniform_func()*xRange;
            parts[i].pos(1) = yMin + get_uniform_func()*yRange;
#if (BL_SPACEDIM == 3)
            parts[i].pos(2) = zMin + get_uniform_func()*zRange;
#endif


            //regenerate velocities
            parts[i].rdata(RealData::vx) = sqrt(particleInfo[0].R*temp)*get_particle_normal_func();
            parts[i].rdata(RealData::vy) = sqrt(particleInfo[0].R*temp)*get_particle_normal_func();
            parts[i].rdata(RealData::vz) = sqrt(particleInfo[0].R*temp)*get_particle_normal_func();


        }
    }

    Redistribute();
    UpdateCellVectors();
    ReBin();
}
