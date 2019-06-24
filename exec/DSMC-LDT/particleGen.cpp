#include "INS_functions.H"
#include "common_functions.H"
#include "FhdParticleContainer.H"

void FhdParticleContainer::InitParticles(species* particleInfo)
{
    
    const int lev = 0;
    const Geometry& geom = Geom(lev);
    const Real* dx = geom.CellSize();
    const Real* plo = geom.ProbLo();
    const Real* phi = geom.ProbHi();  

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
            for (int i_part=0; i_part<particleInfo[i_spec].ppb;i_part++)
            {
                ParticleType p;
                p.id()  = ParticleType::NextID();
                p.cpu() = ParallelDescriptor::MyProc();
                p.idata(IntData::sorted) = 0;

               
                p.pos(0) = plo[0] + get_uniform_func()*(phi[0]-plo[0]);
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

                p.rdata(RealData::vx) = sqrt(particleInfo[i_spec].R*particleInfo[i_spec].T)*get_particle_normal_func();
                p.rdata(RealData::vy) = sqrt(particleInfo[i_spec].R*particleInfo[i_spec].T)*get_particle_normal_func();
                p.rdata(RealData::vz) = sqrt(particleInfo[i_spec].R*particleInfo[i_spec].T)*get_particle_normal_func();

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
                p.rdata(RealData::vy) = sqrt(particleInfo[i_spec].R*0)*get_particle_normal_func();
                p.rdata(RealData::vz) = sqrt(particleInfo[i_spec].R*0)*get_particle_normal_func();

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
                                           surface* surfaces, int ns, Real tL, Real tR) {
    //compute the cell temperatures as variance of velocity. Re-scale to given values

    const int lev = 0;
    const Real Neff = particleInfo->Neff;   

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
                    surfaces, &ns, &tL, &tR);

    }
}

