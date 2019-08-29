#include "INS_functions.H"
#include "common_functions.H"
#include "DsmcParticleContainer.H"
#include <common_namespace.H>

#include "species.H"

using namespace common;
void DsmcParticleContainer::InitParticles(species* particleInfo, const Real* dxp)
{
    
    const int lev = 0;
    const Geometry& geom = Geom(lev);
    const Real* dx = geom.CellSize();
    const Real* plo = geom.ProbLo();
    const Real* phi = geom.ProbHi(); 
        double rad;

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


        if(ParallelDescriptor::MyProc() == 0 && mfi.LocalTileIndex() == 0)
        {

            for(int i_spec=0; i_spec < nspecies; i_spec++)
            {

            for (int i_part=0; i_part<particleInfo[i_spec].total;i_part++)
            {
                ParticleType p;
                p.id()  = ParticleType::NextID();
                p.cpu() = ParallelDescriptor::MyProc();
                p.idata(DSMC_intData::sorted) = 0;

               
//                p.pos(0) = plo[0] + get_uniform_func()*(phi[0]-plo[0]);
//                p.pos(1) = plo[1] + get_uniform_func()*(phi[1]-plo[1]);
//#if (BL_SPACEDIM == 3)
//                p.pos(2) = plo[2] + get_uniform_func()*(phi[2]-plo[2]);
//#endif
                rad = 1;
                while(rad > 0.99*(2.5e-4))
                {                
                    p.pos(0) = plo[0] + get_uniform_func()*(phi[0]-plo[0]);
                    p.pos(1) = plo[1] + get_uniform_func()*(phi[1]-plo[1]);
#if (BL_SPACEDIM == 3)
                    p.pos(2) = plo[2] + get_uniform_func()*(phi[2]-plo[2]);
#endif
                    rad = sqrt((p.pos(0)-0)*(p.pos(0)-0) + (p.pos(1)-0)*(p.pos(1)-0));
                }
                
                p.rdata(DSMC_realData::q) = 0;

                //original position stored for MSD calculations
                p.rdata(DSMC_realData::ox) = p.pos(0);
                p.rdata(DSMC_realData::oy) = p.pos(1);
#if (BL_SPACEDIM == 3)
                p.rdata(DSMC_realData::oz) = p.pos(2);
#endif

                p.rdata(DSMC_realData::vx) = sqrt(particleInfo[i_spec].R*particleInfo[i_spec].T)*get_particle_normal_func();
                p.rdata(DSMC_realData::vy) = sqrt(particleInfo[i_spec].R*particleInfo[i_spec].T)*get_particle_normal_func();
                p.rdata(DSMC_realData::vz) = sqrt(particleInfo[i_spec].R*particleInfo[i_spec].T)*get_particle_normal_func();

                p.rdata(DSMC_realData::fx) = 0;
                p.rdata(DSMC_realData::fy) = 0;
                p.rdata(DSMC_realData::fz) = 0;

                p.rdata(DSMC_realData::ux) = 0;
                p.rdata(DSMC_realData::uy) = 0;
                p.rdata(DSMC_realData::uz) = 0;

                //Print() << "Pos: " << p.pos(0) << ", " << p.pos(1) << ", " << p.pos(2) << p.rdata(DSMC_realData::vx) << ", " << p.rdata(DSMC_realData::vy) << ", " << p.rdata(DSMC_realData::vz) << "\n" ;

                p.rdata(DSMC_realData::mass) = particleInfo[i_spec].m; //mass
                p.rdata(DSMC_realData::R) = particleInfo[i_spec].R; //R
                p.rdata(DSMC_realData::radius) = particleInfo[i_spec].d/2.0; //radius

                p.idata(DSMC_intData::species) = i_spec +1;

                p.rdata(DSMC_realData::ox) = 0;
                p.rdata(DSMC_realData::oy) = 0;
                p.rdata(DSMC_realData::oz) = 0;

                p.rdata(DSMC_realData::ax) = 0;
                p.rdata(DSMC_realData::ay) = 0;
                p.rdata(DSMC_realData::az) = 0;

                p.rdata(DSMC_realData::accelFactor) = 0;
                p.rdata(DSMC_realData::dragFactor) = 0;
                p.rdata(DSMC_realData::travelTime) = 0;
                p.rdata(DSMC_realData::diffAv) = 0;
                p.rdata(DSMC_realData::stepCount) = 0;
                p.rdata(DSMC_realData::multi) = 0;
                p.rdata(DSMC_realData::dryDiff) = 0;
                p.rdata(DSMC_realData::wetDiff) = 0;
                p.rdata(DSMC_realData::totalDiff) = 0;
                p.rdata(DSMC_realData::sigma) = 0;
                p.rdata(DSMC_realData::eepsilon) = 0;
                p.rdata(DSMC_realData::potential) = 0;
                
                particle_tile.push_back(p);

                pcount++;
            }
            }
//           
        }
    }

    UpdateCellVectors();
    ReBin();
    Redistribute();
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
