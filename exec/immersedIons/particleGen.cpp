#include "INS_functions.H"
#include "common_functions.H"
#include "FhdParticleContainer.H"
#include <sstream>
#include <string>
#include <fstream>

void FhdParticleContainer::InitParticles(species* particleInfo, const Real* dxp)
{
    const int lev = 0;
    const Geometry& geom = Geom(lev);

    int pcount = 0;

    bool proc0_enter = true;
    int pinnedParticles = 0;

    for (MFIter mfi = MakeMFIter(lev, true); mfi.isValid(); ++mfi) {

        const Box& tile_box  = mfi.tilebox();
        const RealBox tile_realbox{tile_box, geom.CellSize(), geom.ProbLo()};
        const int grid_id = mfi.index();
        const int tile_id = mfi.LocalTileIndex();
        auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];

        //Assuming tile=box for now, i.e. no tiling.
        IntVect smallEnd = tile_box.smallEnd();
        IntVect bigEnd = tile_box.bigEnd();

        if(ParallelDescriptor::MyProc() == 0 && mfi.LocalTileIndex() == 0 && proc0_enter) {

            proc0_enter = false;

            std::ifstream particleFile("particles.dat");
 //           Print() << "SPEC TOTAL: " << particleInfo[0].total << "\n";
            for(int i_spec=0; i_spec < nspecies; i_spec++) {
                for (int i_part=0; i_part<particleInfo[i_spec].total;i_part++) {
                    ParticleType p;
                    p.id()  = ParticleType::NextID();
                    //std::cout << "ID: " << p.id() << "\n";
                    p.cpu() = ParallelDescriptor::MyProc();
                    p.idata(FHD_intData::sorted) = 0;

                    if(particle_placement == 1)
                    {
                        particleFile >> p.pos(0);
                        particleFile >> p.pos(1);
                        particleFile >> p.pos(2);

                        particleFile >> p.idata(FHD_intData::pinned);

                        if(p.idata(FHD_intData::pinned) != 0)
                        {
                            pinnedParticles++;
                        }
                    }else
                    {

                        p.pos(0) = prob_lo[0] + amrex::Random()*(prob_hi[0]-prob_lo[0]);
                        p.pos(1) = prob_lo[1] + amrex::Random()*(prob_hi[1]-prob_lo[1]);
                        p.pos(2) = prob_lo[2] + amrex::Random()*(prob_hi[2]-prob_lo[2]);

                        p.idata(FHD_intData::pinned) = 0;

                    }

                    p.rdata(FHD_realData::spring) = 0;

                    p.idata(FHD_intData::visible) = 1;

                    p.rdata(FHD_realData::q) = particleInfo[i_spec].q;

 //                    std::cout << "proc " << ParallelDescriptor::MyProc() << " Pos: " << p.pos(0) << ", " << p.pos(1) << ", " << p.pos(2)
 //                              << ", " << p.rdata(FHD_realData::q) << ", " << p.id() << "\n" ;

                    //original position stored for MSD calculations
                    p.rdata(FHD_realData::ox) = p.pos(0);
                    p.rdata(FHD_realData::oy) = p.pos(1);
#if (BL_SPACEDIM == 3)
                    p.rdata(FHD_realData::oz) = p.pos(2);
#endif

                    //p.rdata(FHD_realData::vx) = sqrt(particleInfo.R*particleInfo.T)*get_particle_normal_func();
                    //p.rdata(FHD_realData::vy) = sqrt(particleInfo.R*particleInfo.T)*get_particle_normal_func();
                    //p.rdata(FHD_realData::vz) = sqrt(particleInfo.R*particleInfo.T)*get_particle_normal_func();

                    p.rdata(FHD_realData::pred_posx) = 0;
                    p.rdata(FHD_realData::pred_posy) = 0;
                    p.rdata(FHD_realData::pred_posz) = 0;

                    p.rdata(FHD_realData::pred_velx) = 0;
                    p.rdata(FHD_realData::pred_vely) = 0;
                    p.rdata(FHD_realData::pred_velz) = 0;

                    p.rdata(FHD_realData::pred_forcex) = 0;
                    p.rdata(FHD_realData::pred_forcey) = 0;
                    p.rdata(FHD_realData::pred_forcez) = 0;

                    p.rdata(FHD_realData::fx) = 0;
                    p.rdata(FHD_realData::fy) = 0;
                    p.rdata(FHD_realData::fz) = 0;

                    p.rdata(FHD_realData::vx) = 0;
                    p.rdata(FHD_realData::vy) = 0;
                    p.rdata(FHD_realData::vz) = 0;

                    p.rdata(FHD_realData::ux) = 0;
                    p.rdata(FHD_realData::uy) = 0;
                    p.rdata(FHD_realData::uz) = 0;

                    p.rdata(FHD_realData::ax) = 0;
                    p.rdata(FHD_realData::ay) = 0;
                    p.rdata(FHD_realData::az) = 0;

                    p.rdata(FHD_realData::fx) = 0;
                    p.rdata(FHD_realData::fy) = 0;
                    p.rdata(FHD_realData::fz) = 0;

                    p.rdata(FHD_realData::travelTime) = 0;
                    p.rdata(FHD_realData::diffAv) = 0;
                    p.rdata(FHD_realData::stepCount) = 0;

                    p.rdata(FHD_realData::mass) = particleInfo[i_spec].m; //mass
                    p.rdata(FHD_realData::R) = particleInfo[i_spec].R; //R
                    p.rdata(FHD_realData::radius) = particleInfo[i_spec].d/2.0; //radius
                    p.rdata(FHD_realData::accelFactor) = -6*3.14159265359*p.rdata(FHD_realData::radius)/p.rdata(FHD_realData::mass); //acceleration factor (replace with amrex c++ constant for pi...)
                    p.rdata(FHD_realData::dragFactor) = 6*3.14159265359*p.rdata(FHD_realData::radius); //drag factor
                    //p.rdata(FHD_realData::dragFactor) = 0; //drag factor
                    //p.rdata(FHD_realData::dragFactor) = 6*3.14159265359*dx[0]*1.322; //drag factor

                    p.rdata(FHD_realData::wetDiff) = particleInfo[i_spec].wetDiff;
                    p.rdata(FHD_realData::dryDiff) = particleInfo[i_spec].dryDiff;
                    p.rdata(FHD_realData::totalDiff) = particleInfo[i_spec].totalDiff;

                    p.rdata(FHD_realData::sigma) = particleInfo[i_spec].sigma;
                    p.rdata(FHD_realData::eepsilon) = particleInfo[i_spec].eepsilon;

                    p.idata(FHD_intData::species) = i_spec +1;
                    p.rdata(FHD_realData::potential) = 0;

                    // set distance for which we do direct, short range coulomb force calculation
                    // in p3m to be 6.5*dx_poisson_grid

                    p.rdata(FHD_realData::p3m_radius) = (pkernel_es[p.idata(FHD_intData::species)-1] + 0.5)*dxp[0];

                    particle_tile.push_back(p);

                    pcount++;
                }
            }

            particleFile.close();
        }
    }


    ParallelDescriptor::ReduceIntSum(pinnedParticles);

    Print() << "Loaded " << pinnedParticles << " pinned particles." << std::endl;

    if(pinnedParticles > 0)
    {
        loadPinMatrix(pinnedParticles, "matrixInv.dat");
    }

    totalPinnedMarkers = pinnedParticles;

    Redistribute();
    //clearNeighbors();
    //fillNeighbors();

}



//THIS NEEDS TO BE UPDATED TO HANDLE MORE PARTICLE INFO
void FhdParticleContainer::ReInitParticles()
{
    const int lev = 0;
    const Geometry& geom = Geom(lev);

    int pcount = 0;

    bool proc0_enter = true;
    int pinnedParticles = 0;

    //Note we are resetting the particle ID count here, this is only valid if one rank is doing the generating.
    ParticleType::NextID(1);

    for (MFIter mfi = MakeMFIter(lev, true); mfi.isValid(); ++mfi) {

        const Box& tile_box  = mfi.tilebox();
        const RealBox tile_realbox{tile_box, geom.CellSize(), geom.ProbLo()};
        const int grid_id = mfi.index();
        const int tile_id = mfi.LocalTileIndex();
        auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];
        auto& aos = particle_tile.GetArrayOfStructs();
        int np = aos.numParticles();
        auto* pstruct = aos().dataPtr();

        //Assuming tile=box for now, i.e. no tiling.
        IntVect smallEnd = tile_box.smallEnd();
        IntVect bigEnd = tile_box.bigEnd();

        if (ParallelDescriptor::MyProc() == 0 && mfi.LocalTileIndex() == 0 && proc0_enter) {

            proc0_enter = false;
            for (int i=0; i<np; ++i) {
                auto& p = pstruct[i];
                if(p.idata(FHD_intData::pinned) != 0) {
                    pinnedParticles++;
                }
            }
        }
    }

    ParallelDescriptor::ReduceIntSum(pinnedParticles);

    Print() << "Loaded " << pinnedParticles << " pinned particles." << std::endl;

    totalPinnedMarkers = pinnedParticles;

    if(pinnedParticles > 0)
    {
        loadPinMatrix(pinnedParticles, "matrixInv.dat");
    }

    Redistribute();
    doRedist = 1;

}