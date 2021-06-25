#include "INS_functions.H"
#include "common_functions.H"
#include "DsmcParticleContainer.H"
#include <sstream>
#include <string>
#include <fstream>

void FhdParticleContainer::InitParticles()
{
    const int lev = 0;
    const Geometry& geom = Geom(lev);

    int pcount = 0;

    bool proc0_enter = true;
        
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
                for (int i_part=0; i_part<properties[i_spec].total;i_part++) {
                    ParticleType p;
                    p.id()  = ParticleType::NextID();
 //                   std::cout << "ID: " << p.id() << "\n";
                    p.cpu() = ParallelDescriptor::MyProc();
                    p.idata(FHD_intData::sorted) = -1;

                    p.idata(FHD_intData::species) = i_spec;

                    p.idata(FHD_intData::i) = -100;
                    p.idata(FHD_intData::j) = -100;
                    p.idata(FHD_intData::k) = -100;

                    if(particle_placement == 1)
                    {
                        particleFile >> p.pos(0);                       
                        particleFile >> p.pos(1);
                        particleFile >> p.pos(2);

                    }else
                    {

                        p.pos(0) = prob_lo[0] + amrex::Random()*(prob_hi[0]-prob_lo[0]);
                        p.pos(1) = prob_lo[1] + amrex::Random()*(prob_hi[1]-prob_lo[1]);
                        p.pos(2) = prob_lo[2] + amrex::Random()*(prob_hi[2]-prob_lo[2]);

                    }

                    p.rdata(FHD_realData::velx) = sqrt(properties[i_spec].R*properties[i_spec].T)*amrex::RandomNormal(0.,1.);
                    p.rdata(FHD_realData::vely) = sqrt(properties[i_spec].R*properties[i_spec].T)*amrex::RandomNormal(0.,1.);
                    p.rdata(FHD_realData::velz) = sqrt(properties[i_spec].R*properties[i_spec].T)*amrex::RandomNormal(0.,1.);

                    p.rdata(FHD_realData::boostx) = 0;
                    p.rdata(FHD_realData::boosty) = 0;
                    p.rdata(FHD_realData::boostz) = 0;

                    p.rdata(FHD_realData::timeFrac) = 1;

                    p.rdata(FHD_realData::R) = properties[i_spec].R;

                    particle_tile.push_back(p);

                    pcount++;
                }
            }

            particleFile.close();
        }
    }

    Redistribute();
    SortParticles();
}



