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
    const int nprocs = ParallelDescriptor::NProcs();

    // Inverse cell-size vector => used for determining index corresponding to
    // IBParticle position (pos)
    RealVect inv_dx = RealVect(AMREX_D_DECL(geom.InvCellSize(0),
                                            geom.InvCellSize(1),
                                            geom.InvCellSize(2)   ));
    int speciesAccumulateP[nspecies];
    int speciesPerRank[nspecies][nprocs];
    int pcount = 0;
    int totalParticles = 0;
    int pinnedParticles = 0;
    Real y_lo_wall = 0.0;
    bool proc0_enter = true;

    // calculate accumulate number of particles for species.
    //   E.g. nspecies = 3, particle_count = 4 5 6, then speciesAccumulateP = 4 9 15
    //   Initially thinking about using this method to randomly assign species to particles generated on each rank
    for (int i_spec=0; i_spec < nspecies; i_spec++) {
        totalParticles += particleInfo[i_spec].total;
        if (i_spec == 0) {
            speciesAccumulateP[i_spec] = particleInfo[i_spec].total;
        } else {
            speciesAccumulateP[i_spec] = speciesAccumulateP[i_spec-1]+particleInfo[i_spec].total;
        }
    }

    // load bond connectivity and bond strengths into two separate arrays.
    loadBonds(totalParticles, "bonds.csv");

    for (MFIter mfi = MakeMFIter(lev, true); mfi.isValid(); ++mfi) {

        pinnedParticlesIDGlobal.clear();
        const Box& tile_box  = mfi.tilebox();
        const RealBox tile_realbox{tile_box, geom.CellSize(), geom.ProbLo()};
        const int grid_id = mfi.index();
        const int tile_id = mfi.LocalTileIndex();
        auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];

        //Assuming tile=box for now, i.e. no tiling.
        IntVect smallEnd = tile_box.smallEnd();
        IntVect bigEnd = tile_box.bigEnd();

        //const Real* box_lo = tile_realbox.lo();
        //const Real* box_hi = tile_realbox.hi();


        // for now, only reading particles.dat file will generate particles in parallel
        if(particle_placement == 1)
        {
            std::ifstream particleFile("particles.dat");
 //       Print() << "SPEC TOTAL: " << particleInfo[0].total << "\n";
            for (int i_part=0; i_part<totalParticles; i_part++) {

                int buf_id, buf_species, buf_pinned, buf_groupid;
                Real buf_pos0, buf_pos1, buf_pos2;

                particleFile >> buf_id;
                particleFile >> buf_pos0;
                particleFile >> buf_pos1;
                particleFile >> buf_pos2;
                particleFile >> buf_species;
                particleFile >> buf_pinned;
                particleFile >> buf_groupid;

                if (buf_pinned == 1) {
                    pinnedParticlesIDGlobal.push_back(buf_id);
                }

                RealVect buf_pos = RealVect{buf_pos0, buf_pos1, buf_pos2};

                // IntVect representing particle's position in the tile_box grid.
                RealVect pos_grid = buf_pos;
                pos_grid *= inv_dx;
                IntVect pos_ind = IntVect(AMREX_D_DECL((int) pos_grid[0],
                                                       (int) pos_grid[1],
                                                       (int) pos_grid[2]   ));

                // Add particle at position pos iff it's vector index is contained
                // within tile_box.
                if(tile_box.contains(pos_ind)) {
                    ParticleType p;

                    p.idata(FHD_intData::sorted) = 0;
                    p.id()  = ParticleType::NextID();
                    //std::cout << "ID: " << p.id() << "\n";
                    p.cpu() = ParallelDescriptor::MyProc();
                    p.idata(FHD_intData::id_global) = buf_id;
                    p.pos(0) = buf_pos[0];
                    p.pos(1) = buf_pos[1];
                    p.pos(2) = buf_pos[2];

                    p.idata(FHD_intData::species) = buf_species;

                    p.idata(FHD_intData::pinned) = buf_pinned;

                    p.idata(FHD_intData::groupid) = buf_groupid;
                    //particleFile >> p.idata(FHD_intData::prev);
                    //particleFile >> p.idata(FHD_intData::next);
                    for(int i=0; i<MAX_BONDS; i++)
                    {
                        if (i<num_bond[p.idata(FHD_intData::id_global)]) {
                            p.rdata(FHD_realData::bondCoeff1_1 + i) = bond_coeff1[head_index[p.idata(FHD_intData::id_global)]+i];
                            p.rdata(FHD_realData::bondCoeff2_1 + i) = bond_coeff2[head_index[p.idata(FHD_intData::id_global)]+i];
                            p.idata(FHD_intData::bond1 + i) = bond_atom[head_index[p.idata(FHD_intData::id_global)]+i];
                            std::cout << "id_global=" << p.idata(FHD_intData::id_global) << "; " << "bond=" << p.idata(FHD_intData::bond1 + i) << ", k=" << p.rdata(FHD_realData::bondCoeff1_1+i) << ", x0=" << p.rdata(FHD_realData::bondCoeff2_1+i) << "\n";
                        } else {
                            p.rdata(FHD_realData::bondCoeff1_1 + i) = 0.;
                            p.rdata(FHD_realData::bondCoeff2_1 + i) = 0.;
                            p.idata(FHD_intData::bond1 + i) = -1;
                        }

                    }

                    if(p.idata(FHD_intData::pinned) != 0)
                    {
                        //pinnedParticlesIDGlobal.push_back(p.idata(FHD_intData::id_global));
                        pinnedParticles++;
                        y_lo_wall = p.pos(1);
                    }

                    p.rdata(FHD_realData::spring) = 1;

                    p.idata(FHD_intData::visible) = 1;

                    p.rdata(FHD_realData::q) = particleInfo[p.idata(FHD_intData::species)-1].q;

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

                    p.rdata(FHD_realData::mass) = particleInfo[p.idata(FHD_intData::species)-1].m; //mass
                    p.rdata(FHD_realData::R) = particleInfo[p.idata(FHD_intData::species)-1].R; //R
                    p.rdata(FHD_realData::radius) = particleInfo[p.idata(FHD_intData::species)-1].d/2.0; //radius
                    p.rdata(FHD_realData::accelFactor) = -6*3.14159265359*p.rdata(FHD_realData::radius)/p.rdata(FHD_realData::mass); //acceleration factor (replace with amrex c++ constant for pi...)
                    p.rdata(FHD_realData::dragFactor) = 6*3.14159265359*p.rdata(FHD_realData::radius); //drag factor
                    //p.rdata(FHD_realData::dragFactor) = 0; //drag factor
                    //p.rdata(FHD_realData::dragFactor) = 6*3.14159265359*dx[0]*1.322; //drag factor

                    p.rdata(FHD_realData::wetDiff) = particleInfo[p.idata(FHD_intData::species)-1].wetDiff;
                    p.rdata(FHD_realData::dryDiff) = particleInfo[p.idata(FHD_intData::species)-1].dryDiff;
                    p.rdata(FHD_realData::totalDiff) = particleInfo[p.idata(FHD_intData::species)-1].totalDiff;

                    p.rdata(FHD_realData::sigma) = particleInfo[p.idata(FHD_intData::species)-1].sigma;
                    p.rdata(FHD_realData::eepsilon) = particleInfo[p.idata(FHD_intData::species)-1].eepsilon;

                    p.rdata(FHD_realData::potential) = 0;

                    // set distance for which we do direct, short range coulomb force calculation
                    // in p3m to be 6.5*dx_poisson_grid

                    p.rdata(FHD_realData::p3m_radius) = (pkernel_es[p.idata(FHD_intData::species)-1] + 0.5)*dxp[0];

                    particle_tile.push_back(p);

                    pcount++;


                }else {
                    //// skip the next three entries because this particle is not in the current tile
                    //particleFile >> buf_species;
                    //particleFile >> buf_pinned;
                    //particleFile >> buf_groupid;
                    continue;
                }

            }
            particleFile.close();

        }else
        {
            //TODO: generate random particles in parallel. May need to turn off tiling?
            if(ParallelDescriptor::MyProc() == 0 && mfi.LocalTileIndex() == 0 && proc0_enter) {
                int i = 0;
                proc0_enter = false;
                for(int i_spec=0; i_spec < nspecies; i_spec++) {
                    for(int i_part=0; i_part<particleInfo[i_spec].total;i_part++) {
                        ParticleType p;

                        p.idata(FHD_intData::sorted) = 0;
                        p.id()  = ParticleType::NextID();
 //                       std::cout << "ID: " << p.id() << "\n";
                        p.cpu() = ParallelDescriptor::MyProc();
                        p.idata(FHD_intData::id_global) = i;
                        i++;

                        for(int i=0; i<MAX_BONDS; i++)
                        {
                            p.rdata(FHD_realData::bondCoeff1_1 + i) = 0.;
                            p.rdata(FHD_realData::bondCoeff2_1 + i) = 0.;
                            p.idata(FHD_intData::bond1 + i) = -1;

                        }

                        p.pos(0) = prob_lo[0] + amrex::Random()*(prob_hi[0]-prob_lo[0]);
                        p.pos(1) = prob_lo[1] + amrex::Random()*(prob_hi[1]-prob_lo[1]);
                        p.pos(2) = prob_lo[2] + amrex::Random()*(prob_hi[2]-prob_lo[2]);

                        p.idata(FHD_intData::species) = i_spec +1;
                        p.idata(FHD_intData::pinned) = 0;

                        p.idata(FHD_intData::groupid) = 0;
                        //p.idata(FHD_intData::prev) = -1;
                        //p.idata(FHD_intData::next) = -1;
                        p.rdata(FHD_realData::spring) = 0;

                        p.idata(FHD_intData::visible) = 1;

                        p.rdata(FHD_realData::q) = particleInfo[p.idata(FHD_intData::species)-1].q;

 //                        std::cout << "proc " << ParallelDescriptor::MyProc() << " Pos: " << p.pos(0) << ", " << p.pos(1) << ", " << p.pos(2)
 //                                  << ", " << p.rdata(FHD_realData::q) << ", " << p.id() << "\n" ;

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

                        p.rdata(FHD_realData::mass) = particleInfo[p.idata(FHD_intData::species)-1].m; //mass
                        p.rdata(FHD_realData::R) = particleInfo[p.idata(FHD_intData::species)-1].R; //R
                        p.rdata(FHD_realData::radius) = particleInfo[p.idata(FHD_intData::species)-1].d/2.0; //radius
                        p.rdata(FHD_realData::accelFactor) = -6*3.14159265359*p.rdata(FHD_realData::radius)/p.rdata(FHD_realData::mass); //acceleration factor (replace with amrex c++ constant for pi...)
                        p.rdata(FHD_realData::dragFactor) = 6*3.14159265359*p.rdata(FHD_realData::radius); //drag factor
                        //p.rdata(FHD_realData::dragFactor) = 0; //drag factor
                        //p.rdata(FHD_realData::dragFactor) = 6*3.14159265359*dx[0]*1.322; //drag factor

                        p.rdata(FHD_realData::wetDiff) = particleInfo[p.idata(FHD_intData::species)-1].wetDiff;
                        p.rdata(FHD_realData::dryDiff) = particleInfo[p.idata(FHD_intData::species)-1].dryDiff;
                        p.rdata(FHD_realData::totalDiff) = particleInfo[p.idata(FHD_intData::species)-1].totalDiff;

                        p.rdata(FHD_realData::sigma) = particleInfo[p.idata(FHD_intData::species)-1].sigma;
                        p.rdata(FHD_realData::eepsilon) = particleInfo[p.idata(FHD_intData::species)-1].eepsilon;

                        p.rdata(FHD_realData::potential) = 0;

                        // set distance for which we do direct, short range coulomb force calculation
                        // in p3m to be 6.5*dx_poisson_grid

                        p.rdata(FHD_realData::p3m_radius) = (pkernel_es[p.idata(FHD_intData::species)-1] + 0.5)*dxp[0];

                        particle_tile.push_back(p);

                        pcount++;

                    }
                }
            }
        }
    }


    ParallelDescriptor::ReduceIntSum(pcount);
    if (pcount != totalParticles) {
        Print() << "pcount\ttotalParticles: " << pcount << "    " << totalParticles << endl;
        Abort("Total number of particles mismatch; some particles missing.");
    } else {
        Print() << "Total number of generated particles: " << pcount << std::endl;
    }
    ParallelDescriptor::ReduceIntSum(pinnedParticles);
    ParallelDescriptor::ReduceRealSum(y_lo_wall);

    if (pinnedParticles != pinnedParticlesIDGlobal.size()) {
        Print() << pinnedParticles << " pinned particles loaded." << std::endl;
        Abort("Total number of pinned particles mismatch.");
    } else {
        Print() << "Loaded " << pinnedParticles << " pinned particles." << std::endl;
    }

    if(pinnedParticles > 0)
    {
        loadPinMatrix(pinnedParticles, "matrixInv.dat");
    }

    totalPinnedMarkers = pinnedParticles;
    pinned_y = y_lo_wall;

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
    Real y_lo_wall = 0.0;

    //Note we are resetting the particle ID count here, this is only valid if one rank is doing the generating.
    //ParticleType::NextID(1);

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
                    y_lo_wall = p.pos(1);
                }
            }
        }
    }

    ParallelDescriptor::ReduceIntSum(pinnedParticles);

    Print() << "Loaded " << pinnedParticles << " pinned particles." << std::endl;

    totalPinnedMarkers = pinnedParticles;
    pinned_y = y_lo_wall;

    if(pinnedParticles > 0)
    {
        loadPinMatrix(pinnedParticles, "matrixInv.dat");
    }

    Redistribute();
    doRedist = 1;

}