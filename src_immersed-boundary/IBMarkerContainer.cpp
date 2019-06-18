#include <AMReX.H>
#include <AMReX_AmrParGDB.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Particles.H>
#include <AMReX_NeighborParticles.H>

#include <AMReX_VisMF.H>  // amrex::VisMF::Write(MultiFab)

#include <common_functions.H>
#include <common_namespace.H>

#include <IBMarkerContainer.H>
#include <ib_functions_F.H>



using namespace common;
using namespace amrex;

bool IBMarkerContainer::use_neighbor_list  {true};
bool IBMarkerContainer::sort_neighbor_list {false};



IBMarkerContainer::IBMarkerContainer(const Geometry & geom,
                                         const DistributionMapping & dmap,
                                         const BoxArray & ba,
                                         int n_nbhd)
    : NeighborParticleContainer<IBM_realData::count, IBM_intData::count>(
            geom, dmap, ba, n_nbhd
        ),
    nghost(n_nbhd)
{
    InitInternals(n_nbhd);
}



IBMarkerContainer::IBMarkerContainer(AmrCore * amr_core, int n_nbhd)
    : NeighborParticleContainer<IBM_realData::count, IBM_intData::count>(
            amr_core->GetParGDB(), n_nbhd
        ),
    m_amr_core(amr_core),
    nghost(n_nbhd)
{
    InitInternals(n_nbhd);
}



void IBMarkerContainer::InitList(int lev, const Vector<RealVect> & pos) {

    // Inverse cell-size vector => used for determining index corresponding to
    // IBParticle position (pos)
    RealVect inv_dx = RealVect(AMREX_D_DECL(Geom(lev).InvCellSize(0),
                                            Geom(lev).InvCellSize(1),
                                            Geom(lev).InvCellSize(2)   ));


    int total_np = 0;

    // This uses the particle tile size. Note that the default is to tile so if
    // we remove the true and don't explicitly add false it will still tile.
    for (MFIter mfi = MakeMFIter(lev, true); mfi.isValid(); ++mfi) {

        // This is particles per grid so we reset to 0
        int pcount = 0;

        // Current tile box
        const Box & tile_box = mfi.tilebox();

        // Create a particle container for this grid and add the
        // immersed-boundary particles to it if the particle's position (pos[i])
        // is within current tile box.
        const int grid_id = mfi.index();
        const int tile_id = mfi.LocalTileIndex();
        auto & particles = GetParticles(lev)[std::make_pair(grid_id,tile_id)];

        for(int i = 0; i < pos.size(); i++) {
            // IntVect representing particle's position in the tile_box grid.
            RealVect pos_grid = pos[i];
            pos_grid *= inv_dx;
            IntVect pos_ind = IntVect(AMREX_D_DECL((int) pos_grid[0],
                                                   (int) pos_grid[1],
                                                   (int) pos_grid[2]   ));

            // Add particle at position pos iff it's vector index is contained
            // within tile_box.
            if(tile_box.contains(pos_ind)) {
                pcount ++;

                ParticleType p_new;

                // Set id and cpu for this particle
                p_new.id()  = ParticleType::NextID();
                p_new.cpu() = ParallelDescriptor::MyProc();

                // Set particle position
                p_new.pos(0) = pos[i][0];
                p_new.pos(1) = pos[i][1];
                p_new.pos(2) = pos[i][2];

                // Initialize marker velocity as well as forces to 0
                p_new.rdata(IBM_realData::velx)   = 0.;
                p_new.rdata(IBM_realData::vely)   = 0.;
                p_new.rdata(IBM_realData::velz)   = 0.;

                p_new.rdata(IBM_realData::forcex) = 0.;
                p_new.rdata(IBM_realData::forcey) = 0.;
                p_new.rdata(IBM_realData::forcez) = 0.;

                p_new.rdata(IBM_realData::pred_forcex) = 0.;
                p_new.rdata(IBM_realData::pred_forcey) = 0.;
                p_new.rdata(IBM_realData::pred_forcez) = 0.;

                p_new.idata(IBM_intData::id_0)    = 0;
                p_new.idata(IBM_intData::cpu_0)   = 0;

                p_new.idata(IBM_intData::id_1)    = 0;
                p_new.idata(IBM_intData::cpu_1)   = 0;

                // Add to the data structure
                particles.push_back(p_new);
            }
        }

        const int np = pcount;
        total_np += np;
    }

    ParallelDescriptor::ReduceIntSum(total_np,ParallelDescriptor::IOProcessorNumber());
    std::cout << "Total number of generated particles: " << total_np << std::endl;

    // We shouldn't need this if the particles are tiled with one tile per
    // grid, but otherwise we do need this to move particles from tile 0 to the
    // correct tile.
    Redistribute();
}



void IBMarkerContainer::MoveMarkers(int lev, Real dt) {

    for (IBMarIter pti(* this, lev); pti.isValid(); ++pti) {

        PairIndex index(pti.index(), pti.LocalTileIndex());
        auto & particle_data = GetParticles(lev).at(index);
        long np = particle_data.size();

        AoS & particles = particle_data.GetArrayOfStructs();
        for (int i = 0; i < np; ++i) {
            ParticleType & part = particles[i];

            part.pos(0) += dt * part.rdata(IBM_realData::velx);
            part.pos(1) += dt * part.rdata(IBM_realData::vely);
            part.pos(2) += dt * part.rdata(IBM_realData::velz);
        }
    }
}



void IBMarkerContainer::MovePredictor(int lev, Real dt) {

    for (IBMarIter pti(* this, lev); pti.isValid(); ++pti) {

        PairIndex index(pti.index(), pti.LocalTileIndex());
        auto & particle_data = GetParticles(lev).at(index);
        long np = particle_data.size();

        AoS & particles = particle_data.GetArrayOfStructs();
        for (int i = 0; i < np; ++i) {
            ParticleType & part = particles[i];
            MarkerIndex pindex(part.id(), part.cpu());

            part.rdata(IBM_realData::pred_posx) += dt * part.rdata(IBM_realData::pred_velx);
            part.rdata(IBM_realData::pred_posy) += dt * part.rdata(IBM_realData::pred_vely);
            part.rdata(IBM_realData::pred_posz) += dt * part.rdata(IBM_realData::pred_velz);
        }
    }
}



void IBMarkerContainer::SpreadMarkers(int lev,
                                      const Vector<RealVect> & f_in,
                                      const Vector<RealVect> & f_pos,
                                      std::array<MultiFab, AMREX_SPACEDIM> & f_out,
                                      std::array<MultiFab, AMREX_SPACEDIM> & f_weights) const {


    //___________________________________________________________________________
    // Geometry data
    const Geometry & geom = Geom(0);
    const Real     *   dx = geom.CellSize();

    const int n_marker = f_pos.size();


    //___________________________________________________________________________
    // Cell-centered MultiFab used as a reference for iterating over data
    // WARNING: this will break if IBMarkerContainer is on a different grid
    // than the grid which we're spreading to

    const BoxArray & ba            = ParticleBoxArray(lev);
    const DistributionMapping & dm = ParticleDistributionMap(lev);

    MultiFab dummy(ba, dm, 1, f_out[0].nGrow());


    for (MFIter mfi(dummy); mfi.isValid(); ++mfi) {

        Box bx = mfi.growntilebox();

        spread_markers(BL_TO_FORTRAN_BOX(bx),
                       BL_TO_FORTRAN_ANYD(f_out[0][mfi]),
#if (AMREX_SPACEDIM > 1)
                       BL_TO_FORTRAN_ANYD(f_out[1][mfi]),
#endif
#if (AMREX_SPACEDIM > 2)
                       BL_TO_FORTRAN_ANYD(f_out[2][mfi]),
#endif
                       BL_TO_FORTRAN_ANYD(f_weights[0][mfi]),
#if (AMREX_SPACEDIM > 1)
                       BL_TO_FORTRAN_ANYD(f_weights[1][mfi]),
#endif
#if (AMREX_SPACEDIM > 2)
                       BL_TO_FORTRAN_ANYD(f_weights[2][mfi]),
#endif
                       BL_TO_FORTRAN_ANYD(face_coords[lev][0][mfi]),
#if (AMREX_SPACEDIM > 1)
                       BL_TO_FORTRAN_ANYD(face_coords[lev][1][mfi]),
#endif
#if (AMREX_SPACEDIM > 2)
                       BL_TO_FORTRAN_ANYD(face_coords[lev][2][mfi]),
#endif
                       f_pos.dataPtr(),
                       f_in.dataPtr(),
                       & n_marker,
                       dx );
    }
}



void IBMarkerContainer::SpreadMarkers(int lev,
                                      const Vector<RealVect> & f_in,
                                      std::array<MultiFab, AMREX_SPACEDIM> & f_out,
                                      std::array<MultiFab, AMREX_SPACEDIM> & f_weights) const {

    //___________________________________________________________________________
    // Fill vector of marker positions (for current level)
    Vector<IBM_info> marker_info = IBMarkerInfo(lev);
    Vector<RealVect> marker_positions(marker_info.size());
    for (int i=0; i<marker_info.size(); ++i)
        marker_positions[i] = marker_info[i].pos;


    SpreadMarkers(lev, f_in, marker_positions, f_out, f_weights);

}

void IBMarkerContainer::SpreadMarkers(int lev,
                                      const Vector<RealVect> & f_in,
                                      std::array<MultiFab, AMREX_SPACEDIM> & f_out) const {

    //___________________________________________________________________________
    // We don't need these spreading weights => create a dummy MF
    std::array<MultiFab, AMREX_SPACEDIM> f_weights;
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        f_weights[d].define(f_out[d].boxArray(), f_out[d].DistributionMap(),
                            1, f_out[d].nGrow());
        f_weights[d].setVal(0.);
    }

    SpreadMarkers(lev, f_in, f_out, f_weights);

}



void IBMarkerContainer::SpreadMarkers(int lev,
                                      std::array<MultiFab, AMREX_SPACEDIM> & f_out) const {


    //___________________________________________________________________________
    // We don't need these spreading weights => create a dummy MF
    std::array<MultiFab, AMREX_SPACEDIM> f_weights;
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        f_weights[d].define(f_out[d].boxArray(), f_out[d].DistributionMap(),
                            1, f_out[d].nGrow());
        f_weights[d].setVal(0.);
    }


    for (MFIter pti = MakeMFIter(lev, true); pti.isValid(); ++pti) {

        // Marker (non-neighbor particle) data for current tile
        PairIndex index(pti.index(), pti.LocalTileIndex());
        const auto & particle_data = GetParticles(lev).at(index);


        //_______________________________________________________________________
        // Fill vector of marker positions and forces (for current level)
        long np = particle_data.size();

        Vector<RealVect> marker_positions(np), marker_forces(np);

        const AoS & particles = particle_data.GetArrayOfStructs();
        for (int i = 0; i < np; ++i) {
            const ParticleType & part = particles[i];

            RealVect ppos, pfor;
            for (int d=0; d<AMREX_SPACEDIM; ++d)
                ppos[d] = part.pos(d);

            pfor[0] = part.rdata(IBM_realData::forcex);
#if (AMREX_SPACEDIM > 1)
            pfor[1] = part.rdata(IBM_realData::forcey);
#endif
#if (AMREX_SPACEDIM > 2)
            pfor[2] = part.rdata(IBM_realData::forcez);
#endif

            marker_positions[i] = ppos;
            marker_forces[i]    = pfor;
        }

        //_______________________________________________________________________
        // Spread the non-neighbor particles (markers)
        SpreadMarkers(lev, marker_forces, marker_positions, f_out, f_weights);

        // Clear vectors => to be filled with neighbor data now
        marker_positions.clear();
        marker_forces.clear();


        // Neighbor marker data for current tile
        const ParticleType * nbhd_data = (ParticleType *) neighbors[lev].at(index).dataPtr();

        //_______________________________________________________________________
        // Fill vector of neighbor marker positions and forces
        int ng = neighbors[lev].at(index).size();

        marker_positions.resize(ng);
        marker_forces.resize(ng);

        for (int i = 0; i < ng; ++i) {
            const ParticleType & part = nbhd_data[i];

            RealVect ppos, pfor;
            for (int d=0; d<AMREX_SPACEDIM; ++d)
                ppos[d] = part.pos(d);

            pfor[0] = part.rdata(IBM_realData::forcex);
#if (AMREX_SPACEDIM > 1)
            pfor[1] = part.rdata(IBM_realData::forcey);
#endif
#if (AMREX_SPACEDIM > 2)
            pfor[2] = part.rdata(IBM_realData::forcez);
#endif

            marker_positions[i] = ppos;
            marker_forces[i]    = pfor;
        }

        //_______________________________________________________________________
        // Spread the neighbor particles (markers)
        SpreadMarkers(lev, marker_forces, marker_positions, f_out, f_weights);
    }
}



void IBMarkerContainer::SpreadPredictor(int lev,
                                        std::array<MultiFab, AMREX_SPACEDIM> & f_out) const {


    //___________________________________________________________________________
    // We don't need these spreading weights => create a dummy MF
    std::array<MultiFab, AMREX_SPACEDIM> f_weights;
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        f_weights[d].define(f_out[d].boxArray(), f_out[d].DistributionMap(),
                            1, f_out[d].nGrow());
        f_weights[d].setVal(0.);
    }


    for (MFIter pti = MakeMFIter(lev, true); pti.isValid(); ++pti) {

        // Marker (non-neighbor particle) data for current tile
        PairIndex index(pti.index(), pti.LocalTileIndex());
        const auto & particle_data = GetParticles(lev).at(index);

        //_______________________________________________________________________
        // Fill vector of marker positions and forces (for current level)
        long np = particle_data.size();

        Vector<RealVect> marker_positions(np), marker_forces(np);

        const AoS & particles = particle_data.GetArrayOfStructs();
        for (int i = 0; i < np; ++i) {
            const ParticleType & part = particles[i];

            RealVect ppos, pfor;

            ppos[0] = part.rdata(IBM_realData::pred_posx);
#if (AMREX_SPACEDIM > 1)
            ppos[1] = part.rdata(IBM_realData::pred_posy);
#endif
#if (AMREX_SPACEDIM > 2)
            ppos[2] = part.rdata(IBM_realData::pred_posz);
#endif

            pfor[0] = part.rdata(IBM_realData::pred_forcex);
#if (AMREX_SPACEDIM > 1)
            pfor[1] = part.rdata(IBM_realData::pred_forcey);
#endif
#if (AMREX_SPACEDIM > 2)
            pfor[2] = part.rdata(IBM_realData::pred_forcez);
#endif

            marker_positions[i] = ppos;
            marker_forces[i]    = pfor;
        }

        //_______________________________________________________________________
        // Spread the non-neighbor particles (markers)
        SpreadMarkers(lev, marker_forces, marker_positions, f_out, f_weights);

        // Clear vectors => to be filled with neighbor data now
        marker_positions.clear();
        marker_forces.clear();


        // Neighbor marker data for current tile
        const ParticleType * nbhd_data = (ParticleType *) neighbors[lev].at(index).dataPtr();


        //_______________________________________________________________________
        // Fill vector of neighbor marker positions and forces
        int ng = neighbors[lev].at(index).size();

        marker_positions.resize(ng);
        marker_forces.resize(ng);

        for (int i = 0; i < ng; ++i) {
            const ParticleType & part = nbhd_data[i];

            RealVect ppos, pfor;

            ppos[0] = part.rdata(IBM_realData::pred_posx);
#if (AMREX_SPACEDIM > 1)
            ppos[1] = part.rdata(IBM_realData::pred_posy);
#endif
#if (AMREX_SPACEDIM > 2)
            ppos[2] = part.rdata(IBM_realData::pred_posz);
#endif

            pfor[0] = part.rdata(IBM_realData::pred_forcex);
#if (AMREX_SPACEDIM > 1)
            pfor[1] = part.rdata(IBM_realData::pred_forcey);
#endif
#if (AMREX_SPACEDIM > 2)
            pfor[2] = part.rdata(IBM_realData::pred_forcez);
#endif

            marker_positions[i] = ppos;
            marker_forces[i]    = pfor;
        }

        //_______________________________________________________________________
        // Spread the eighbor particles (markers)
        SpreadMarkers(lev, marker_forces, marker_positions, f_out, f_weights);
    }
}



void IBMarkerContainer::InterpolateMarkers(int lev,
                                           Vector<RealVect> & f_out,
                                           const Vector<RealVect> & f_pos,
                                           const std::array<MultiFab, AMREX_SPACEDIM> & f_in,
                                           const std::array<MultiFab, AMREX_SPACEDIM> & f_weights) const {


    //___________________________________________________________________________
    // Geometry data
    const Geometry & geom = Geom(0);
    const Real     *   dx = geom.CellSize();

    const int n_marker = f_pos.size();


    //___________________________________________________________________________
    // Cell-centered MultiFab used as a reference for iterating over data
    // WARNING: this will break if IBMarkerContainer is on a differnt grid
    // than the grid which we're interpolating from

    const BoxArray & ba            = ParticleBoxArray(lev);
    const DistributionMapping & dm = ParticleDistributionMap(lev);

    MultiFab dummy(ba, dm, 1, f_in[0].nGrow());


    for (MFIter mfi(dummy); mfi.isValid(); ++mfi) {

        Box bx = mfi.growntilebox();

        interpolate_markers(BL_TO_FORTRAN_BOX(bx),
                            BL_TO_FORTRAN_ANYD(f_in[0][mfi]),
#if (AMREX_SPACEDIM > 1)
                            BL_TO_FORTRAN_ANYD(f_in[1][mfi]),
#endif
#if (AMREX_SPACEDIM > 2)
                            BL_TO_FORTRAN_ANYD(f_in[2][mfi]),
#endif
                            BL_TO_FORTRAN_ANYD(f_weights[0][mfi]),
#if (AMREX_SPACEDIM > 1)
                            BL_TO_FORTRAN_ANYD(f_weights[1][mfi]),
#endif
#if (AMREX_SPACEDIM > 2)
                            BL_TO_FORTRAN_ANYD(f_weights[2][mfi]),
#endif
                            BL_TO_FORTRAN_ANYD(face_coords[lev][0][mfi]),
#if (AMREX_SPACEDIM > 1)
                            BL_TO_FORTRAN_ANYD(face_coords[lev][1][mfi]),
#endif
#if (AMREX_SPACEDIM > 2)
                            BL_TO_FORTRAN_ANYD(face_coords[lev][2][mfi]),
#endif
                            f_pos.dataPtr(),
                            f_out.dataPtr(),
                            & n_marker,
                            dx );
    }

}




void IBMarkerContainer::InterpolateMarkers(int lev,
                                           Vector<RealVect> & f_out,
                                           const std::array<MultiFab, AMREX_SPACEDIM> & f_in,
                                           const std::array<MultiFab, AMREX_SPACEDIM> & f_weights) const {


    //___________________________________________________________________________
    // Fill vector of marker positions (for current level)
    Vector<IBM_info> marker_info = IBMarkerInfo(lev);
    Vector<RealVect> marker_positions(marker_info.size());
    for (int i=0; i<marker_info.size(); ++i)
        marker_positions[i] = marker_info[i].pos;

    InterpolateMarkers(lev, f_out, marker_positions, f_in, f_weights);
}



void IBMarkerContainer::InterpolateMarkers(int lev,
        Vector<RealVect> & f_out, const std::array<MultiFab, AMREX_SPACEDIM> & f_in) const {

    //___________________________________________________________________________
    // We don't need these spreading weights => create a dummy MF
    std::array<MultiFab, AMREX_SPACEDIM> f_weights;
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        f_weights[d].define(f_in[d].boxArray(), f_in[d].DistributionMap(),
                            1, f_in[d].nGrow());
        f_weights[d].setVal(-1.); // Set to <0 to guarantee that weights are ignored
    }

    InterpolateMarkers(lev, f_out, f_in, f_weights);

}



void IBMarkerContainer::InterpolateMarkers(int lev,
                                           const std::array<MultiFab, AMREX_SPACEDIM> & f_out) {


    //___________________________________________________________________________
    // We don't need these spreading weights => create a dummy MF
    std::array<MultiFab, AMREX_SPACEDIM> f_weights;
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        f_weights[d].define(f_out[d].boxArray(), f_out[d].DistributionMap(),
                            1, f_out[d].nGrow());
        f_weights[d].setVal(-1.);
    }


    for (MFIter pti = MakeMFIter(lev, true); pti.isValid(); ++pti) {

        // Marker (non-neighbor particle) data for current tile
        PairIndex index(pti.index(), pti.LocalTileIndex());
        auto & particle_data = GetParticles(lev).at(index);

        //_______________________________________________________________________
        // Fill vector of marker positions and forces (for current level)
        long np = particle_data.size();

        Vector<RealVect> marker_positions(np), marker_forces(np);

        AoS & particles = particle_data.GetArrayOfStructs();
        for (int i = 0; i < np; ++i) {
            const ParticleType & part = particles[i];

            RealVect ppos, pfor;
            for (int d=0; d<AMREX_SPACEDIM; ++d)
                ppos[d] = part.pos(d);

            // pfor should be (0, .., 0)

            marker_positions[i] = ppos;
            marker_forces[i]    = pfor;
        }


        //_______________________________________________________________________
        // Interpolate the non-neighbor particles (markers)
        InterpolateMarkers(lev, marker_forces, marker_positions, f_out, f_weights);

        // Add interpolated markers back to the particles (markers)
        for (int i = 0; i < np; ++i) {
            ParticleType & part = particles[i];

            part.rdata(IBM_realData::velx) += marker_forces[i][0];
#if (AMREX_SPACEDIM > 1)
            part.rdata(IBM_realData::vely) += marker_forces[i][1];
#endif
#if (AMREX_SPACEDIM > 2)
            part.rdata(IBM_realData::velz) += marker_forces[i][2];
#endif
        }

        // Clear vectors => to be filled with neighbor data now
        marker_positions.clear();
        marker_forces.clear();


        // Neighbor marker data for current tile
        ParticleType * nbhd_data = (ParticleType *) neighbors[lev].at(index).dataPtr();


        //_______________________________________________________________________
        // Fill vector of neighbor marker positions and forces
        int ng = neighbors[lev].at(index).size();

        marker_positions.resize(ng);
        marker_forces.resize(ng);

        for (int i = 0; i < ng; ++i) {
            ParticleType & part = nbhd_data[i];

            RealVect ppos, pfor;
            for (int d=0; d<AMREX_SPACEDIM; ++d)
                ppos[d] = part.pos(d);

            // pfor should be (0, .., 0)

            marker_positions[i] = ppos;
            marker_forces[i]    = pfor;
        }

        //_______________________________________________________________________
        // Interpolate the neighbor particles (markers)
        InterpolateMarkers(lev, marker_forces, marker_positions, f_out, f_weights);

        // Add interpolated markers back to the particles (markers)
        for (int i = 0; i < ng; ++i) {
            ParticleType & part = nbhd_data[i];

            part.rdata(IBM_realData::velx) += marker_forces[i][0];
#if (AMREX_SPACEDIM > 1)
            part.rdata(IBM_realData::vely) += marker_forces[i][1];
#endif
#if (AMREX_SPACEDIM > 2)
            part.rdata(IBM_realData::velz) += marker_forces[i][2];
#endif

        }

        // TODO: sync neighbors?
    }
}




void IBMarkerContainer::InterpolatePredictor(int lev,
                                             const std::array<MultiFab, AMREX_SPACEDIM> & f_in) {

    //___________________________________________________________________________
    // Fill vector of marker positions and forces (for current level)
    Vector<RealVect> marker_positions, marker_forces;

    for (MFIter pti = MakeMFIter(lev, true); pti.isValid(); ++pti) {

        PairIndex index(pti.index(), pti.LocalTileIndex());
        const auto & particle_data = GetParticles(lev).at(index);
        long np = particle_data.size();

        const AoS & particles = particle_data.GetArrayOfStructs();
        for (int i = 0; i < np; ++i) {
            const ParticleType & part = particles[i];

            RealVect ppos, pfor;

            ppos[0] = part.rdata(IBM_realData::pred_posx);
#if (AMREX_SPACEDIM > 1)
            ppos[1] = part.rdata(IBM_realData::pred_posy);
#endif
#if (AMREX_SPACEDIM > 2)
            ppos[2] = part.rdata(IBM_realData::pred_posz);
#endif

            // pfor should be (0, .., 0)
#ifdef _OPENMP
#pragma omp critical
#endif
            {
                marker_positions.push_back(ppos);
                marker_forces.push_back(pfor);
            }
        }

        int ng = neighbors[lev].at(index).size();
        const ParticleType * nbhd_data = (ParticleType *) neighbors[lev].at(index).dataPtr();
        for (int i = 0; i < np; ++i) {
            const ParticleType & part = nbhd_data[i];

            RealVect ppos, pfor;

            ppos[0] = part.rdata(IBM_realData::pred_posx);
#if (AMREX_SPACEDIM > 1)
            ppos[1] = part.rdata(IBM_realData::pred_posy);
#endif
#if (AMREX_SPACEDIM > 2)
            ppos[2] = part.rdata(IBM_realData::pred_posz);
#endif


            // pfor should be (0, .., 0)

#ifdef _OPENMP
#pragma omp critical
#endif
            {
                marker_positions.push_back(ppos);
                marker_forces.push_back(pfor);
            }
        }
    }


    //___________________________________________________________________________
    // We don't need these spreading weights => create a dummy MF
    std::array<MultiFab, AMREX_SPACEDIM> f_weights;
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        f_weights[d].define(f_in[d].boxArray(), f_in[d].DistributionMap(),
                            1, f_in[d].nGrow());
        f_weights[d].setVal(-1.);
    }


    //___________________________________________________________________________
    // Geometry data
    const Geometry & geom = Geom(0);
    const Real     *   dx = geom.CellSize();

    const int n_marker = marker_positions.size();


    //___________________________________________________________________________
    // Cell-centered MultiFab used as a reference for iterating over data
    // WARNING: this will break if IBMarkerContainer is on a differnt grid
    // than the grid which we're spreading to

    const BoxArray & ba            = ParticleBoxArray(lev);
    const DistributionMapping & dm = ParticleDistributionMap(lev);

    MultiFab dummy(ba, dm, 1, f_in[0].nGrow());


    for (MFIter mfi(dummy); mfi.isValid(); ++mfi) {

        Box bx = mfi.growntilebox();

        interpolate_markers(BL_TO_FORTRAN_BOX(bx),
                            BL_TO_FORTRAN_ANYD(f_in[0][mfi]),
#if (AMREX_SPACEDIM > 1)
                            BL_TO_FORTRAN_ANYD(f_in[1][mfi]),
#endif
#if (AMREX_SPACEDIM > 2)
                            BL_TO_FORTRAN_ANYD(f_in[2][mfi]),
#endif
                            BL_TO_FORTRAN_ANYD(f_weights[0][mfi]),
#if (AMREX_SPACEDIM > 1)
                            BL_TO_FORTRAN_ANYD(f_weights[1][mfi]),
#endif
#if (AMREX_SPACEDIM > 2)
                            BL_TO_FORTRAN_ANYD(f_weights[2][mfi]),
#endif
                            BL_TO_FORTRAN_ANYD(face_coords[lev][0][mfi]),
#if (AMREX_SPACEDIM > 1)
                            BL_TO_FORTRAN_ANYD(face_coords[lev][1][mfi]),
#endif
#if (AMREX_SPACEDIM > 2)
                            BL_TO_FORTRAN_ANYD(face_coords[lev][2][mfi]),
#endif
                            marker_positions.dataPtr(),
                            marker_forces.dataPtr(),
                            & n_marker,
                            dx );
    }


    //___________________________________________________________________________
    // Add interpolated markers back to the particles (markers)
    for (MFIter pti = MakeMFIter(lev, true); pti.isValid(); ++pti) {

        PairIndex index(pti.index(), pti.LocalTileIndex());
        auto & particle_data = GetParticles(lev).at(index);
        long np = particle_data.size();

        AoS & particles = particle_data.GetArrayOfStructs();
        for (int i = 0; i < np; ++i) {
            ParticleType & part = particles[i];

            part.rdata(IBM_realData::pred_velx) += marker_forces[i][0];
#if (AMREX_SPACEDIM > 1)
            part.rdata(IBM_realData::pred_vely) += marker_forces[i][1];
#endif
#if (AMREX_SPACEDIM > 2)
            part.rdata(IBM_realData::pred_velz) += marker_forces[i][2];
#endif
        }

        int ng = neighbors[lev].at(index).size();
        ParticleType * nbhd_data = (ParticleType *) neighbors[lev].at(index).dataPtr();
        for (int i = 0; i < np; ++i) {
            ParticleType & part = nbhd_data[i];

            part.rdata(IBM_realData::pred_velx) += marker_forces[i][0];
#if (AMREX_SPACEDIM > 1)
            part.rdata(IBM_realData::pred_vely) += marker_forces[i][1];
#endif
#if (AMREX_SPACEDIM > 2)
            part.rdata(IBM_realData::pred_velz) += marker_forces[i][2];
#endif

        }

        // TODO: sync neighbors?
    }



}


void IBMarkerContainer::PrintMarkerData(int lev) const {

    // Inverse cell-size vector => max is used for determining IBParticle
    // radius in units of cell size
    Vector<Real> inv_dx = {AMREX_D_DECL(Geom(lev).InvCellSize(0),
                                        Geom(lev).InvCellSize(1),
                                        Geom(lev).InvCellSize(2)   )};

    // Find max inv_dx (in case we have an anisotropic grid)
    Real mx_inv_dx = * std::max_element(inv_dx.begin(), inv_dx.end());


    amrex::AllPrintToFile("ib_marker_data") << "Particles on each box:" << std::endl;


    long local_count = 0;

    // ParIter skips tiles without particles => Iterate over MultiFab instead of
    // ParticleIter. Note also that AmrexParticleContainer uses wired tiling =>
    // turn tiling off
    for(MFIter pti = MakeMFIter(lev, false); pti.isValid(); ++pti) {
        // MuliFabs are indexed using a pair: (BoxArray index, tile index):
        PairIndex index(pti.index(), pti.LocalTileIndex());

        // Neighbours are stored as raw data (see below)
        int ng = neighbors[lev].at(index).size();

        auto & particle_data = GetParticles(lev).at(index);
        long np = particle_data.size();

        local_count += np;

        // Print current box info
        AllPrintToFile("ib_marker_data") << "Box:"         << pti.index()
                                         << " "            << pti.tilebox()
                                         << ", count: "    << np
                                         << ", nb count: " << ng
                                         << std::endl;

        // Print IBMarker
        AllPrintToFile("ib_marker_data") << " * IBMarkers:" << std::endl;

        const AoS & particles = particle_data.GetArrayOfStructs();
        for(int i = 0; i < np; i++){
            const ParticleType & part = particles[i];

            AllPrintToFile("ib_marker_data") << "   +--> " << part << std::endl;
        }

        // Print neighbour IBParticles
        AllPrintToFile("ib_marker_data") << " * Grown IBParticles:" << std::endl;

        // TODO: HAXOR!!! This should be fixed ASAP: if I understand this
        // correctly, the neighbor data contains the particle data as a binary
        // array (char). By casting to ParticleType, what we're doing is
        // interpreting the data in neighbours[index] as valid particle data.
        // Also we stride the neighbors[index] array in units of
        // sizeof(ParticleData). All of this is a little too dangerous for my
        // taste: never hide what you're doing from your compiler!!!
        const ParticleType * nbhd_data = (ParticleType *) neighbors[lev].at(index).dataPtr();
        for(int i = 0; i < ng; i++){
            const ParticleType & part = nbhd_data[i];

            AllPrintToFile("ib_marker_data") << "   +--> " << part << std::endl;
        }
    }

    AllPrintToFile("ib_marker_data") << "Total for this process: "
                                       << local_count << std::endl << std::endl;
}



void IBMarkerContainer::LocalIBMarkerInfo(Vector<IBM_info> & info,
                                          int lev, PairIndex index,
                                          bool unique) const {

    // Inverse cell-size vector => used for determining index corresponding to
    // IBParticle position (pos)
    RealVect inv_dx = RealVect(
            AMREX_D_DECL(
                Geom(lev).InvCellSize(0),
                Geom(lev).InvCellSize(1),
                Geom(lev).InvCellSize(2)
            )
        );


    auto & particle_data = GetParticles(lev).at(index);
    long np = particle_data.size();

    // Iterate over local particle data
    const AoS & particles = particle_data.GetArrayOfStructs();
    for(int i = 0; i < np; i++){
        const ParticleType & part = particles[i];

        // Position of IBParticle
        RealVect pos = RealVect(
                AMREX_D_DECL(part.pos(0), part.pos(1), part.pos(2))
            );

        // Velocity of IBParticle
        RealVect vel = RealVect(
                AMREX_D_DECL(part.rdata(IBM_realData::velx),
                             part.rdata(IBM_realData::vely),
                             part.rdata(IBM_realData::velz)   )
            );

        // Velocity of IBParticle
        RealVect force = RealVect(
                AMREX_D_DECL(part.rdata(IBM_realData::forcex),
                             part.rdata(IBM_realData::forcey),
                             part.rdata(IBM_realData::forcez)   )
            );

        // Construct info struct
        IBM_info part_info;
        part_info.pos    = pos;
        part_info.vel    = vel;
        part_info.force  = force;
        part_info.id     = part.id();
        part_info.cpu    = part.cpu();
        part_info.real   = 1; // 1 => real (non-neighbor particle)

        // Add to list

        if (unique) {
            // If in unique-mode: Don't add unless `part_info` is not already in `info`
            const auto & search = std::find(std::begin(info), std::end(info), part_info);
            if ( search == std::end(info)) {
#ifdef _OPENMP
#pragma omp critical
#endif
                { info.push_back(part_info); }
            }
        } else {
#ifdef _OPENMP
#pragma omp critical
#endif
            { info.push_back(part_info); }
        }
    }
}



Vector<IBM_info> IBMarkerContainer::LocalIBMarkerInfo(int lev, PairIndex index,
                                                      bool unique) const {

    // Allocate Marker Info vector
    Vector<IBM_info> info;

    //___________________________________________________________________________
    // Fill Marker Info vector with local (non-neighbour) data
    LocalIBMarkerInfo(info, lev, index, unique);


    return info;
}



Vector<IBM_info> IBMarkerContainer::LocalIBMarkerInfo(int lev, bool unique) const {

    // Allocate Marker Info vector
    Vector<IBM_info> info;

    //___________________________________________________________________________
    // Cell-centered MultiFab used as a reference for iterating over data
    // WARNING: this will break if IBMarkerContainer is on a differnt grid
    // than the grid which we're searching for particles (this should usually
    // be fine though)

    const BoxArray & ba            = ParticleBoxArray(lev);
    const DistributionMapping & dm = ParticleDistributionMap(lev);

    MultiFab dummy(ba, dm, 1, 1);

    //___________________________________________________________________________
    // Iterate over `dummy` looking for particles. NOTE: use the
    // IBMarkerContainer tile size
    for (MFIter mfi(dummy, tile_size); mfi.isValid(); ++mfi){
        PairIndex index(mfi.index(), mfi.LocalTileIndex());
        LocalIBMarkerInfo(info, lev, index, true);
    }


    return info;
}



void IBMarkerContainer::NeighborIBMarkerInfo(Vector<IBM_info> & info,
                                             int lev, PairIndex index,
                                             bool unique) const {

    RealVect inv_dx = RealVect(
            AMREX_D_DECL(
                Geom(lev).InvCellSize(0),
                Geom(lev).InvCellSize(1),
                Geom(lev).InvCellSize(2)
            )
        );

    int ng = neighbors[lev].at(index).size();

    // Iterate over neighbour particles:
    // TODO: HAXOR!!! This should be fixed ASAP: if I understand this correctly,
    // the neighbor data contains the particle data as a binary array (char). By
    // casting to ParticleType, what we're doing is interpreting the data in
    // neighbours[index] as valid particle data. Also we stride the
    // neighbors[index] array in units of sizeof(ParticleData). All of this is a
    // little too dangerous for my taste: never hide what you're doing from your
    // compiler!!!
    const ParticleType * nbhd_data = (ParticleType *) neighbors[lev].at(index).dataPtr();
    for(int i = 0; i < ng; i++){
        const ParticleType & part = nbhd_data[i];

        // Position of neighbour IBParticle
        RealVect pos = RealVect(
                AMREX_D_DECL(part.pos(0), part.pos(1), part.pos(2))
            );

        // Velocity of IBParticle
        RealVect vel = RealVect(
                AMREX_D_DECL(part.rdata(IBM_realData::velx),
                             part.rdata(IBM_realData::vely),
                             part.rdata(IBM_realData::velz)   )
            );

        // Velocity of IBParticle
        RealVect force = RealVect(
                AMREX_D_DECL(part.rdata(IBM_realData::forcex),
                             part.rdata(IBM_realData::forcey),
                             part.rdata(IBM_realData::forcez)   )
            );


        // Construct info struct
        IBM_info part_info;
        part_info.pos    = pos;
        part_info.vel    = vel;
        part_info.force  = force;
        part_info.id     = part.id();
        part_info.cpu    = part.cpu();
        part_info.real   = 0; // 0 => neighbor particle

        // Add to list

        if (unique) {
            // If in unique-mode: Don't add unless `part_info` is not already in `info`
            const auto & search = std::find(std::begin(info), std::end(info), part_info);
            if ( search == std::end(info)) {
#ifdef _OPENMP
#pragma omp critical
#endif
                { info.push_back(part_info); }
            }
        } else {
#ifdef _OPENMP
#pragma omp critical
#endif
            { info.push_back(part_info); }
        }
    }
}



Vector<IBM_info> IBMarkerContainer::NeighborIBMarkerInfo(int lev, PairIndex index,
                                                         bool unique) const {

    // Allocate Marker Info vector
    Vector<IBM_info> info;

    //___________________________________________________________________________
    // Fill Marker Info vector with neighbour data
    NeighborIBMarkerInfo(info, lev, index, unique);


    return info;
}



Vector<IBM_info> IBMarkerContainer::NeighborIBMarkerInfo(int lev, bool unique) const {

    // Allocate Marker Info vector
    Vector<IBM_info> info;

    //___________________________________________________________________________
    // Cell-centered MultiFab used as a reference for iterating over data
    // WARNING: this will break if IBMarkerContainer is on a differnt grid
    // than the grid which we're searching for particles (this should usually
    // be fine though)

    const BoxArray & ba            = ParticleBoxArray(lev);
    const DistributionMapping & dm = ParticleDistributionMap(lev);

    MultiFab dummy(ba, dm, 1, 1);

    //___________________________________________________________________________
    // Iterate over `dummy` looking for particles. NOTE: use the
    // IBMarkerContainer tile size
    for (MFIter mfi(dummy, tile_size); mfi.isValid(); ++mfi){
        PairIndex index(mfi.index(), mfi.LocalTileIndex());
        NeighborIBMarkerInfo(info, lev, index, true);
    }


    return info;
}



void IBMarkerContainer::IBMarkerInfo(Vector<IBM_info> & info, int lev, PairIndex index,
                                     bool unique) const {

    //___________________________________________________________________________
    // Fill Marker Info vector with local (non-neighbour) and neighbour data
       LocalIBMarkerInfo(info, lev, index, unique);
    NeighborIBMarkerInfo(info, lev, index, unique);
}



Vector<IBM_info> IBMarkerContainer::IBMarkerInfo(int lev, PairIndex index,
                                                 bool unique) const {

    // Allocate Marker Info vector
    Vector<IBM_info> info;

    //___________________________________________________________________________
    // Fill Marker Info vector with local (non-neighbour) and neighbour data
       LocalIBMarkerInfo(info, lev, index, unique);
    NeighborIBMarkerInfo(info, lev, index, unique);


    return info;
}



Vector<IBM_info> IBMarkerContainer::IBMarkerInfo(int lev, bool unique) const {

    // Allocate Marker Info vector
    Vector<IBM_info> info;

    //___________________________________________________________________________
    // Cell-centered MultiFab used as a reference for iterating over data
    // WARNING: this will break if IBMarkerContainer is on a differnt grid
    // than the grid which we're searching for particles (this should usually
    // be fine though)

    const BoxArray & ba            = ParticleBoxArray(lev);
    const DistributionMapping & dm = ParticleDistributionMap(lev);

    MultiFab dummy(ba, dm, 1, 1);

    //___________________________________________________________________________
    // Iterate over `dummy` looking for particles. NOTE: use the
    // IBMarkerContainer tile size
    for (MFIter mfi(dummy, tile_size); mfi.isValid(); ++mfi){
        PairIndex index(mfi.index(), mfi.LocalTileIndex());
        IBMarkerInfo(info, lev, index, true);
    }


    return info;
}



void IBMarkerContainer::InitInternals(int ngrow) {
    ReadStaticParameters();

    this->SetVerbose(0);

    // Turn off certain components for ghost particle communication
    // Field numbers: {0, 1, 2} => {x, y, z} particle coordinates
    //      => 3 corresponds to the start of IBM_realData
    setRealCommComp(4,  true);  // IBM_realData.velx
    setRealCommComp(5,  true);  // IBM_realData.vely
    setRealCommComp(6,  true);  // IBM_realData.velz
    setRealCommComp(7,  true);  // IBM_realData.forcex
    setRealCommComp(8,  true);  // IBM_realData.forcey
    setRealCommComp(9,  true);  // IBM_realData.forcez
    setRealCommComp(10, true);  // IBM_realData.pred_posx
    setRealCommComp(11, true);  // IBM_realData.pred_posy
    setRealCommComp(12, true);  // IBM_realData.pred_posz
    setRealCommComp(13, true);  // IBM_realData.pred_forcex
    setRealCommComp(14, true);  // IBM_realData.pred_forcey
    setRealCommComp(15, true);  // IBM_realData.pred_forcez

    // Field numbers: {0, 1} => {ID, CPU}
    //      => 2 corresponds to the start of IBM_intData
    // We _do_ want the the neighbour particles to have ID and cpu init data.
    setIntCommComp(2, true);  // IBM_intData.id_0
    setIntCommComp(3, true);  // IBM_intData.cpu_0
    setIntCommComp(4, true);  // IBM_intData.id_1
    setIntCommComp(5, true);  // IBM_intData.cpu_1



    /****************************************************************************
     *                                                                          *
     * Fill auxiallry data used by interpolsation                               *
     *   -> face_coords: the face-centered coordinates used by the fluid grids  *
     *                                                                          *
     ***************************************************************************/

    // TODO: this is only assuming 1 fluid level (level 0)
    int lev = 0;

    face_coords.resize(lev + 1);
    const BoxArray & ba            = ParticleBoxArray(lev);
    const DistributionMapping & dm = ParticleDistributionMap(lev);
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        const BoxArray ba_fc = convert(ba, nodal_flag_dir[d]);
        face_coords[lev][d].define(ba_fc, dm, AMREX_SPACEDIM, ngrow);
    }

    const Geometry & geom = Geom(lev);
    FindFaceCoords(face_coords[lev], geom);
}



void IBMarkerContainer::ReadStaticParameters() {
    static bool initialized = false;

    if (!initialized) {
        ParmParse pp("particles");

        // AMReX default is false => enable by default
        do_tiling = true;
        // Allow user to overwrite
        pp.query("do_tiling",  do_tiling);

        // If tiling is enabled, make sure that the tile size is at least the
        // number of ghost cells (otherwise strange things happen)
        if (do_tiling)
            tile_size = IntVect{AMREX_D_DECL(nghost, nghost, nghost)};
        // User can overwrite
        Vector<int> ts(BL_SPACEDIM);
        if (pp.queryarr("tile_size", ts))
            tile_size = IntVect(ts);

        pp.query("use_neighbor_list", use_neighbor_list);
        pp.query("sort_neighbor_list", sort_neighbor_list);

        initialized = true;
    }
}
