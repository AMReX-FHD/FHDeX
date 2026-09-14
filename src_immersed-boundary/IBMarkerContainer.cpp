#include <AMReX.H>
#include <AMReX_AmrParGDB.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Particles.H>
#include <AMReX_NeighborParticles.H>

#include <AMReX_VisMF.H>  // amrex::VisMF::Write(MultiFab)

#include <common_functions.H>

#include <IBMarkerContainer.H>
#include <ib_functions_F.H>

#include <iostream>

using namespace amrex;


bool IBMarkerContainer::use_neighbor_list  {true};
bool IBMarkerContainer::sort_neighbor_list {false};



IBMarkerContainer::IBMarkerContainer(const Geometry & geom,
                                         const DistributionMapping & dmap,
                                         const BoxArray & ba,
                                         int n_nbhd)
    : IBMarkerContainerBase<IBMReal, IBMInt>(
            geom, dmap, ba, n_nbhd
        )
      , n_list(0)
{
    InitInternals(n_nbhd);
    nghost = n_nbhd;
}



IBMarkerContainer::IBMarkerContainer(AmrCore * amr_core, int n_nbhd)
    : IBMarkerContainerBase<IBMReal, IBMInt>(
            amr_core->GetParGDB(), n_nbhd
        )
      , n_list(0)
{
    InitInternals(n_nbhd);
    nghost     = n_nbhd;
    m_amr_core = amr_core;
}



void IBMarkerContainer::InitList(int lev,
                                 const Vector<Real> & radius,
                                 const Vector<RealVect> & pos,
                                 int i_ib) {

    // Inverse cell-size vector => used for determining index corresponding to
    // IBParticle position (pos)
    RealVect inv_dx = RealVect(AMREX_D_DECL(Geom(lev).InvCellSize(0),
                                            Geom(lev).InvCellSize(1),
                                            Geom(lev).InvCellSize(2)   ));


    int total_np = 0;


    //__________________________________________________________________________
    // First sweep: generate particles in their respective cores. Do not link
    // yet because linkage target could be a neighbor particle (which aren't
    // filled yet)

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


        for(int i=0; i<pos.size(); i++) {
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

                for (int d=0; d<AMREX_SPACEDIM; ++d) {
                    // Set particle (marker) position
                    p_new.pos(d) = pos[i][d];

                    // Initialize marker velocity as well as forces to 0
                    p_new.rdata(IBMReal::velx + d)   = 0.;
                    p_new.rdata(IBMReal::forcex + d) = 0.;

                    p_new.rdata(IBMReal::pred_posx + d)   = 0.;
                    p_new.rdata(IBMReal::pred_velx + d)   = 0.;
                    p_new.rdata(IBMReal::pred_forcex + d) = 0.;
                }

                // Marker metadata (pt #2 are filled in the next sweep):
                // 1. Marker search radius
                p_new.rdata(IBMReal::radius) = radius[i];

                // 2. Marker contexual (and connectivity) metadata
                p_new.idata(IBMInt::id_0)  = -1;
                p_new.idata(IBMInt::cpu_0) = -1;

                // ID_1 remembers the particle's position in the linked list
                p_new.idata(IBMInt::id_1)  = i;
                p_new.idata(IBMInt::cpu_1) = i_ib; // label immersed boundaries

                // Add to the data structure
                particles.push_back(p_new);
            }
        }

        const int np = pcount;
        total_np += np;
    }

    ParallelDescriptor::ReduceIntSum(total_np, ParallelDescriptor::IOProcessorNumber());
    Print() << "Total number of generated particles: " << total_np << std::endl;



    //__________________________________________________________________________
    // Second Sweep: At this point, each particle's position in the original
    // list is its linkage identifier => sweep list to link up particles
    // properly

    Redistribute();
    clearNeighbors();
    // fillNeighbors();

    // for (IBMarIter pti(*this, lev); pti.isValid(); ++pti) {

    //     // Get marker data (local to current thread)
    //     PairIndex index(pti.index(), pti.LocalTileIndex());
    //     AoS & markers = GetParticles(lev).at(index).GetArrayOfStructs();

    //     // Get neighbor marker data (from neighboring threads)
    //     ParticleVector & nbhd_data = GetNeighbors(lev, pti.index(), pti.LocalTileIndex()).GetArrayOfStructs()();

    //     long np = pti.numParticles();
    //     long nn = nbhd_data.size();

    //     // Sweep over particles, check N^2 candidates for previous list member
    //     for (int i=0; i<np; ++i) {

    //         ParticleType & mark = markers[i];

    //         // Check other (real) particles in tile
    //         for (int j=0; j<np; ++j) {

    //             const ParticleType & other = markers[j];
    //             if ((mark.idata(IBMInt::id_1) == other.idata(IBMInt::id_1) + 1)
    //                 && (mark.idata(IBMInt::cpu_1) == other.idata(IBMInt::cpu_1)))
    //             {
    //                 mark.idata(IBMInt::id_0)  = other.id();
    //                 mark.idata(IBMInt::cpu_0) = other.cpu();
    //             }
    //         }

    //         // Check neighbor particles
    //         for (int j=0; j<nn; ++j) {

    //             const ParticleType & other = nbhd_data[j];
    //             if ((mark.idata(IBMInt::id_1) == other.idata(IBMInt::id_1) + 1)
    //                 && (mark.idata(IBMInt::cpu_1) == other.idata(IBMInt::cpu_1)))
    //             {
    //                 mark.idata(IBMInt::id_0)  = other.id();
    //                 mark.idata(IBMInt::cpu_0) = other.cpu();
    //             }
    //         }
    //     }
    // }
}



void IBMarkerContainer::InitSingle(int lev, Real radius, const RealVect & pos,
                                   int id, int cpu, int i_ref) {

    // Inverse cell-size vector => used for determining index corresponding to
    // IBParticle position (pos)
    RealVect inv_dx = RealVect(AMREX_D_DECL(Geom(lev).InvCellSize(0),
                                            Geom(lev).InvCellSize(1),
                                            Geom(lev).InvCellSize(2) ));


    // This uses the particle tile size. Note that the default is to tile so if
    // we remove the true and don't explicitly add false it will still tile.
    for (MFIter mfi = MakeMFIter(lev, true); mfi.isValid(); ++mfi) {

        // Current tile box
        const Box & tile_box = mfi.tilebox();

        // Create a particle container for this grid and add the
        // immersed-boundary particles to it if the particle's position (pos)
        // is within current tile box.
        const int grid_id = mfi.index();
        const int tile_id = mfi.LocalTileIndex();
        auto & particles = GetParticles(lev)[std::make_pair(grid_id,tile_id)];


        // IntVect representing particle's position in the tile_box grid.
        RealVect pos_grid = pos; // Important: need to initialize on same CPU
        pos_grid *= inv_dx;
        IntVect pos_ind = IntVect(AMREX_D_DECL((int) pos_grid[0],
                                               (int) pos_grid[1],
                                               (int) pos_grid[2] ));

        // Add particle at position pos iff it's vector index is contained
        // within tile_box.
        if(tile_box.contains(pos_ind)) {

            ParticleType p_new;

            // Set id and cpu for this particle
            p_new.id()  = ParticleType::NextID();
            p_new.cpu() = ParallelDescriptor::MyProc();

            for (int d=0; d<AMREX_SPACEDIM; ++d) {
                // Set particle (marker) position
                p_new.pos(d) = pos[d];

                // Initialize marker velocity as well as forces to 0
                p_new.rdata(IBMReal::velx + d)   = 0.;
                p_new.rdata(IBMReal::forcex + d) = 0.;

                p_new.rdata(IBMReal::pred_posx + d)   = 0.;
                p_new.rdata(IBMReal::pred_velx + d)   = 0.;
                p_new.rdata(IBMReal::pred_forcex + d) = 0.;
            }

            // Marker metadata:
            // 1. Marker search radius
            p_new.rdata(IBMReal::radius) = radius;

            // 2. Marker contexual (and connectivity) metadata
            p_new.idata(IBMInt::id_0)  = id;
            p_new.idata(IBMInt::cpu_0) = cpu;

            p_new.idata(IBMInt::id_1)  = i_ref;
            p_new.idata(IBMInt::cpu_1) = -1;

            // Add to the data structure
            particles.push_back(p_new);
        }
    }

    // We shouldn't need this if the particles are tiled with one tile per
    // grid, but otherwise we do need this to move particles from tile 0 to the
    // correct tile.
    Redistribute();
}



void IBMarkerContainer::SpreadMarkers(int lev,
                                      const Vector<RealVect> & f_in,
                                      std::array<MultiFab, AMREX_SPACEDIM> & f_out,
                                      std::array<MultiFab, AMREX_SPACEDIM> & f_weights) const {

    //___________________________________________________________________________
    // Geometry data
    const Geometry & geom = this->Geom(0);
    const Real     *   dx = geom.CellSize();

    //___________________________________________________________________________
    // Fill vector of marker positions (for current level)
    Vector<IBM_info> marker_info = IBMarkerInfo(lev);
    Vector<RealVect> marker_positions(marker_info.size());
    for (int i=0; i<marker_info.size(); ++i)
        marker_positions[i] = marker_info[i].pos;


    SpreadMarkers(f_in, marker_positions, f_out, f_weights,
                  get_face_coords(lev), dx, 0);
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



void IBMarkerContainer::InterpolateMarkers(int lev,
                                           Vector<RealVect> & f_out,
                                           const std::array<MultiFab, AMREX_SPACEDIM> & f_in,
                                           const std::array<MultiFab, AMREX_SPACEDIM> & f_weights) const {

    //___________________________________________________________________________
    // Geometry data
    const Geometry & geom = this->Geom(0);
    const Real     *   dx = geom.CellSize();

    //___________________________________________________________________________
    // Fill vector of marker positions (for current level)
    Vector<IBM_info> marker_info = IBMarkerInfo(lev);
    Vector<RealVect> marker_positions(marker_info.size());
    for (int i=0; i<marker_info.size(); ++i)
        marker_positions[i] = marker_info[i].pos;

    InterpolateMarkers(f_out, marker_positions, f_in, f_weights,
                       get_face_coords(lev), dx, 0);
}



void IBMarkerContainer::InterpolateMarkers(int lev,
                       Vector<RealVect> & f_out,
                       const std::array<MultiFab, AMREX_SPACEDIM> & f_in) const {

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



int IBMarkerContainer::ConnectedMarkers(
            int lev, const TileIndex & tile, MarkerListIndex & part_index,
            ParticleType *& prev_marker,     ParticleType *& next_marker
        ) {

    BL_PROFILE("IBMarkerContainer::ConnectedMarkers");

    // Get marker data
    AoS & particles = GetParticles(lev).at(tile).GetArrayOfStructs();
    ParticleType & part = particles[part_index.first];
    long np = GetParticles(lev).at(tile).numParticles();

    // Get neighbor marker data (from neighboring threads)
    ParticleVector & nbhd_data = GetNeighbors(lev, tile.first, tile.second).GetArrayOfStructs()();

    // Get neighbor list (for collision checking)
    std::pair<int, int> index = std::make_pair(tile.first, tile.second);
    auto nbor_data = m_neighbor_list[lev][index].data();

    // Point to the right spot in the neighbor list
    bool prev_set = false;
    bool next_set = false;

    int nn = part_index.second;
    for (auto & p2 : nbor_data.getNeighbors(part_index.first)) {
        ParticleType * npart = & p2; // Get pointer to neighbor particle

        // Check if the neighbor candidate is the previous/minus neighbor
        if (        (npart->id() == part.idata(IBMInt::id_0))
                && (npart->cpu() == part.idata(IBMInt::cpu_0)) ) {

            prev_marker = npart;
            prev_set = true;
        }

        // Check if the neighbor candidate is the next/plus neighbor
        if (        (part.id() == npart->idata(IBMInt::id_0))
                && (part.cpu() == npart->idata(IBMInt::cpu_0)) ) {

            next_marker = npart;
            next_set = true;
        }

        nn ++;
    }

    // return new index in neighbor list
    part_index.second += nn + 1;

    if (prev_set && next_set) {
        return 0;
    } else if (! prev_set &&   next_set ) {
        return 1;
    } else if (  prev_set && ! next_set)  {
        return 2;
    } else {
        return -1;
    }
//    BL_PROFILE_VAR_STOP(FindNeighbors);
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
    long np = particle_data.numParticles();

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
                AMREX_D_DECL(part.rdata(IBMReal::velx),
                             part.rdata(IBMReal::vely),
                             part.rdata(IBMReal::velz)   )
            );

        // Velocity of IBParticle
        RealVect force = RealVect(
                AMREX_D_DECL(part.rdata(IBMReal::forcex),
                             part.rdata(IBMReal::forcey),
                             part.rdata(IBMReal::forcez)   )
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

    RealVect inv_dx = RealVect(AMREX_D_DECL(Geom(lev).InvCellSize(0),
                                            Geom(lev).InvCellSize(1),
                                            Geom(lev).InvCellSize(2)  ));

    int ng = neighbors[lev].at(index).size();

    // Iterate over neighbour particles:
    // TODO: HAXOR!!! This should be fixed ASAP: if I understand this correctly,
    // the neighbor data contains the particle data as a binary array (char). By
    // casting to ParticleType, what we're doing is interpreting the data in
    // neighbours[index] as valid particle data. Also we stride the
    // neighbors[index] array in units of sizeof(ParticleData). All of this is a
    // little too dangerous for my taste: never hide what you're doing from your
    // compiler!!!
    const ParticleType * nbhd_data = (ParticleType *) neighbors[lev].at(index).GetArrayOfStructs().dataPtr();
    for(int i = 0; i < ng; i++){
        const ParticleType & part = nbhd_data[i];

        // Position of neighbour IBParticle
        RealVect pos = RealVect(
                AMREX_D_DECL(part.pos(0), part.pos(1), part.pos(2))
            );

        // Velocity of IBParticle
        RealVect vel = RealVect(
                AMREX_D_DECL(part.rdata(IBMReal::velx),
                             part.rdata(IBMReal::vely),
                             part.rdata(IBMReal::velz)   )
            );

        // Velocity of IBParticle
        RealVect force = RealVect(
                AMREX_D_DECL(part.rdata(IBMReal::forcex),
                             part.rdata(IBMReal::forcey),
                             part.rdata(IBMReal::forcez)   )
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



void IBMarkerContainer::UpdatePIDMap() {
    // Calls UpdatePIDMap from base class (IBMarkerContainerBase)
    IBMarkerContainerBase<IBMReal, IBMInt>::UpdatePIDMap();

    Vector<int> ibs(getTotalNumIDs());
    PullDownInt(0, ibs, IBMInt::cpu_1);

    Vector<int> ids(getTotalNumIDs());
    PullDownInt(0, ids, IBMInt::id_1);

    for (int i = 0; i < ibs.size(); ++i) {
        sorted_map.push_back(std::make_tuple(ibs[i], ids[i], i));
    }

    sort(sorted_map.begin(), sorted_map.end());

    Print() << "Flagellum number\t" << "index in PullDown Vector" << std::endl;
    for (int i = 0; i < ibs.size(); i++) {
        Print() << immbdy_idx(sorted_map[i])
                << "\t" << marker_idx(sorted_map[i])
                << "\t" << storage_idx(sorted_map[i])
                << std::endl;
    }

    //Create a reduced vector storing only the beginning index of each flagellum in the sorted map
    int n_ibs =1 + ibs[std::distance(ibs.begin(), std::max_element(ibs.begin(), ibs.end()))];

    for (int i=0; i < n_ibs; i++) {
        reduced_map.push_back(
            std::distance(
                sorted_map.begin(),
                std::find_if(
                    sorted_map.begin(),
                    sorted_map.end(),
                    [&](const auto & t) { return immbdy_idx(t) == i; }
                )
            )
        );
    }

    Print() << "Flagellum number\t" << "Beginning index in the sorted map" << std::endl;
    for (int i = 0; i < n_ibs; i++) {
        Print() << i << "\t" << reduced_map[i] << std::endl;
    }
}
