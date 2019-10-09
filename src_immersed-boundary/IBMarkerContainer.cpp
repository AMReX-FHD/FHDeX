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

#include <iostream>

using namespace common;
using namespace amrex;


bool IBMarkerContainer::use_neighbor_list  {true};
bool IBMarkerContainer::sort_neighbor_list {false};



IBMarkerContainer::IBMarkerContainer(const Geometry & geom,
                                         const DistributionMapping & dmap,
                                         const BoxArray & ba,
                                         int n_nbhd)
    : IBMarkerContainerBase<IBM_realData, IBM_intData>(
            geom, dmap, ba, n_nbhd
        )
      , n_list(0)
{
    InitInternals(n_nbhd);
    nghost = n_nbhd;
}



IBMarkerContainer::IBMarkerContainer(AmrCore * amr_core, int n_nbhd)
    : IBMarkerContainerBase<IBM_realData, IBM_intData>(
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

                // Set particle position
                p_new.pos(0) = pos[i][0];
                p_new.pos(1) = pos[i][1];
                p_new.pos(2) = pos[i][2];

                p_new.rdata(IBM_realData::radius) = radius[i];

                // Initialize marker velocity as well as forces to 0
                p_new.rdata(IBM_realData::velx)   = 0.;
                p_new.rdata(IBM_realData::vely)   = 0.;
                p_new.rdata(IBM_realData::velz)   = 0.;

                p_new.rdata(IBM_realData::forcex) = 0.;
                p_new.rdata(IBM_realData::forcey) = 0.;
                p_new.rdata(IBM_realData::forcez) = 0.;

                p_new.rdata(IBM_realData::pred_posx)   = 0.;
                p_new.rdata(IBM_realData::pred_posy)   = 0.;
                p_new.rdata(IBM_realData::pred_posz)   = 0.;

                p_new.rdata(IBM_realData::pred_velx)   = 0.;
                p_new.rdata(IBM_realData::pred_vely)   = 0.;
                p_new.rdata(IBM_realData::pred_velz)   = 0.;

                p_new.rdata(IBM_realData::pred_forcex) = 0.;
                p_new.rdata(IBM_realData::pred_forcey) = 0.;
                p_new.rdata(IBM_realData::pred_forcez) = 0.;

                // These are filled in the next sweep:
                p_new.idata(IBM_intData::id_0)  = -1;
                p_new.idata(IBM_intData::cpu_0) = -1;

                // ID_1 remembers the particle's position in the linked list
                p_new.idata(IBM_intData::id_1)  = i;
                p_new.idata(IBM_intData::cpu_1) = i_ib; // label immersed boundaries

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
    fillNeighbors();

    for (IBMarIter pti(*this, lev); pti.isValid(); ++pti) {

        // Get marker data (local to current thread)
        PairIndex index(pti.index(), pti.LocalTileIndex());
        AoS & markers = GetParticles(lev).at(index).GetArrayOfStructs();

        // Get neighbor marker data (from neighboring threads)
        ParticleVector & nbhd_data = GetNeighbors(lev, pti.index(), pti.LocalTileIndex());

        long np = markers.size();
        long nn = nbhd_data.size();

        // Sweep over particles, check N^2 candidates for previous list member
        for (int i=0; i<np; ++i) {

            ParticleType & mark = markers[i];

            // Check other (real) particles in tile
            for (int j=0; j<np; ++j) {

                const ParticleType & other = markers[j];
                if ((mark.idata(IBM_intData::id_1) == other.idata(IBM_intData::id_1) + 1)
                    && (mark.idata(IBM_intData::cpu_1) == other.idata(IBM_intData::cpu_1)))
                {
                    mark.idata(IBM_intData::id_0)  = other.id();
                    mark.idata(IBM_intData::cpu_0) = other.cpu();
                }
            }

            // Check neighbor particles
            for (int j=0; j<nn; ++j) {

                const ParticleType & other = nbhd_data[j];
                if ((mark.idata(IBM_intData::id_1) == other.idata(IBM_intData::id_1) + 1)
                    && (mark.idata(IBM_intData::cpu_1) == other.idata(IBM_intData::cpu_1)))
                {
                    mark.idata(IBM_intData::id_0)  = other.id();
                    mark.idata(IBM_intData::cpu_0) = other.cpu();
                }
            }
        }
    }
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

            // Set particle position
            p_new.pos(0) = pos[0];
            p_new.pos(1) = pos[1];
            p_new.pos(2) = pos[2];

            p_new.rdata(IBM_realData::radius) = radius;

            // Initialize marker velocity as well as forces to 0
            p_new.rdata(IBM_realData::velx)   = 0.;
            p_new.rdata(IBM_realData::vely)   = 0.;
            p_new.rdata(IBM_realData::velz)   = 0.;

            p_new.rdata(IBM_realData::forcex) = 0.;
            p_new.rdata(IBM_realData::forcey) = 0.;
            p_new.rdata(IBM_realData::forcez) = 0.;

            p_new.rdata(IBM_realData::pred_posx)   = 0.;
            p_new.rdata(IBM_realData::pred_posy)   = 0.;
            p_new.rdata(IBM_realData::pred_posz)   = 0.;

            p_new.rdata(IBM_realData::pred_velx)   = 0.;
            p_new.rdata(IBM_realData::pred_vely)   = 0.;
            p_new.rdata(IBM_realData::pred_velz)   = 0.;

            p_new.rdata(IBM_realData::pred_forcex) = 0.;
            p_new.rdata(IBM_realData::pred_forcey) = 0.;
            p_new.rdata(IBM_realData::pred_forcez) = 0.;

            p_new.idata(IBM_intData::id_0)  = id;
            p_new.idata(IBM_intData::cpu_0) = cpu;

            p_new.idata(IBM_intData::id_1)  = i_ref;
            p_new.idata(IBM_intData::cpu_1) = -1;

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



int IBMarkerContainer::FindConnectedMarkers(      AoS & particles,
                                            const ParticleType & part,
                                                  ParticleVector & nbhd_data,
                                            const IntVector & nbhd,
                                            int nbhd_index,
                                            ParticleType *& prev_marker,
                                            ParticleType *& next_marker) const {

    BL_PROFILE_VAR("IBMarkerContainer::FindConnectedMarkers", FindNeighbors);

    long np = particles.size();
    int nn  = nbhd[nbhd_index]; // number of neighbors for particle at nbhd_index
    nbhd_index ++; // pointing at first neighbor


    bool prev_set = false;
    bool next_set = false;

    // Loops over neighbor list
    for (int j=0; j < nn; ++j) {
        int ni = nbhd[nbhd_index] - 1; // -1 <= neighbor list uses Fortran indexing

        ParticleType * npart;
        if (ni >= np) {
            ni = ni - np;
            npart = & nbhd_data[ni];
        } else {
            npart = & particles[ni];
        }

        // Check if the neighbor candidate is the previous/minus neighbor
        if (        (npart->id() == part.idata(IBM_intData::id_0))
                && (npart->cpu() == part.idata(IBM_intData::cpu_0)) ) {

            prev_marker = npart;
            prev_set = true;
        }

        // Check if the neighbor candidate is the next/plus neighbor
        if (        (part.id() == npart->idata(IBM_intData::id_0))
                && (part.cpu() == npart->idata(IBM_intData::cpu_0)) ) {

            next_marker = npart;
            next_set = true;
        }

        nbhd_index ++;
    }

    if (prev_set && next_set) {
        return 0;
    } else if (! prev_set &&   next_set ) {
        return 1;
    } else if (  prev_set && ! next_set)  {
        return 2;
    } else {
        return -1;
    }

    BL_PROFILE_VAR_STOP(FindNeighbors);
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
