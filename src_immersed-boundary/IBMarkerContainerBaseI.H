#ifndef _IBMarkerContainerBaseI_H_
#define _IBMarkerContainerBaseI_H_



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



template <int NStructReal, int NStructInt>
IBMarkerContainerBase<NStructReal, NStructInt>::IBMarkerContainerBase(
            const Geometry & geom,
            const DistributionMapping & dmap,
            const BoxArray & ba,
            int n_nbhd
        ) : NeighborParticleContainer<NStructReal, NStructInt>(
            geom, dmap, ba, n_nbhd
        ),
    nghost(n_nbhd)
{
    InitInternals(n_nbhd);
}



template <int NStructReal, int NStructInt>
IBMarkerContainerBase<NStructReal, NStructInt>::IBMarkerContainerBase(
            AmrCore * amr_core, int n_nbhd
        ) : NeighborParticleContainer<NStructReal, NStructReal>(
            amr_core->GetParGDB(), n_nbhd
        ),
    m_amr_core(amr_core),
    nghost(n_nbhd)
{
    InitInternals(n_nbhd);
}



template <int NStructReal, int NStructInt>
void IBMarkerContainerBase<NStructReal, NStructInt>::MoveMarkers(int lev, Real dt) {

    for (IBMarIterBase<NStructReal, NStructInt> pti(* this, lev); pti.isValid(); ++pti) {

        PairIndex index(pti.index(), pti.LocalTileIndex());

        AoS & particles = GetParticles(lev).at(index).GetArrayOfStructs();
        long np = particles.size();

        for (int i=0; i<np; ++i) {
            ParticleType & part = particles[i];

            part.pos(0) += dt * part.rdata(IBM_realData::velx);
            part.pos(1) += dt * part.rdata(IBM_realData::vely);
            part.pos(2) += dt * part.rdata(IBM_realData::velz);
        }


        //ParticleVector & nbhd = GetNeighbors(lev, pti.index(), pti.LocalTileIndex());
        //long nn = nbhd.size();

        //for (int i=0; i<nn; ++i) {
        //    ParticleType & part = nbhd[i];

        //    part.pos(0) += dt * part.rdata(IBM_realData::velx);
        //    part.pos(1) += dt * part.rdata(IBM_realData::vely);
        //    part.pos(2) += dt * part.rdata(IBM_realData::velz);

        //    std::cout << "neighbor corrector velz is " << part.rdata(IBM_realData::velz);
        //    std::cout << "neighbor corrector posz is " << part.pos(2);

        //}
    }
}



template <int NStructReal, int NStructInt>
void IBMarkerContainerBase<NStructReal, NStructInt>::MovePredictor(int lev, Real dt) {

    for (IBMarIterBase<NStructReal, NStructInt> pti(* this, lev); pti.isValid(); ++pti) {
    // for (MFIter pti = MakeMFIter(lev, false); pti.isValid(); ++pti) {


        PairIndex index(pti.index(), pti.LocalTileIndex());

        AoS & particles = GetParticles(lev).at(index).GetArrayOfStructs();
        long np = particles.size();

        for (int i = 0; i < np; ++i) {
            ParticleType & part = particles[i];

            // // update predictor to match the position
            // part.rdata(IBM_realData::pred_posx) = part.pos(0);
            // part.rdata(IBM_realData::pred_posy) = part.pos(1);
            // part.rdata(IBM_realData::pred_posz) = part.pos(2);

            // part.rdata(IBM_realData::pred_posx) += dt * part.rdata(IBM_realData::pred_velx);
            // part.rdata(IBM_realData::pred_posy) += dt * part.rdata(IBM_realData::pred_vely);
            // part.rdata(IBM_realData::pred_posz) += dt * part.rdata(IBM_realData::pred_velz);
            part.rdata(IBM_realData::pred_posx) = dt * part.rdata(IBM_realData::pred_velx);
            part.rdata(IBM_realData::pred_posy) = dt * part.rdata(IBM_realData::pred_vely);
            part.rdata(IBM_realData::pred_posz) = dt * part.rdata(IBM_realData::pred_velz);
        }


        // ParticleVector & nbhd = GetNeighbors(lev, pti.index(), pti.LocalTileIndex());
        // long nn = nbhd.size();

        // for (int i=0; i<nn; ++i) {
        //     ParticleType & part = nbhd[i];

        //     // part.rdata(IBM_realData::pred_posx) = part.pos(0);
        //     // part.rdata(IBM_realData::pred_posy) = part.pos(1);
        //     // part.rdata(IBM_realData::pred_posz) = part.pos(2);

        //     // part.rdata(IBM_realData::pred_posx) += dt * part.rdata(IBM_realData::pred_velx);
        //     // part.rdata(IBM_realData::pred_posy) += dt * part.rdata(IBM_realData::pred_vely);
        //     // part.rdata(IBM_realData::pred_posz) += dt * part.rdata(IBM_realData::pred_velz);
        //     part.rdata(IBM_realData::pred_posx) = dt * part.rdata(IBM_realData::pred_velx);
        //     part.rdata(IBM_realData::pred_posy) = dt * part.rdata(IBM_realData::pred_vely);
        //     part.rdata(IBM_realData::pred_posz) = dt * part.rdata(IBM_realData::pred_velz);
        // }
    }
}



template <int NStructReal, int NStructInt>
void IBMarkerContainerBase<NStructReal, NStructInt>::ResetMarkers(int lev) {

    for (IBMarIterBase<NStructReal, NStructInt> pti(* this, lev); pti.isValid(); ++pti) {

        PairIndex index(pti.index(), pti.LocalTileIndex());
        auto & particle_data = GetParticles(lev).at(index);
        long np = particle_data.size();

        AoS & particles = particle_data.GetArrayOfStructs();
        for (int i = 0; i < np; ++i) {
            ParticleType & part = particles[i];

            part.rdata(IBM_realData::velx) = 0.;
            part.rdata(IBM_realData::vely) = 0.;
            part.rdata(IBM_realData::velz) = 0.;

            part.rdata(IBM_realData::forcex) = 0.;
            part.rdata(IBM_realData::forcey) = 0.;
            part.rdata(IBM_realData::forcez) = 0.;
        }
    }
}



template <int NStructReal, int NStructInt>
void IBMarkerContainerBase<NStructReal, NStructInt>::ResetPredictor(int lev) {

    for (IBMarIterBase<NStructInt, NStructInt> pti(* this, lev); pti.isValid(); ++pti) {

        PairIndex index(pti.index(), pti.LocalTileIndex());
        auto & particle_data = GetParticles(lev).at(index);
        long np = particle_data.size();

        AoS & particles = particle_data.GetArrayOfStructs();
        for (int i = 0; i < np; ++i) {
            ParticleType & part = particles[i];

            part.rdata(IBM_realData::pred_velx) = 0.;
            part.rdata(IBM_realData::pred_vely) = 0.;
            part.rdata(IBM_realData::pred_velz) = 0.;

            part.rdata(IBM_realData::pred_forcex) = 0.;
            part.rdata(IBM_realData::pred_forcey) = 0.;
            part.rdata(IBM_realData::pred_forcez) = 0.;
        }
    }
}



template <int NStructReal, int NStructInt>
void IBMarkerContainerBase<NStructReal, NStructInt>::SpreadMarkers(
            int lev,
            const Vector<RealVect> & f_in,
            const Vector<RealVect> & f_pos,
            std::array<MultiFab, AMREX_SPACEDIM> & f_out,
            std::array<MultiFab, AMREX_SPACEDIM> & f_weights
        ) const {


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

        Box bx       = mfi.growntilebox();
        Box tile_box = mfi.tilebox();

        spread_markers(BL_TO_FORTRAN_BOX(bx),
                       BL_TO_FORTRAN_BOX(tile_box),
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



template <int NStructReal, int NStructInt>
void IBMarkerContainerBase<NStructReal, NStructInt>::SpreadMarkers(
            int lev,
            std::array<MultiFab, AMREX_SPACEDIM> & f_out
        ) const {


    //___________________________________________________________________________
    // We don't need these spreading weights => create a dummy MF
    std::array<MultiFab, AMREX_SPACEDIM> f_weights;
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        f_weights[d].define(f_out[d].boxArray(), f_out[d].DistributionMap(),
                            1, f_out[d].nGrow());
        f_weights[d].setVal(0.);
    }


    // for (MFIter pti = MakeMFIter(lev, true); pti.isValid(); ++pti) {
    for (MyConstParIter pti(* this, lev); pti.isValid(); ++pti) {

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
    }
}



template <int NStructReal, int NStructInt>
void IBMarkerContainerBase<NStructReal, NStructInt>::SpreadPredictor(
            int lev,
            std::array<MultiFab, AMREX_SPACEDIM> & f_out
        ) const {


    //___________________________________________________________________________
    // We don't need these spreading weights => create a dummy MF
    std::array<MultiFab, AMREX_SPACEDIM> f_weights;
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        f_weights[d].define(f_out[d].boxArray(), f_out[d].DistributionMap(),
                            1, f_out[d].nGrow());
        f_weights[d].setVal(0.);
    }


    // for (MFIter pti = MakeMFIter(lev, true); pti.isValid(); ++pti) {
    for (MyConstParIter pti(* this, lev); pti.isValid(); ++pti) {

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

            ppos[0] = part.pos(0) + part.rdata(IBM_realData::pred_posx);
#if (AMREX_SPACEDIM > 1)
            ppos[1] = part.pos(1) + part.rdata(IBM_realData::pred_posy);
#endif
#if (AMREX_SPACEDIM > 2)
            ppos[2] = part.pos(2) + part.rdata(IBM_realData::pred_posz);
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
    }
}



template <int NStructReal, int NStructInt>
void IBMarkerContainerBase<NStructReal, NStructInt>::InterpolateMarkers(
            int lev,
            Vector<RealVect> & f_out,
            const Vector<RealVect> & f_pos,
            const std::array<MultiFab, AMREX_SPACEDIM> & f_in,
            const std::array<MultiFab, AMREX_SPACEDIM> & f_weights
        ) const {


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

    Box bx       = mfi.growntilebox();
    Box tile_box = mfi.tilebox();

    interpolate_markers(BL_TO_FORTRAN_BOX(bx),
                        BL_TO_FORTRAN_BOX(tile_box),
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



template <int NStructReal, int NStructInt>
void IBMarkerContainerBase<NStructReal, NStructInt>::InterpolateMarkers(
            int lev,
            const std::array<MultiFab, AMREX_SPACEDIM> & f_out
        ) {


    //___________________________________________________________________________
    // We don't need these spreading weights => create a dummy MF
    std::array<MultiFab, AMREX_SPACEDIM> f_weights;
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        f_weights[d].define(f_out[d].boxArray(), f_out[d].DistributionMap(),
                            1, f_out[d].nGrow());
        f_weights[d].setVal(-1.);
    }


    // for (MFIter pti = MakeMFIter(lev, true); pti.isValid(); ++pti) {
    for (MyConstParIter pti(* this, lev); pti.isValid(); ++pti) {

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
    }
}



template <int NStructReal, int NStructInt>
void IBMarkerContainerBase<NStructReal, NStructInt>::InterpolatePredictor(
            int lev,
            const std::array<MultiFab, AMREX_SPACEDIM> & f_out
        ) {


    //___________________________________________________________________________
    // We don't need these spreading weights => create a dummy MF
    std::array<MultiFab, AMREX_SPACEDIM> f_weights;
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        f_weights[d].define(f_out[d].boxArray(), f_out[d].DistributionMap(),
                            1, f_out[d].nGrow());
        f_weights[d].setVal(-1.);
    }


    // for (MFIter pti = MakeMFIter(lev, false); pti.isValid(); ++pti) {
    for (MyConstParIter pti(* this, lev); pti.isValid(); ++pti) {

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

            ppos[0] = part.pos(0) + part.rdata(IBM_realData::pred_posx);
#if (AMREX_SPACEDIM > 1)
            ppos[1] = part.pos(1) + part.rdata(IBM_realData::pred_posy);
#endif
#if (AMREX_SPACEDIM > 2)
            ppos[2] = part.pos(2) + part.rdata(IBM_realData::pred_posz);
#endif

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

            part.rdata(IBM_realData::pred_velx) += marker_forces[i][0];
#if (AMREX_SPACEDIM > 1)
            part.rdata(IBM_realData::pred_vely) += marker_forces[i][1];
#endif
#if (AMREX_SPACEDIM > 2)
            part.rdata(IBM_realData::pred_velz) += marker_forces[i][2];
#endif
        }
    }
}



template <int NStructReal, int NStructInt>
void IBMarkerContainerBase<NStructReal, NStructInt>::PrintMarkerData(int lev) const {

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



template <int NStructReal, int NStructInt>
void IBMarkerContainerBase<NStructReal, NStructInt>::InitInternals(
            int ngrow, int real_comp, int int_comp
        ) {

    ReadStaticParameters();

    this->SetVerbose(0);

    // Turn off certain components for ghost particle communication
    // Field numbers: {0, 1, 2} => {x, y, z} particle coordinates
    //      => 3 corresponds to the start of IBM_realData
    for (int i=3; i < real_comp + 3; ++i)
        setRealCommComp(i,  true);

    // Field numbers: {0, 1} => {ID, CPU}
    //      => 2 corresponds to the start of IBM_intData
    // We _do_ want the the neighbour particles to have ID and cpu init data.
    for (int i = 2; i < int_comp + 2; ++i)
        setIntCommComp(i, true);



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



template <int NStructReal, int NStructInt>
void IBMarkerContainerBase<NStructReal, NStructInt>::ReadStaticParameters() {
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

        // pp.query("use_neighbor_list", use_neighbor_list);
        // pp.query("sort_neighbor_list", sort_neighbor_list);

        initialized = true;
    }
}


#endif
