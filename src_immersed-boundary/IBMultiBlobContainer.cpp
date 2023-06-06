#include <AMReX.H>
#include <AMReX_AmrCore.H>
#include <AMReX_AmrParGDB.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Particles.H>
#include <AMReX_NeighborParticles.H>

#include <AMReX_VisMF.H>  // amrex::VisMF::Write(MultiFab)

#include <common_functions.H>

#include <IBMultiBlobContainer.H>
#include <ib_functions_F.H>



using namespace amrex;



BlobContainer::BlobContainer(const Geometry & geom,
                                         const DistributionMapping & dmap,
                                         const BoxArray & ba,
                                         int n_nbhd)
    : IBMarkerContainerBase<IBBReal, IBBInt>(
            geom, dmap, ba, n_nbhd
        )
{
    InitInternals(n_nbhd);
    nghost = n_nbhd;
}



BlobContainer::BlobContainer(AmrCore * amr_core, int n_nbhd)
    : IBMarkerContainerBase<IBBReal, IBBInt>(
            amr_core->GetParGDB(), n_nbhd
        )
{
    InitInternals(n_nbhd);
    nghost     = n_nbhd;
    m_amr_core = amr_core;
}



void BlobContainer::AddSingle(int lev, const TileIndex & tile, Real radius,
                              Real k_spring, const RealVect & pos, int id,
                              int cpu, int i_ref) {


        // Create a particle container for this grid and add the
        // immersed-boundary particles to it. NOTE: use the square-bracket
        // operator here, because we might need to create this tile (tiles are
        // created as needed)
        auto & particles = GetParticles(lev)[tile];

        ParticleType p_new;

        // Set id and cpu for this particle
        p_new.id()  = ParticleType::NextID();
        p_new.cpu() = ParallelDescriptor::MyProc();

        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            // Set particle (blob) and reference position
            p_new.pos(d) = pos[d];
            p_new.rdata(IBBReal::ref_delx + d) = 0;

            // Initialize marker velocity as well as forces to 0
            p_new.rdata(IBBReal::velx + d)   = 0.;
            p_new.rdata(IBBReal::forcex + d) = 0.;

            p_new.rdata(IBBReal::pred_posx + d)   = 0.;
            p_new.rdata(IBBReal::pred_velx + d)   = 0.;
            p_new.rdata(IBBReal::pred_forcex + d) = 0.;
        }

        // Blob metadata
        // 1. Blob search radius
        p_new.rdata(IBBReal::radius) = radius;

        // 2. Blob (anchoring) spring stiffness
        p_new.rdata(IBBReal::k_spring) = k_spring;

        // 2. Blob contexual metadata
        p_new.idata(IBBInt::id_0)  = id;
        p_new.idata(IBBInt::cpu_0) = cpu;

        p_new.idata(IBBInt::id_1)  = i_ref;
        p_new.idata(IBBInt::cpu_1) = -1;

        // Add to the data structure
        particles.push_back(p_new);
}



void BlobContainer::AddSingle(int lev, Real radius, Real k_spring,
                              const RealVect & pos, int id, int cpu, int i_ref) {

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


        // IntVect representing particle's position in the tile_box grid.
        RealVect pos_grid = pos; // Important: need to initialize on same CPU
        pos_grid *= inv_dx;
        IntVect pos_ind = IntVect(AMREX_D_DECL((int) pos_grid[0],
                                               (int) pos_grid[1],
                                               (int) pos_grid[2] ));

        // Add particle at position pos iff it's vector index is contained
        // within tile_box.
        if(tile_box.contains(pos_ind)) {

            AddSingle(lev, std::make_pair(grid_id, tile_id), radius, k_spring,
                      pos, id, cpu, i_ref);
        }
    }

    // Don't call `Redistribute` => if you want residstributing, use
    // `BlobContainer::InitSingle`
}



void BlobContainer::InitSingle(int lev, Real radius, Real k_spring,
                               const RealVect & pos, int id, int cpu, int i_ref) {

    AddSingle(lev, radius, k_spring, pos, id, cpu, i_ref);
    // We shouldn't need this if the particles are tiled with one tile per
    // grid, but otherwise we do need this to move particles from tile 0 to the
    // correct tile.
    Redistribute();
}



void BlobContainer::MovePredictor(int lev, Real dt) {

    // MovePredictor preseves the property that `pred_pos[x,y,z]` is a
    // difference to the reference position

    //___________________________________________________________________________
    // First set Set `pred_pos[x,y,z]` to the reference distances `ref_del[x,y,z]

    for (MyIBMarIter pti(* this, lev); pti.isValid(); ++pti) {

        TileIndex index(pti.index(), pti.LocalTileIndex());

        AoS & particles = this->GetParticles(lev).at(index).GetArrayOfStructs();
        long np = this->GetParticles(lev).at(index).numParticles();

        for (int i=0; i<np; ++i) {
            ParticleType & part = particles[i];

            for (int d=0; d<AMREX_SPACEDIM; ++d) {
                part.rdata(IBBReal::pred_posx + d) =
                    part.rdata(IBBReal::ref_delx + d);
            }
        }
    }

    //___________________________________________________________________________
    // Now add dt*pred_vel[x,y,z]

    IBMarkerContainerBase<IBBReal, IBBInt>::MovePredictor(lev, dt);
}



void BlobContainer::MoveMarkers(int lev, Real dt) {

    // MoveMarkers preseves the property that `ref_del[x,y,z]` is a
    // difference to the reference position


    //___________________________________________________________________________
    // First update marker position using `vel[x,y,z]

    IBMarkerContainerBase<IBBReal, IBBInt>::MoveMarkers(lev, dt);

    //___________________________________________________________________________
    // Now do the same with ref_del[x,y,z]

    for (MyIBMarIter pti(* this, lev); pti.isValid(); ++pti) {

        TileIndex index(pti.index(), pti.LocalTileIndex());

        AoS & particles = this->GetParticles(lev).at(index).GetArrayOfStructs();
        long np = this->GetParticles(lev).at(index).numParticles();

        for (int i=0; i<np; ++i) {
            ParticleType & part = particles[i];

            for (int d=0; d<AMREX_SPACEDIM; ++d) {
                part.rdata(IBBReal::ref_delx + d) +=
                    dt * part.rdata(IBBReal::velx + d);
            }
        }
    }
}



void BlobContainer::PredictorForces(int lev, Real k) {

    for (MyIBMarIter pti(* this, lev); pti.isValid(); ++pti) {

        // Get marker data (local to current thread)
        TileIndex index(pti.index(), pti.LocalTileIndex());
        AoS & markers = this->GetParticles(lev).at(index).GetArrayOfStructs();
        long np = this->GetParticles(lev).at(index).numParticles();

        // m_index.second is used to keep track of the neighbor list
        // currently we don't use the neighbor list, but we might in future
        for (MarkerListIndex m_index(0, 0); m_index.first<np; ++m_index.first) {

            ParticleType & mark = markers[m_index.first];

            for (int d=0; d<AMREX_SPACEDIM; ++d) {
                mark.rdata(IBBReal::pred_forcex + d) +=
                    - k * mark.rdata(IBBReal::pred_posx + d);
            }
        }
    }
}



void BlobContainer::PredictorForces(int lev) {

    for (MyIBMarIter pti(* this, lev); pti.isValid(); ++pti) {

        // Get marker data (local to current thread)
        TileIndex index(pti.index(), pti.LocalTileIndex());
        AoS & markers = this->GetParticles(lev).at(index).GetArrayOfStructs();
        long np = this->GetParticles(lev).at(index).numParticles();

        // m_index.second is used to keep track of the neighbor list
        // currently we don't use the neighbor list, but we might in future
        for (MarkerListIndex m_index(0, 0); m_index.first<np; ++m_index.first) {

            ParticleType & mark = markers[m_index.first];

            for (int d=0; d<AMREX_SPACEDIM; ++d) {
                mark.rdata(IBBReal::pred_forcex + d) +=
                    - mark.rdata(IBBReal::k_spring) * mark.rdata(IBBReal::pred_posx + d);
            }
        }
    }
}



void BlobContainer::MarkerForces(int lev, Real k) {

    for (MyIBMarIter pti(* this, lev); pti.isValid(); ++pti) {

        // Get marker data (local to current thread)
        TileIndex index(pti.index(), pti.LocalTileIndex());
        AoS & markers = this->GetParticles(lev).at(index).GetArrayOfStructs();
        long np = this->GetParticles(lev).at(index).numParticles();

        // m_index.second is used to keep track of the neighbor list
        // currently we don't use the neighbor list, but we might in future
        for (MarkerListIndex m_index(0, 0); m_index.first<np; ++m_index.first) {

            ParticleType & mark = markers[m_index.first];

            for (int d=0; d<AMREX_SPACEDIM; ++d) {
                mark.rdata(IBBReal::forcex + d) +=
                    - k * mark.rdata(IBBReal::ref_delx + d);
            }
        }
    }
}



void BlobContainer::MarkerForces(int lev) {

    for (MyIBMarIter pti(* this, lev); pti.isValid(); ++pti) {

        // Get marker data (local to current thread)
        TileIndex index(pti.index(), pti.LocalTileIndex());
        AoS & markers = this->GetParticles(lev).at(index).GetArrayOfStructs();
        long np = this->GetParticles(lev).at(index).numParticles();

        // m_index.second is used to keep track of the neighbor list
        // currently we don't use the neighbor list, but we might in future
        for (MarkerListIndex m_index(0, 0); m_index.first<np; ++m_index.first) {

            ParticleType & mark = markers[m_index.first];

            for (int d=0; d<AMREX_SPACEDIM; ++d) {
                mark.rdata(IBBReal::forcex + d) +=
                    - mark.rdata(IBBReal::k_spring) * mark.rdata(IBBReal::ref_delx + d);
            }
        }
    }
}




bool IBMultiBlobContainer::use_neighbor_list  {true};
bool IBMultiBlobContainer::sort_neighbor_list {false};



IBMultiBlobContainer::IBMultiBlobContainer(const Geometry & geom,
                                           const DistributionMapping & dmap,
                                           const BoxArray & ba,
                                           int n_nbhd, int blob_nbhd)
    : NeighborParticleContainer<IBMBReal::count, IBMBInt::count>(
            geom, dmap, ba, n_nbhd
        ),
      nghost(n_nbhd),
      markers(geom, dmap, ba, blob_nbhd)
{
    InitInternals(n_nbhd);
}



IBMultiBlobContainer::IBMultiBlobContainer(AmrCore * amr_core, int n_nbhd,
                                           int blob_nbhd)
    : NeighborParticleContainer<IBMBReal::count, IBMBInt::count>(
            amr_core->GetParGDB(), n_nbhd
        ),
      m_amr_core(amr_core),
      nghost(n_nbhd),
      markers(amr_core, blob_nbhd)

{
    InitInternals(n_nbhd);
}



void IBMultiBlobContainer::InitSingle(int lev, const RealVect & pos, Real r,
                                      Real rho, int n_marker, Real k_sping) {

    // Inverse cell-size vector => used for determining index corresponding to
    // IBParticle position (pos)
    RealVect inv_dx = RealVect(AMREX_D_DECL(Geom(lev).InvCellSize(0),
                                            Geom(lev).InvCellSize(1),
                                            Geom(lev).InvCellSize(2)  ));


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

        // IntVect representing particle's position in the tile_box grid.
        RealVect pos_grid = pos;
        pos_grid *= inv_dx;
        IntVect pos_ind = IntVect(AMREX_D_DECL((int) pos_grid[0],
                                               (int) pos_grid[1],
                                               (int) pos_grid[2]  ));

        // Add particle at position pos iff it's vector index is contained
        // within tile_box.
        if(tile_box.contains(pos_ind)) {
            pcount ++;

            ParticleType p_new;

            // Set id and cpu for this multiblob
            p_new.id()  = ParticleType::NextID();
            p_new.cpu() = ParallelDescriptor::MyProc();

            for (int d=0; d<AMREX_SPACEDIM; ++d){
                // Set multiblob position
                p_new.pos(d) = pos[d];

                // Initialize particle velocity (and angular velocity) as well
                // as drag to 0
                p_new.rdata(IBMBReal::velx + d)   = 0.;
                p_new.rdata(IBMBReal::omegax + d) = 0.;
                p_new.rdata(IBMBReal::dragx + d)  = 0.;

                // Initialize external forcing
                p_new.rdata(IBMBReal::pred_forcex + d) = 0.;
                p_new.rdata(IBMBReal::forcex + d)      = 0.;
            }

            // Physical radius of multiblob
            p_new.rdata(IBMBReal::radius)   = r;
            p_new.rdata(IBMBReal::rho)      = rho;
            p_new.rdata(IBMBReal::k_spring) = k_sping;

            // Number of markers to put on surface
            p_new.idata(IBMBInt::n_marker) = n_marker;
            // TODO: Audit
            p_new.idata(IBMBInt::phase) = -1;
            p_new.idata(IBMBInt::state) = -1;

            // Add to the data structure
            particles.push_back(p_new);
        }

        const int np = pcount;
        total_np += np;
    }

    // ParallelDescriptor::ReduceIntSum(total_np, ParallelDescriptor::IOProcessorNumber());
    // amrex::Print() << "Total number of generated particles: " << total_np << std::endl;

    // We shouldn't need this if the particles are tiled with one tile per
    // grid, but otherwise we do need this to move particles from tile 0 to the
    // correct tile.
    Redistribute();
}



void IBMultiBlobContainer::FillMarkerPositions(int lev) {


    // fillNeighbors();


    /****************************************************************************
     *                                                                          *
     * Fill Markert list                                                        *
     *                                                                          *
     ***************************************************************************/

    int total_np = 0;

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (IBMBIter pti(* this, lev); pti.isValid(); ++pti) {

        TileIndex index(pti.index(), pti.LocalTileIndex());

        AoS & particles = GetParticles(lev).at(index).GetArrayOfStructs();
        long np = GetParticles(lev).at(index).numParticles();

        for (int i = 0; i < np; ++i) {
            ParticleType & part = particles[i];
            MarkerIndex pindex(part.id(), part.cpu());

            //___________________________________________________________________
            // Based on the paper:
            // >*Distributing many points on a sphere*, E. B. Saff, A. B. J.
            // Kuijlaars, *The Mathematical Intelligencer*, **19** (1997)

            // elt = (particle ID, particle data)
            //        ^^ first ^^, ^^ second  ^^

            // HACK: put markers slightly inside TODO: fix
            Real     r        = part.rdata(IBMBReal::radius)*0.8;
            Real     k_spring = part.rdata(IBMBReal::k_spring);
            RealVect pos_0    = {AMREX_D_DECL(part.pos(0),
                                              part.pos(1),
                                              part.pos(2))};


            //___________________________________________________________________
            // Fill marker using Saff spiral
            int n_marker      = part.idata(IBMBInt::n_marker);
            double inv_sqrt_n = 1./std::sqrt(n_marker);
            double phi = 0.;
            for (int i=0; i<n_marker; ++i) {

                // Compute polar coordinates of marker positions
                double ck    = -1. + (2.*i)/(n_marker-1);
                double theta = std::acos(ck);

                if ( (i==0) || (i==n_marker-1) ) phi = 0;
                else phi = std::fmod(phi + 3.6*inv_sqrt_n/std::sqrt(1-ck*ck), 2*M_PI);

                // Convert to cartesian coordinates
                RealVect pos;
#if   (AMREX_SPACEDIM == 2)
                pos[0] = pos_0[0] + r*std::sin(theta);
                pos[1] = pos_0[1] + r*std::cos(theta);
#elif (AMREX_SPACEDIM == 3)
                pos[0] = pos_0[0] + r*std::sin(theta)*std::cos(phi);
                pos[1] = pos_0[1] + r*std::sin(theta)*std::sin(phi);
                pos[2] = pos_0[2] + r*std::cos(theta);
#endif

#ifdef _OPENMP
#pragma omp critical
#endif
                {   // Add to list (use the `BlobContainer::AddSingle` function
                    // and call `BlobContainer::Redistribute` **outside** the
                    // `IBMBIter` loop)
                    markers.AddSingle(lev, index, 1., k_spring, pos, part.id(),
                                      part.cpu(), i);
                }
            }

            total_np += n_marker;
        }
    }


    ParallelDescriptor::ReduceIntSum(total_np,
                                     ParallelDescriptor::IOProcessorNumber());

    amrex::Print() << "Total number of generated markers: "
                   << total_np << std::endl;

    //___________________________________________________________________________
    // Redistribute markers to correct tiles
    markers.Redistribute();
}



void IBMultiBlobContainer::ResetPredictor(int lev) {
    markers.ResetPredictor(lev);
}



void IBMultiBlobContainer::ResetMarkers(int lev) {
    markers.ResetMarkers(lev);
}



void IBMultiBlobContainer::ResetDrag(int lev) {

    for (IBMBIter pti(* this, lev); pti.isValid(); ++pti) {

        TileIndex index(pti.index(), pti.LocalTileIndex());
        auto & particle_data = this->GetParticles(lev).at(index);
        long np = this->GetParticles(lev).at(index).numParticles();

        AoS & particles = particle_data.GetArrayOfStructs();
        for (int i = 0; i < np; ++i) {
            ParticleType & part = particles[i];

            for (int d=0; d<AMREX_SPACEDIM; ++d)
                part.rdata(IBMBReal::dragx + d) = 0.;
        }
    }
}



void IBMultiBlobContainer::SpreadMarkers(int lev,
                                         std::array<MultiFab, AMREX_SPACEDIM> & f_out) const {

    // Since the Blob-Container already contains all the markers, we hand this
    // call off to the `BlobContainer markers` member.
    markers.SpreadMarkers(lev, f_out);
}



void IBMultiBlobContainer::SpreadPredictor(int lev,
                                           std::array<MultiFab, AMREX_SPACEDIM> & f_out) const {

    // Since the Blob-Container already contains all the markers, we hand this
    // call off to the `BlobContainer markers` member.
    markers.SpreadPredictor(lev, f_out);
}



void IBMultiBlobContainer::InterpolateMarkers(int lev,
                                              const std::array<MultiFab, AMREX_SPACEDIM> & f_in) {

    // Since the Blob-Container already contains all the markers, we hand this
    // call off to the `BlobContainer markers` member.
    markers.InterpolateMarkers(lev, f_in);
}



void IBMultiBlobContainer::InterpolatePredictor(int lev,
                                                const std::array<MultiFab, AMREX_SPACEDIM> & f_in) {

    // Since the Blob-Container already contains all the markers, we hand this
    // call off to the `BlobContainer markers` member.
    markers.InterpolatePredictor(lev, f_in);
}



void IBMultiBlobContainer::MoveMarkers(int lev, Real dt) {

    //___________________________________________________________________________
    // First update marker position using `vel[x,y,z]

    markers.MoveMarkers(lev, dt);


    //___________________________________________________________________________
    // Now update ref_del[x,y,z] according to blob COM velocity

    for (BlobIter pti(markers, lev); pti.isValid(); ++pti) {

        // Get marker data (local to current thread)
        TileIndex index(pti.index(), pti.LocalTileIndex());
        std::map<MarkerIndex, ParticleType *> particle_dict = GetParticleDict(lev, index);

        BlobContainer::AoS & marker_data =
            markers.GetParticles(lev).at(index).GetArrayOfStructs();
        long np = markers.GetParticles(lev).at(index).numParticles();

        // m_index.second is used to keep track of the neighbor list
        // currently we don't use the neighbor list, but we might in future
        for (MarkerListIndex m_index(0, 0); m_index.first<np; ++m_index.first) {

            BlobContainer::ParticleType & mark = marker_data[m_index.first];
            MarkerIndex parent = std::make_pair(mark.idata(IBBInt::id_0),
                                                mark.idata(IBBInt::cpu_0));

            ParticleType * blob = particle_dict.at(parent);
            for (int d=0; d<AMREX_SPACEDIM; ++d) {
                mark.rdata(IBBReal::ref_delx + d) -=
                    dt * blob->rdata(IBMBReal::velx + d);
            }
        }
    }
}



void IBMultiBlobContainer::MovePredictor(int lev, Real dt) {

    //___________________________________________________________________________
    // First update marker predictor position using `vel[x,y,z]

    markers.MovePredictor(lev, dt);


    //___________________________________________________________________________
    // Now update pred_pos[x,y,z] according to blob COM velocity

    for (BlobIter pti(markers, lev); pti.isValid(); ++pti) {

        // Get marker data (local to current thread)
        TileIndex index(pti.index(), pti.LocalTileIndex());
        std::map<MarkerIndex, ParticleType *> particle_dict = GetParticleDict(lev, index);

        BlobContainer::AoS & marker_data =
            markers.GetParticles(lev).at(index).GetArrayOfStructs();
        long np = markers.GetParticles(lev).at(index).numParticles();

        // m_index.second is used to keep track of the neighbor list
        // currently we don't use the neighbor list, but we might in future
        for (MarkerListIndex m_index(0, 0); m_index.first<np; ++m_index.first) {

            BlobContainer::ParticleType & mark = marker_data[m_index.first];
            MarkerIndex parent = std::make_pair(mark.idata(IBBInt::id_0),
                                                mark.idata(IBBInt::cpu_0));

            ParticleType * blob = particle_dict.at(parent);
            for (int d=0; d<AMREX_SPACEDIM; ++d) {
                mark.rdata(IBBReal::pred_posx + d) -=
                    dt * blob->rdata(IBMBReal::velx + d);
            }
        }
    }
}



void IBMultiBlobContainer::RedistributeMarkers() {
    markers.Redistribute();
}




void IBMultiBlobContainer::PredictorForces(int lev) {

    markers.PredictorForces(lev);

    for (BlobIter pti(markers, lev); pti.isValid(); ++pti) {

        TileIndex index(pti.index(), pti.LocalTileIndex());
        std::map<MarkerIndex, ParticleType *> particle_dict = GetParticleDict(lev, index);

        // Get marker data (local to current thread)
        BlobContainer::AoS & marker_data =
            markers.GetParticles(lev).at(index).GetArrayOfStructs();
        long np = markers.GetParticles(lev).at(index).numParticles();

        // m_index.second is used to keep track of the neighbor list
        // currently we don't use the neighbor list, but we might in future
        for (MarkerListIndex m_index(0, 0); m_index.first<np; ++m_index.first) {

            BlobContainer::ParticleType & mark = marker_data[m_index.first];
            MarkerIndex parent = std::make_pair(mark.idata(IBBInt::id_0),
                                                mark.idata(IBBInt::cpu_0));

            ParticleType * blob = particle_dict.at(parent);
            for (int d=0; d<AMREX_SPACEDIM; ++d) {
                mark.rdata(IBBReal::pred_forcex + d) +=
                    blob->rdata(IBMBReal::pred_forcex + d)/blob->idata(IBMBInt::n_marker);
            }
        }
    }
}



void IBMultiBlobContainer::PredictorForces(int lev, Real k) {

    markers.PredictorForces(lev, k);

    for (BlobIter pti(markers, lev); pti.isValid(); ++pti) {

        TileIndex index(pti.index(), pti.LocalTileIndex());
        std::map<MarkerIndex, ParticleType *> particle_dict = GetParticleDict(lev, index);

        // Get marker data (local to current thread)
        BlobContainer::AoS & marker_data =
            markers.GetParticles(lev).at(index).GetArrayOfStructs();
        long np = markers.GetParticles(lev).at(index).numParticles();

        // m_index.second is used to keep track of the neighbor list
        // currently we don't use the neighbor list, but we might in future
        for (MarkerListIndex m_index(0, 0); m_index.first<np; ++m_index.first) {

            BlobContainer::ParticleType & mark = marker_data[m_index.first];
            MarkerIndex parent = std::make_pair(mark.idata(IBBInt::id_0),
                                                mark.idata(IBBInt::cpu_0));

            ParticleType * blob = particle_dict.at(parent);
            for (int d=0; d<AMREX_SPACEDIM; ++d) {
                mark.rdata(IBBReal::pred_forcex + d) +=
                    blob->rdata(IBMBReal::pred_forcex + d)/blob->idata(IBMBInt::n_marker);
            }
        }
    }

}



void IBMultiBlobContainer::MarkerForces(int lev, Real k) {

    markers.MarkerForces(lev, k);

    for (BlobIter pti(markers, lev); pti.isValid(); ++pti) {

        TileIndex index(pti.index(), pti.LocalTileIndex());
        std::map<MarkerIndex, ParticleType *> particle_dict = GetParticleDict(lev, index);

        // Get marker data (local to current thread)
        BlobContainer::AoS & marker_data =
            markers.GetParticles(lev).at(index).GetArrayOfStructs();
        long np = markers.GetParticles(lev).at(index).numParticles();

        // m_index.second is used to keep track of the neighbor list
        // currently we don't use the neighbor list, but we might in future
        for (MarkerListIndex m_index(0, 0); m_index.first<np; ++m_index.first) {

            BlobContainer::ParticleType & mark = marker_data[m_index.first];
            MarkerIndex parent = std::make_pair(mark.idata(IBBInt::id_0),
                                                mark.idata(IBBInt::cpu_0));

            ParticleType * blob = particle_dict.at(parent);
            for (int d=0; d<AMREX_SPACEDIM; ++d) {
                mark.rdata(IBBReal::forcex + d) +=
                    blob->rdata(IBMBReal::forcex + d)/blob->idata(IBMBInt::n_marker);
            }
        }
    }
}



void IBMultiBlobContainer::MarkerForces(int lev) {

    markers.MarkerForces(lev);

    for (BlobIter pti(markers, lev); pti.isValid(); ++pti) {

        TileIndex index(pti.index(), pti.LocalTileIndex());
        std::map<MarkerIndex, ParticleType *> particle_dict = GetParticleDict(lev, index);

        // Get marker data (local to current thread)
        BlobContainer::AoS & marker_data =
            markers.GetParticles(lev).at(index).GetArrayOfStructs();
        long np = markers.GetParticles(lev).at(index).numParticles();

        // m_index.second is used to keep track of the neighbor list
        // currently we don't use the neighbor list, but we might in future
        for (MarkerListIndex m_index(0, 0); m_index.first<np; ++m_index.first) {

            BlobContainer::ParticleType & mark = marker_data[m_index.first];
            MarkerIndex parent = std::make_pair(mark.idata(IBBInt::id_0),
                                                mark.idata(IBBInt::cpu_0));

            ParticleType * blob = particle_dict.at(parent);
            for (int d=0; d<AMREX_SPACEDIM; ++d) {
                mark.rdata(IBBReal::forcex + d) +=
                    blob->rdata(IBMBReal::forcex + d)/blob->idata(IBMBInt::n_marker);
            }
        }
    }
}



std::map<IBMultiBlobContainer::MarkerIndex, IBMultiBlobContainer::ParticleType *>
IBMultiBlobContainer::GetParticleDict(int lev, const TileIndex & index) {

    std::map<MarkerIndex, ParticleType *> particle_dict;

    ParticleVector & nbhd_data = GetNeighbors(lev, index.first, index.second).GetArrayOfStructs()();
    long nn = nbhd_data.size();
    for (int j=0; j<nn; ++j) {

        ParticleType & part = nbhd_data[j];

        MarkerIndex parent = std::make_pair(part.id(), part.cpu());

        // check if already in ParticleDict NOTE: c++20 has contains()
        auto search = particle_dict.find(parent);
        particle_dict[parent] = & part;
    }


    auto & particle_data = this->GetParticles(lev).at(index);
    long np = this->GetParticles(lev).at(index).numParticles();

    AoS & particles = particle_data.GetArrayOfStructs();
    for (int i = 0; i < np; ++i) {
        ParticleType & part = particles[i];

        MarkerIndex parent = std::make_pair(part.id(), part.cpu());

        // Overwrite neighbor data with "real" particle data => don't check
        // if (ID, CPU) is already in dict
        particle_dict[parent] = & part;
    }

    return particle_dict;
}




void IBMultiBlobContainer::AccumulateDrag(int lev) {

    int total_np = 0;
    Real max_del = 0;

    for (BlobIter pti(markers, lev); pti.isValid(); ++pti) {

        TileIndex index(pti.index(), pti.LocalTileIndex());
        std::map<MarkerIndex, ParticleType *> particle_dict = GetParticleDict(lev, index);

        // Get marker data (local to current thread)
        BlobContainer::AoS & marker_data =
            markers.GetParticles(lev).at(index).GetArrayOfStructs();
        long np = markers.GetParticles(lev).at(index).numParticles();

        // m_index.second is used to keep track of the neighbor list
        // currently we don't use the neighbor list, but we might in future
        for (MarkerListIndex m_index(0, 0); m_index.first<np; ++m_index.first) {

            BlobContainer::ParticleType & mark = marker_data[m_index.first];
            MarkerIndex parent = std::make_pair(mark.idata(IBBInt::id_0),
                                                mark.idata(IBBInt::cpu_0));

            ParticleType * target = particle_dict.at(parent);
            Real mag_del = 0;
            for (int d=0; d<AMREX_SPACEDIM; ++d) {
                target->rdata(IBMBReal::dragx + d) -=
                    mark.rdata(IBBReal::forcex + d);

                mag_del += mark.rdata(IBBReal::ref_delx + d)*mark.rdata(IBBReal::ref_delx + d);
            }

            mag_del = std::sqrt(mag_del);
            max_del = amrex::max(max_del, mag_del);
        }

        total_np += np;
    }


    ParallelDescriptor::ReduceIntSum(total_np,
                                     ParallelDescriptor::IOProcessorNumber());

    amrex::Print() << "Total number of markers contributing to drag: "
                   << total_np << std::endl;

    ParallelDescriptor::ReduceRealSum(max_del,
                                     ParallelDescriptor::IOProcessorNumber());

    amrex::Print() << "Maximum |ref_del| for markers contributing to drag: "
                   << max_del << std::endl;
}



void IBMultiBlobContainer::MoveBlob(int lev, Real dt) {

    for (IBMBIter pti(* this, lev); pti.isValid(); ++pti) {

        TileIndex index(pti.index(), pti.LocalTileIndex());
        std::map<MarkerIndex, ParticleType *> particle_dict = GetParticleDict(lev, index);

        AoS & particles = this->GetParticles(lev).at(index).GetArrayOfStructs();
        long np = this->GetParticles(lev).at(index).numParticles();

        for (int i = 0; i < np; ++ i) {
            ParticleType & part = particles[i];

            for (int d=0; d<AMREX_SPACEDIM; ++d) {
                part.rdata(IBMBReal::velx + d) +=
                    dt *( part.rdata(IBMBReal::dragx + d) +
                          part.rdata(IBMBReal::forcex + d)  );
                part.pos(d) += dt * part.rdata(IBMBReal::velx + d);
            }
        }
    }


    for (BlobIter pti(markers, lev); pti.isValid(); ++pti) {

        TileIndex index(pti.index(), pti.LocalTileIndex());
        std::map<MarkerIndex, ParticleType *> particle_dict = GetParticleDict(lev, index);

        BlobContainer::AoS & marker_data =
            markers.GetParticles(lev).at(index).GetArrayOfStructs();
        long np = markers.GetParticles(lev).at(index).numParticles();

        // m_index.second is used to keep track of the neighbor list
        // currently we don't use the neighbor list, but we might in future
        for (MarkerListIndex m_index(0, 0); m_index.first<np; ++m_index.first) {

            BlobContainer::ParticleType & mark = marker_data[m_index.first];
            MarkerIndex parent = std::make_pair(mark.idata(IBBInt::id_0),
                                                mark.idata(IBBInt::cpu_0));

            ParticleType * blob = particle_dict.at(parent);
            for (int d=0; d<AMREX_SPACEDIM; ++d) {
                mark.rdata(IBBReal::ref_delx + d) -=
                    dt*blob->rdata(IBMBReal::velx + d);
            }
        }
    }
}



void IBMultiBlobContainer::WritePlotFile(const std::string & dir,
                                         const std::string & name,
                                         const std::string & blob_name) const {

    // save multi-blob data
    NeighborParticleContainer<IBMBReal::count, IBMBInt::count>::WritePlotFile(
                dir, name, IBMBReal::names(), IBMBInt::names()
            );

    // save marker data
    markers.WritePlotFile(dir, blob_name, IBBReal::names(), IBBInt::names());
}



void IBMultiBlobContainer::PrintMarkerData(int lev) const {

    // Inverse cell-size vector => max is used for determining IBParticle
    // radius in units of cell size
    Vector<Real> inv_dx = {AMREX_D_DECL(this->Geom(lev).InvCellSize(0),
                                        this->Geom(lev).InvCellSize(1),
                                        this->Geom(lev).InvCellSize(2)   )};

    // Find max inv_dx (in case we have an anisotropic grid)
    Real mx_inv_dx = * std::max_element(inv_dx.begin(), inv_dx.end());


    amrex::AllPrintToFile("ib_marker_data") << "Particles on each box:" << std::endl;


    long local_count = 0;

    // ParIter skips tiles without particles => Iterate over MultiFab instead
    // of ParticleIter. Note also that AmrexParticleContainer uses strange
    // tiling => don't turn off tiling (particles are stored in tile)
    for(MFIter pti = this->MakeMFIter(lev, true); pti.isValid(); ++pti) {
        // MuliFabs are indexed using a pair: (BoxArray index, tile index):
        TileIndex index(pti.index(), pti.LocalTileIndex());

        // Neighbours are stored as raw data (see below)
        int ng = this->neighbors[lev].at(index).size();

        auto & particle_data = this->GetParticles(lev).at(index);
        long np = this->GetParticles(lev).at(index).numParticles();

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
        const ParticleType * nbhd_data = (ParticleType *) this->neighbors[lev].at(index).GetArrayOfStructs().dataPtr();
        for(int i = 0; i < ng; i++){
            const ParticleType & part = nbhd_data[i];

            AllPrintToFile("ib_marker_data") << "   +--> " << part << std::endl;
        }
    }

    AllPrintToFile("ib_marker_data") << "Total for this process: "
                                     << local_count << std::endl << std::endl;
}



void IBMultiBlobContainer::InitInternals(int ngrow) {

    ReadStaticParameters();

    this->SetVerbose(0);

    // Turn off certain components for ghost particle communication
    // Field numbers: {0, 1, 2} => {x, y, z} particle coordinates
    //      => 3 corresponds to the start of IBM_realData
    for (int i=3; i < IBMBReal::count + 3; ++i)
        this->setRealCommComp(i,  true);

    // Field numbers: {0, 1} => {ID, CPU}
    //      => 2 corresponds to the start of IBM_intData
    // We _do_ want the the neighbour particles to have ID and cpu init data.
    for (int i = 2; i < IBMBInt::count + 2; ++i)
        this->setIntCommComp(i, true);

    // Needed to copy force data back to owner
    this->setEnableInverse(true);
}



void IBMultiBlobContainer::ReadStaticParameters() {

    static bool initialized = false;

    if (!initialized) {
        ParmParse pp("particles");

//        // AMReX default is false => enable by default
//        this->do_tiling = true;
//        // Allow user to overwrite
//        pp.query("do_tiling",  this->do_tiling);

        // If tiling is enabled, make sure that the tile size is at least the
        // number of ghost cells (otherwise strange things happen)
        if (this->do_tiling)
            this->tile_size = IntVect{AMREX_D_DECL(max_grid_size[0],
                                                   max_grid_size[1],
                                                   max_grid_size[2])};
            //this->tile_size = IntVect{AMREX_D_DECL(nghost, nghost, nghost)};

        // User can overwrite
        Vector<int> ts(BL_SPACEDIM);
        if (pp.queryarr("tile_size", ts))
            this->tile_size = IntVect(ts);

        // pp.query("use_neighbor_list", use_neighbor_list);
        // pp.query("sort_neighbor_list", sort_neighbor_list);

        initialized = true;
    }
}
