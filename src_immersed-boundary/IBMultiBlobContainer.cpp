#include <AMReX.H>
#include <AMReX_AmrCore.H>
#include <AMReX_AmrParGDB.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Particles.H>
#include <AMReX_NeighborParticles.H>

#include <AMReX_VisMF.H>  // amrex::VisMF::Write(MultiFab)

#include <common_functions.H>
#include <common_namespace.H>

#include <IBMultiBlobContainer.H>
#include <ib_functions_F.H>



using namespace common;
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



void BlobContainer::AddSingle(int lev,
                              Real radius, const RealVect & pos,
                              int id, int cpu, int i_ref ) {

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

            // 2. Blob contexual metadata
            p_new.idata(IBBInt::id_0)  = id;
            p_new.idata(IBBInt::cpu_0) = cpu;

            p_new.idata(IBBInt::id_1)  = i_ref;
            p_new.idata(IBBInt::cpu_1) = -1;

            // Add to the data structure
            particles.push_back(p_new);
        }
    }

    // Don't call `Redistribute` => if you want residstributing, use
    // `BlobContainer::InitSingle`
}



void BlobContainer::InitSingle(int lev,
                               Real radius, const RealVect & pos,
                               int id, int cpu, int i_ref ) {

    AddSingle(lev, radius, pos, id, cpu, i_ref);
    // We shouldn't need this if the particles are tiled with one tile per
    // grid, but otherwise we do need this to move particles from tile 0 to the
    // correct tile.
    Redistribute();
}



void BlobContainer::MoveMarkers(int lev, Real dt) {

    IBMarkerContainerBase<IBBReal, IBBInt>::MoveMarkers(lev, dt);

    for (MyIBMarIter pti(* this, lev); pti.isValid(); ++pti) {

        TileIndex index(pti.index(), pti.LocalTileIndex());

        AoS & particles = this->GetParticles(lev).at(index).GetArrayOfStructs();
        long np = particles.size();

        for (int i=0; i<np; ++i) {
            ParticleType & part = particles[i];

            for (int d=0; d<AMREX_SPACEDIM; ++d) {
                part.rdata(IBBReal::ref_delx + d) =
                    dt * part.rdata(IBBReal::velx + d);
            }
        }
    }
}



bool IBMultiBlobContainer::use_neighbor_list  {true};
bool IBMultiBlobContainer::sort_neighbor_list {false};



IBMultiBlobContainer::IBMultiBlobContainer(const Geometry & geom,
                                           const DistributionMapping & dmap,
                                           const BoxArray & ba,
                                           int n_nbhd)
    : NeighborParticleContainer<IBMB_realData::count, IBMB_intData::count>(
            geom, dmap, ba, n_nbhd
        ),
      nghost(n_nbhd),
      markers(geom, dmap, ba, n_nbhd)
{
    InitInternals(n_nbhd);
}



IBMultiBlobContainer::IBMultiBlobContainer(AmrCore * amr_core, int n_nbhd)
    : NeighborParticleContainer<IBMB_realData::count, IBMB_intData::count>(
            amr_core->GetParGDB(), n_nbhd
        ),
      m_amr_core(amr_core),
      nghost(n_nbhd),
      markers(amr_core, n_nbhd)

{
    InitInternals(n_nbhd);
}



void IBMultiBlobContainer::InitSingle(int lev, const RealVect & pos, Real r, Real rho) {

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
                p_new.rdata(IBMB_realData::velx + d)   = 0.;
                p_new.rdata(IBMB_realData::omegax + d) = 0.;
                p_new.rdata(IBMB_realData::dragx + d)  = 0.;
            }

            // Physical radius of multiblob
            p_new.rdata(IBMB_realData::radius) = r;

            // TODO: Audit
            p_new.idata(IBMB_intData::phase) = -1;
            p_new.idata(IBMB_intData::state) = -1;

            // Add to the data structure
            particles.push_back(p_new);
        }

        const int np = pcount;
        total_np += np;
    }

    ParallelDescriptor::ReduceIntSum(total_np,ParallelDescriptor::IOProcessorNumber());
    amrex::Print() << "Total number of generated particles: " << total_np << std::endl;

    // We shouldn't need this if the particles are tiled with one tile per
    // grid, but otherwise we do need this to move particles from tile 0 to the
    // correct tile.
    Redistribute();
}



void IBMultiBlobContainer::FillMarkerPositions(int lev, int n_marker) {


    double inv_sqrt_n = 1./std::sqrt(n_marker);


    fillNeighbors();


    /****************************************************************************
     *                                                                          *
     * Fill Markert list                                                        *
     *                                                                          *
     ***************************************************************************/


#ifdef _OPENMP
#pragma omp parallel
#endif
    for (IBMBIter pti(* this, lev); pti.isValid(); ++pti) {

        PairIndex index(pti.index(), pti.LocalTileIndex());

        AoS & particles = GetParticles(lev).at(index).GetArrayOfStructs();
        long np = particles.size();

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
            double   r     = part.rdata(IBMB_realData::radius)*0.8;
            RealVect pos_0 = {AMREX_D_DECL(part.pos(0), part.pos(1), part.pos(2))};


            //___________________________________________________________________
            // Fill marker using Saff spiral
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
                    markers.AddSingle(lev, 1., pos, part.id(), part.cpu(), i);
                }
            }
        }
    }
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
    markers.MoveMarkers(lev, dt);
}



void IBMultiBlobContainer::MovePredictor(int lev, Real dt) {
    markers.MovePredictor(lev, dt);
}



void IBMultiBlobContainer::InitInternals(int ngrow) {
    ReadStaticParameters();

    this->SetVerbose(0);

    // Turn off certain components for ghost particle communication
    // Field numbers: {0, 1, 2} => {x, y, z} particle coordinates
    //      => 3 corresponds to the start of IBP_realData
    // setRealCommComp(4, false);   // IBP_realData.volume

    // Field numbers: {0, 1} => {ID, CPU}
    //      => 2 corresponds to the start of IBP_intData
    // We _do_ want the the neighbour particles to have ID and cpu init data.
    // setIntCommComp(0, false);  // IBP_intData.phase
    // setIntCommComp(1, false);  // IBP_intData.state
}



void IBMultiBlobContainer::ReadStaticParameters() {
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

