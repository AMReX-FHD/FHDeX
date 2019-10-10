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



void BlobContainer::InitSingle(int lev,
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

            // Set particle position
            p_new.pos(0) = pos[0];
            p_new.pos(1) = pos[1];
            p_new.pos(2) = pos[2];

            p_new.rdata(IBBReal::radius) = radius;

            // Initialize marker velocity as well as forces to 0
            p_new.rdata(IBBReal::velx)   = 0.;
            p_new.rdata(IBBReal::vely)   = 0.;
            p_new.rdata(IBBReal::velz)   = 0.;

            p_new.rdata(IBBReal::forcex) = 0.;
            p_new.rdata(IBBReal::forcey) = 0.;
            p_new.rdata(IBBReal::forcez) = 0.;

            p_new.rdata(IBBReal::pred_posx)   = 0.;
            p_new.rdata(IBBReal::pred_posy)   = 0.;
            p_new.rdata(IBBReal::pred_posz)   = 0.;

            p_new.rdata(IBBReal::pred_velx)   = 0.;
            p_new.rdata(IBBReal::pred_vely)   = 0.;
            p_new.rdata(IBBReal::pred_velz)   = 0.;

            p_new.rdata(IBBReal::pred_forcex) = 0.;
            p_new.rdata(IBBReal::pred_forcey) = 0.;
            p_new.rdata(IBBReal::pred_forcez) = 0.;

            p_new.idata(IBBInt::id_0)  = id;
            p_new.idata(IBBInt::cpu_0) = cpu;

            p_new.idata(IBBInt::id_1)  = i_ref;
            p_new.idata(IBBInt::cpu_1) = -1;

            // Add to the data structure
            particles.push_back(p_new);
        }
    }

    // We shouldn't need this if the particles are tiled with one tile per
    // grid, but otherwise we do need this to move particles from tile 0 to the
    // correct tile.
    Redistribute();

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



void IBMultiBlobContainer::InitList(int lev,
                                   const Vector<RealVect> & pos,
                                   const Vector<Real> & r,
                                   const Vector<Real> & rho) {

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

        for(int i = 0; i < pos.size(); i++) {
            // IntVect representing particle's position in the tile_box grid.
            RealVect pos_grid = pos[i];
            pos_grid *= inv_dx;
            IntVect pos_ind = IntVect(AMREX_D_DECL((int) pos_grid[0],
                                                   (int) pos_grid[1],
                                                   (int) pos_grid[2]  ));

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

                p_new.rdata(IBMB_realData::radius) = r[i];

                // Initialize particle velocity (and angular velocity) as well
                // as drag to 0
                p_new.rdata(IBMB_realData::velx)   = 0.;
                p_new.rdata(IBMB_realData::vely)   = 0.;
                p_new.rdata(IBMB_realData::velz)   = 0.;

                p_new.rdata(IBMB_realData::omegax) = 0.;
                p_new.rdata(IBMB_realData::omegay) = 0.;
                p_new.rdata(IBMB_realData::omegaz) = 0.;

                p_new.rdata(IBMB_realData::dragx)  = 0.;
                p_new.rdata(IBMB_realData::dragy)  = 0.;
                p_new.rdata(IBMB_realData::dragz)  = 0.;

                // TODO: Audit
                p_new.idata(IBMB_intData::phase) = -1;
                p_new.idata(IBMB_intData::state) = -1;

                // Add to the data structure
                particles.push_back(p_new);
            }
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

    //___________________________________________________________________________
    // Ensure that the marker lists have enough levels, and clear previous ones
    if (marker_ref_pos.size() <= lev)
        marker_ref_pos.resize(lev+1);


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
            ParticleIndex pindex(part.id(), part.cpu());

            //___________________________________________________________________
            // Based on the paper:
            // >*Distributing many points on a sphere*, E. B. Saff, A. B. J.
            // Kuijlaars, *The Mathematical Intelligencer*, **19** (1997)

            // elt = (particle ID, particle data)
            //        ^^ first ^^, ^^ second  ^^

            //___________________________________________________________________
            // Create blank marker list, and access particle data
            // ... initialized to (0..0) by default constructor
            marker_ref_pos[lev][pindex].resize(n_marker);


            double   r     = part.rdata(IBMB_realData::radius)*0.8; // HACK: put markers slightly inside
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
                {
                    // Add to list
                    marker_ref_pos[lev][pindex][i] = pos;
                    markers.InitSingle(lev, 1., pos, part.id(), part.cpu(), i);
                }
            }
        }
    }


#ifdef _OPENMP
#pragma omp parallel
#endif
    for (IBMBIter pti(* this, lev); pti.isValid(); ++pti) {

        ParticleVector & particles = GetNeighbors(lev, pti.index(),
                                                  pti.LocalTileIndex());
        long np = particles.size();

        for (int i = 0; i < np; ++i) {
            ParticleType & part = particles[i];
            ParticleIndex pindex(part.id(), part.cpu());

            //___________________________________________________________________
            // Based on the paper:
            // >*Distributing many points on a sphere*, E. B. Saff, A. B. J.
            // Kuijlaars, *The Mathematical Intelligencer*, **19** (1997)

            // elt = (particle ID, particle data)
            //        ^^ first ^^, ^^ second  ^^

            //___________________________________________________________________
            // Create blank marker list, and access particle data
            // ... initialized to (0..0) by default constructor
            marker_ref_pos[lev][pindex].resize(n_marker);


            double   r     = part.rdata(IBMB_realData::radius)*0.8; // HACK: put markers slightly inside
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
                {
                    // Add to list
                    marker_ref_pos[lev][pindex][i] = pos;
                    // Don't do this for the neighbor-markers (already being
                    // done by the owner's core)
                    // markers.InitSingle(lev, 1., pos, part.id(), part.cpu(), i);
                }
            }
        }
    }
}



void IBMultiBlobContainer::SpreadMarkers(int lev,
                                         std::array<MultiFab, AMREX_SPACEDIM> & f_out) const {


}



void IBMultiBlobContainer::SpreadPredictor(int lev,
                                           std::array<MultiFab, AMREX_SPACEDIM> & f_out) const {

}



void IBMultiBlobContainer::InterpolateMarkers(int lev,
                                              const std::array<MultiFab, AMREX_SPACEDIM> & f_in) {


}



void IBMultiBlobContainer::InterpolatePredictor(int lev,
                                                const std::array<MultiFab, AMREX_SPACEDIM> & f_in) {


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

