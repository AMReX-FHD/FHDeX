#include <AMReX.H>
#include <AMReX_AmrParGDB.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Particles.H>
#include <AMReX_NeighborParticles.H>

#include <AMReX_VisMF.H>  // amrex::VisMF::Write(MultiFab)

#include <common_functions.H>

#include <IBMarkerContainer.H>
#include <ib_functions_F.H>



using namespace amrex;

bool IBMarkerContainer::use_neighbor_list  {true};
bool IBMarkerContainer::sort_neighbor_list {false};



IBMarkerContainer::IBMarkerContainer(const Geometry & geom,
                                         const DistributionMapping & dmap,
                                         const BoxArray & ba,
                                         int n_nbhd)
    : NeighborParticleContainer<IBP_realData::count, IBP_intData::count>(
            geom, dmap, ba, n_nbhd
        ),
    nghost(n_nbhd)
{
    InitInternals(n_nbhd);
}



IBMarkerContainer::IBMarkerContainer(AmrCore * amr_core, int n_nbhd)
    : NeighborParticleContainer<IBP_realData::count, IBP_intData::count>(
            amr_core->GetParGDB(), n_nbhd
        ),
    m_amr_core(amr_core),
    nghost(n_nbhd)
{
    InitInternals(n_nbhd);
}



void IBParticleContainer::InitList(int lev,
                                   const Vector<RealVect> & pos,
                                   const Vector<Real> & r,
                                   const Vector<Real> & rho) {

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
                p_new.rdata(IBP_realData::velx)   = 0.;
                p_new.rdata(IBP_realData::vely)   = 0.;
                p_new.rdata(IBP_realData::velz)   = 0.;

                p_new.rdata(IBP_realData::forcex) = 0.;
                p_new.rdata(IBP_realData::forcey) = 0.;
                p_new.rdata(IBP_realData::forcez) = 0.;

                // TODO: Audit
                p_new.idata(IBP_intData::phase)     = 1;
                p_new.idata(IBP_intData::state)     = pstate;

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




void IBParticleContainer::SpreadMarkers(int lev, const ParticleIndex & pindex,
        const Vector<RealVect> & f_in, std::array<MultiFab, AMREX_SPACEDIM> & f_out,
        std::array<MultiFab, AMREX_SPACEDIM> & f_weights) const {


    //___________________________________________________________________________
    // Don't do anything if pindex isn't on this rank
    auto part_it = marker_positions[lev].find(pindex);
    if (part_it == marker_positions[lev].end())
        return;


    //___________________________________________________________________________
    // Geometry data
    const Geometry & geom = Geom(0);
    const Real     *   dx = geom.CellSize();

    const int n_marker = marker_positions[lev].at(pindex).size();


    //___________________________________________________________________________
    // Cell-centered MultiFab used as a reference for iterating over data
    // WARNING: this will break if IBParticleContainer is on a differnt grid
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
                       marker_positions[lev].at(pindex).dataPtr(),
                       f_in.dataPtr(),
                       & n_marker,
                       dx );
    }
}







void IBMarkerContainer::InitInternals(int ngrow) {
    ReadStaticParameters();

    this->SetVerbose(0);

    // Turn off certain components for ghost particle communication
    // Field numbers: {0, 1, 2} => {x, y, z} particle coordinates
    //      => 3 corresponds to the start of IBP_realData
    setRealCommComp(4, true);  // IBM_realData.velx
    setRealCommComp(5, true);  // IBM_realData.vely
    setRealCommComp(6, true);  // IBM_realData.velz
    setRealCommComp(7, true);  // IBM_realData.forcex
    setRealCommComp(8, true);  // IBM_realData.forcey
    setRealCommComp(9, true);  // IBM_realData.forcez

    // Field numbers: {0, 1} => {ID, CPU}
    //      => 2 corresponds to the start of IBP_intData
    // We _do_ want the the neighbour particles to have ID and cpu init data.
    //setIntCommComp(0, false);  // IBM_intData.phase
    //setIntCommComp(1, false);  // IBM_intData.state


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
