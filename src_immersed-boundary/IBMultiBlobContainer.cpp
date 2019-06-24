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

bool IBMultiBlobContainer::use_neighbor_list  {true};
bool IBMultiBlobContainer::sort_neighbor_list {false};



IBMultiBlobContainer::IBMultiBlobContainer(const Geometry & geom,
                                           const DistributionMapping & dmap,
                                           const BoxArray & ba,
                                           int n_nbhd)
    : NeighborParticleContainer<IBMB_realData::count, IBMB_intData::count>(
            geom, dmap, ba, n_nbhd
        ),
    nghost(n_nbhd)
{
    InitInternals(n_nbhd);
}



IBMultiBlobContainer::IBMultiBlobContainer(AmrCore * amr_core, int n_nbhd)
    : NeighborParticleContainer<IBMB_realData::count, IBMB_intData::count>(
            amr_core->GetParGDB(), n_nbhd
        ),
    m_amr_core(amr_core),
    nghost(n_nbhd)
{
    InitInternals(n_nbhd);
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

