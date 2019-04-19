#include <AMReX.H>
#include <AMReX_NeighborParticles.H>

#include <AMReX_VisMF.H>


#include <IBCore.H>
#include <ib_functions_F.H>



IBCore::IBCore () {

}



IBCore::~IBCore () {

}



void IBCore::InitData () {

}



// TODO: Implement refinement. For now this does nothing.
void IBCore::ErrorEst (int lev, TagBoxArray & tags, Real time,
                       int ngrow) {

}



void IBCore::MakeNewLevelFromScratch (int lev, Real time,
                                      const BoxArray & a_ba, const DistributionMapping & a_dm) {

    // Copy target BA and DM data into internal storage
    dmap[lev]  = a_dm;
    grids[lev] = amrex::convert(a_ba,       IntVect::TheCellVector());
    ba_nd      = amrex::convert(grids[lev], IntVect::TheNodeVector());


    // Define a new level-set as well as the identifier iMultiFab and the
    // interface tags
    ls            = std::unique_ptr<MultiFab> (new  MultiFab);
    ls_id         = std::unique_ptr<iMultiFab>(new iMultiFab);
    tag_interface = std::unique_ptr<iMultiFab>(new iMultiFab);
    ls_vel        = std::unique_ptr<MultiFab> (new  MultiFab);

    const int n_pad = 1; // Pad level-set by at least 1 ghost cell =>
                         // interpolation requires at least 1 neighbor node.



   /****************************************************************************
    * Allocate new MultiFabs:                                                  *
    *       1. ls (level-set)                                                  *
    *       2. ls_id (nearest IBParticle tags)                                 *
    *       3. tag_interface (points near interface)                           *
    ****************************************************************************/

    ls->define(ba_nd, dmap[lev], 1, n_pad);
    ls->setVal(0);

    // The ls_id iMultiFab stores the identifier information of the _nearest_
    // IBParticle. The components are:  ID, CPU, INDEX
    //      unique particle indetifier --+---+     |
    //      index in local AoS --------------------+
    ls_id->define(ba_nd, dmap[lev], 3, n_pad);
    ls_id->setVal(-1);


    // Tag those cells that are exactly 1 from an interface (ls = 0)
    tag_interface->define(grids[lev], dmap[lev], 1, n_pad);
    tag_interface->setVal(0);


    ls_vel->define(ba_nd, dmap[lev], 3, n_pad);
    ls_vel->setVal(0);



    /****************************************************************************
     * Fill neighbors: IB particles can overlap from nearby grids               *
     ****************************************************************************/

    ib_pc->fillNeighbors();



   /****************************************************************************
    * Collect local Immersed-Boundary Particle Data                            *
    ****************************************************************************/

    n_ibm_loc = 0;
    part_loc.clear();
    part_index_loc.clear();


    // ParIter skips tiles without particles => Iterate over MultiFab instead
    // of ParticleIter for the following.
#ifdef _OPENMP
#pragma omp parallel reduction(+:n_ibm_loc)
#endif
    for (MFIter mfi(* ls, ib_pc->tile_size); mfi.isValid(); ++mfi) {
        // NOTE: mfi's tile size must match the ParticleContainer tile size
        // MuliFabs are indexed using a pair: (BoxArray index, tile index).


        //_______________________________________________________________________
        // Get immersed-boundary data from IBParticleContainer
        PairIndex index(mfi.index(), mfi.LocalTileIndex());

        Vector<IBP_info> info = ib_pc->IBParticleInfo(lev, index);
        n_ibm_loc = info.size();
        test_interface(info.dataPtr(), & n_ibm_loc);


        //_______________________________________________________________________
        // Collect local (to memory) unique copies of particle data. Use the
        // `critical` pragma to prevent race conditions.
#ifdef _OPENMP
#pragma omp critical
#endif
        {
            for (IBP_info pinfo : info) {
                std::pair<int,int> pindex(pinfo.id, pinfo.cpu);

                // Check if particle's index (pindex) is already contained in
                // the local list. Neighbor particles can be duplicates.
                if ( std::find(part_index_loc.begin(), part_index_loc.end(), pindex)
                        == part_index_loc.end()) {
                    part_index_loc.push_back(pindex);
                    part_loc.push_back(pinfo);
                }
            }
        }
    }

    // Make sure that the `n_ibm_loc` integer does not count duplicates
    if (n_ibm_loc > 0) n_ibm_loc = part_loc.size();



   /****************************************************************************
    * Fill Level Set and Tags                                                  *
    ****************************************************************************/

    RealVect dx = RealVect(
        AMREX_D_DECL( Geom(lev).CellSize(0), Geom(lev).CellSize(1), Geom(lev).CellSize(2) )
        );


#ifdef _OPENMP
#pragma omp critical
#endif
    for (MFIter mfi(* ls, true); mfi.isValid(); ++mfi) {
        //_______________________________________________________________________
        // Fill the global level-set data
        const Box & tile_box = mfi.tilebox();
        auto & phi_tile      = (* ls)[mfi];
        auto & tag_tile      = (* ls_id)[mfi];
        auto & vel_tile      = (* ls_vel)[mfi];

        fill_levelset_ib (BL_TO_FORTRAN_BOX(tile_box),
                          part_loc.dataPtr(), & n_ibm_loc,
                          BL_TO_FORTRAN_3D(phi_tile),
                          BL_TO_FORTRAN_3D(tag_tile),
                          BL_TO_FORTRAN_3D(vel_tile),
                          dx.dataPtr());
    }

    ls->FillBoundary(Geom(lev).periodicity());
    ls_id->FillBoundary(Geom(lev).periodicity());
    ls_vel->FillBoundary(Geom(lev).periodicity());



   /****************************************************************************
    * Tag IB particle Interfaces                                               *
    ****************************************************************************/

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(* tag_interface, true); mfi.isValid(); ++ mfi) {
        const FArrayBox & phi_tile   = (* ls)[mfi];
        // Using the ID tags to test if current box contains immersed boundaries
        const IArrayBox & tag_tile   = (* ls_id)[mfi];
              IArrayBox & iface_tile = (* tag_interface)[mfi];

        tag_interface_ib (BL_TO_FORTRAN_3D(iface_tile),
                          BL_TO_FORTRAN_3D(phi_tile),
                          BL_TO_FORTRAN_3D(tag_tile));
    }

    tag_interface->FillBoundary(Geom(lev).periodicity());




    /****************************************************************************
     * Construct Local Immersed-Boundary Data                                   *
     ***************************************************************************/

    //___________________________________________________________________________
    // Allocate MultiFabs for local levelset data
    // TODO: currently these have ib_pc->get_nghost() many ghost cells, this is
    // way too much. Once Andrew implements back-commuincation for neighbor
    // particles, we'll only need 1 ghost cell.
    if (n_ibm_loc > 0) {
        level_sets_loc.resize(n_ibm_loc);
        iface_tags_loc.resize(n_ibm_loc);
        level_set_valid.resize(n_ibm_loc);

       for (int i=0; i<n_ibm_loc; ++i) {
            level_sets_loc[i].define(ba_nd, dmap[lev], 1, ib_pc->get_nghost());
            // Tag those cells that are exactly 1 from an interface (ls = 0)
            iface_tags_loc[i].define(grids[lev], dmap[lev], 1, ib_pc->get_nghost());
            // Tag those boxes which are being touched
            level_set_valid[i].define(ba_nd, dmap[lev], 1, ib_pc->get_nghost());

            // The default is important as not every box will be touched
            // (because there might be no corresponding particle/neighbor
            // particle in this core domain). TODO: use local level-set
            // approach
            level_sets_loc[i].setVal(0.);
            iface_tags_loc[i].setVal(0);
            level_set_valid[i].setVal(0);
        }
    }


    //___________________________________________________________________________
    // Fill each level-set and interface tag MultiFab
    for (int i=0; i<n_ibm_loc; ++i) {
#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(level_sets_loc[i], true); mfi.isValid(); ++mfi) {
            const Box & tile_box = mfi.tilebox();
            auto & phi_tile      = level_sets_loc[i][mfi];

            fill_levelset_sphere (BL_TO_FORTRAN_BOX(tile_box),
                                  & part_loc[i],
                                  BL_TO_FORTRAN_3D(phi_tile),
                                  dx.dataPtr());

            // `level_sets_loc` and `iface_tags_loc` are local to each MPI rank
            // => we tag all boxes that get touched as "valid".
            level_set_valid[i].setVal(1, tile_box, 1);
        }

        level_sets_loc[i].FillBoundary(Geom(lev).periodicity());
        level_set_valid[i].FillBoundary(Geom(lev).periodicity());


#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(level_sets_loc[i], true); mfi.isValid(); ++ mfi) {
            const FArrayBox & phi_tile   = level_sets_loc[i][mfi];
            const IArrayBox & tag_tile   = level_set_valid[i][mfi];
                  IArrayBox & iface_tile = iface_tags_loc[i][mfi];

            tag_interface_ib (BL_TO_FORTRAN_3D(iface_tile),
                              BL_TO_FORTRAN_3D(phi_tile),
                              BL_TO_FORTRAN_3D(tag_tile));
        }

        iface_tags_loc[i].FillBoundary(Geom(lev).periodicity());

    }


    save_levelset_data();
}



// TODO: Implement refinement. For now this does nothing.
void IBCore::MakeNewLevelFromCoarse (int lev, Real time,
                                     const BoxArray & ba, const DistributionMapping & dm ) {

}



void IBCore::RemakeLevel (int lev, Real time,
                          const BoxArray & a_ba, const DistributionMapping & a_dm) {

    Periodicity periodicity = Geom(lev).periodicity();
    BoxArray ba_cc          = convert(a_ba, IntVect::TheCellVector());

    bool ba_changed = (grids[lev] != ba_cc);
    bool dm_changed = (dmap[lev]  != a_dm);

    if (ba_changed || dm_changed) {
        update_loc = true;

        grids[lev] = ba_cc;
        ba_nd      = convert(grids[lev], IntVect::TheNodeVector());

        dmap[lev] = a_dm;
    }


    // Only overwrite level if BA/DM have changed
    if (update_loc) {

        std::unique_ptr<MultiFab> ls_new = IBMFUtil::duplicate<MultiFab>(
                ba_nd, dmap[lev], * ls, periodicity
            );
        ls = std::move(ls_new);

        std::unique_ptr<iMultiFab> ls_id_new = IBMFUtil::duplicate<iMultiFab>(
                ba_nd, dmap[lev], * ls_id, periodicity
            );
        ls_id = std::move(ls_id_new);

        std::unique_ptr<iMultiFab> tag_interface_new = IBMFUtil::duplicate<iMultiFab>(
                grids[lev], dmap[lev], * tag_interface, periodicity
            );
        tag_interface = std::move(tag_interface_new);

        std::unique_ptr<MultiFab> ls_vel_new = IBMFUtil::duplicate<MultiFab>(
                ba_nd, dmap[lev], * ls_vel, periodicity
            );
        ls_vel = std::move(ls_vel_new);
    }

    update_loc = false;
}



// TODO: Implement refinement. For now this does nothing.
void IBCore::ClearLevel (int lev){

}



void IBCore::IBForceDeposition (       MultiFab & f_u,       MultiFab & f_v,       MultiFab & f_w,
                                       MultiFab & u_d,       MultiFab & v_d,       MultiFab & w_d,
                                 const MultiFab & u_s, const MultiFab & v_s, const MultiFab & w_s,
                                 int lev, Real dt) {

    f_u.setVal(0.); f_v.setVal(0.); f_w.setVal(0.);
    u_d.setVal(0.); v_d.setVal(0.); w_d.setVal(0.);


    const Real strttime = ParallelDescriptor::second();


    /****************************************************************************
     * Create local copies => which have the same BA/DM as f and u              *
     ***************************************************************************/

    std::unique_ptr<iMultiFab> et = std::unique_ptr<iMultiFab>(
	new iMultiFab(u_s.boxArray(), u_s.DistributionMap(), 1, u_s.nGrow())
	);

    et->setVal(0);

    //Ensure that the target BA is cell-centered:
    BoxArray ba_cc_target = convert(u_s.boxArray(), IntVect{0, 0, 0});
    RemakeLevel (lev, 0., ba_cc_target, u_s.DistributionMap());


#ifdef _OPENMP
#pragma omp parallel
#endif
    for(MFIter mfi(v_d, true); mfi.isValid(); ++ mfi) {
        const Box & tile_box = mfi.tilebox();

        FArrayBox & u_d_tile = u_d[mfi];
        FArrayBox & v_d_tile = v_d[mfi];
        FArrayBox & w_d_tile = w_d[mfi];

        auto & et_tile = (* et)[mfi];

        const FArrayBox & u_s_tile = u_s[mfi];
        const FArrayBox & v_s_tile = v_s[mfi];
        const FArrayBox & w_s_tile = w_s[mfi];

        const FArrayBox & ls_tile    = (* ls)[mfi];
        const IArrayBox & iface_tile = (* tag_interface)[mfi];
        const FArrayBox & vel_tile   = (* ls_vel)[mfi];

        interpolate_ib_staggered (BL_TO_FORTRAN_BOX(tile_box),
                                  BL_TO_FORTRAN_3D(u_d_tile),
                                  BL_TO_FORTRAN_3D(v_d_tile),
                                  BL_TO_FORTRAN_3D(w_d_tile),
                                  BL_TO_FORTRAN_3D(u_s_tile),
                                  BL_TO_FORTRAN_3D(v_s_tile),
                                  BL_TO_FORTRAN_3D(w_s_tile),
                                  BL_TO_FORTRAN_3D(et_tile),
                                  BL_TO_FORTRAN_3D(ls_tile),
                                  BL_TO_FORTRAN_3D(iface_tile),
                                  BL_TO_FORTRAN_3D(vel_tile));
    }


#ifdef _OPENMP
#pragma omp parallel
#endif
    for(MFIter mfi(f_u, true); mfi.isValid(); ++ mfi) {
        const Box & tile_box = mfi.tilebox();

        FArrayBox & f_u_tile = f_u[mfi];
        FArrayBox & f_v_tile = f_v[mfi];
        FArrayBox & f_w_tile = f_w[mfi];

        const auto & et_tile = (* et)[mfi];

        const FArrayBox & u_s_tile = u_s[mfi];
        const FArrayBox & v_s_tile = v_s[mfi];
        const FArrayBox & w_s_tile = w_s[mfi];

        const FArrayBox & u_d_tile = u_d[mfi];
        const FArrayBox & v_d_tile = v_d[mfi];
        const FArrayBox & w_d_tile = w_d[mfi];


        fill_force_ib_staggered (BL_TO_FORTRAN_BOX(tile_box),
                                 BL_TO_FORTRAN_3D(f_u_tile),
                                 BL_TO_FORTRAN_3D(f_v_tile),
                                 BL_TO_FORTRAN_3D(f_w_tile),
                                 BL_TO_FORTRAN_3D(u_s_tile),
                                 BL_TO_FORTRAN_3D(v_s_tile),
                                 BL_TO_FORTRAN_3D(w_s_tile),
                                 BL_TO_FORTRAN_3D(u_d_tile),
                                 BL_TO_FORTRAN_3D(v_d_tile),
                                 BL_TO_FORTRAN_3D(w_d_tile),
                                 BL_TO_FORTRAN_3D(et_tile),
                                 & dt);
    }

    if (m_verbose > 1) {
        Real stoptime = ParallelDescriptor::second() - strttime;

        ParallelDescriptor::ReduceRealMax(stoptime,ParallelDescriptor::IOProcessorNumber());

        amrex::Print() << "MFIXParticleContainer::PICMultiDeposition time: " << stoptime << '\n';
    }
}



void IBCore::IBForceDeposition ( MultiFab & force, MultiFab & vel_d, const MultiFab & vel_s, int lev, Real dt) {

    force.setVal(0.);
    vel_d.setVal(0.);


    const Real strttime = ParallelDescriptor::second();


    /****************************************************************************
     * Create local copies => which have the same BA/DM as f and u              *
     ***************************************************************************/

    std::unique_ptr<iMultiFab> et = std::unique_ptr<iMultiFab>(
	new iMultiFab(vel_s.boxArray(), vel_s.DistributionMap(), vel_s.nComp(), vel_s.nGrow())
	);

    et->setVal(0);

    //Ensure that the target BA is cell-centered:
    BoxArray ba_cc_target = convert(vel_s.boxArray(), IntVect{0, 0, 0});
    RemakeLevel (lev, 0., ba_cc_target, vel_s.DistributionMap());


#ifdef _OPENMP
#pragma omp parallel
#endif
    for(MFIter mfi(vel_d, true); mfi.isValid(); ++ mfi) {
        const Box & tile_box = mfi.tilebox();

              FArrayBox & vel_d_tile = vel_d[mfi];
        const FArrayBox & vel_s_tile = vel_s[mfi];

              auto & et_tile = (* et)[mfi];

        const FArrayBox & ls_tile    = (* ls)[mfi];
        const IArrayBox & iface_tile = (* tag_interface)[mfi];
        const FArrayBox & vel_tile   = (* ls_vel)[mfi];

        interpolate_ib_cc (BL_TO_FORTRAN_BOX(tile_box),
                           BL_TO_FORTRAN_3D(vel_d_tile),
                           BL_TO_FORTRAN_3D(vel_s_tile),
                           BL_TO_FORTRAN_3D(et_tile),
                           BL_TO_FORTRAN_3D(ls_tile),
                           BL_TO_FORTRAN_3D(iface_tile),
                           BL_TO_FORTRAN_3D(vel_tile));
    }


#ifdef _OPENMP
#pragma omp parallel
#endif
    for(MFIter mfi(force, true); mfi.isValid(); ++ mfi) {
        const Box & tile_box = mfi.tilebox();

              FArrayBox & force_tile = force[mfi];

        const auto & et_tile = (* et)[mfi];

        const FArrayBox & vel_s_tile = vel_s[mfi];
        const FArrayBox & vel_d_tile = vel_d[mfi];


        fill_force_ib_cc (BL_TO_FORTRAN_BOX(tile_box),
                          BL_TO_FORTRAN_3D(force_tile),
                          BL_TO_FORTRAN_3D(vel_s_tile),
                          BL_TO_FORTRAN_3D(vel_d_tile),
                          BL_TO_FORTRAN_3D(et_tile),
                          & dt);

    }

    if (m_verbose > 1) {
        Real stoptime = ParallelDescriptor::second() - strttime;

        ParallelDescriptor::ReduceRealMax(stoptime,ParallelDescriptor::IOProcessorNumber());

        amrex::Print() << "MFIXParticleContainer::PICMultiDeposition time: " << stoptime << '\n';
    }
}



void IBCore::ImplicitDeposition (      MultiFab & f_u,       MultiFab & f_v,       MultiFab & f_w,
                                       MultiFab & u_d,       MultiFab & v_d,       MultiFab & w_d,
                                 const MultiFab & u_s, const MultiFab & v_s, const MultiFab & w_s,
                                 int lev, Real dt) {

    f_u.setVal(0.); f_v.setVal(0.); f_w.setVal(0.);
    u_d.setVal(0.); v_d.setVal(0.); w_d.setVal(0.);


    const Real strttime = ParallelDescriptor::second();


    /****************************************************************************
     * Create local copies => which have the same BA/DM as f and u              *
     ***************************************************************************/

    // Ensure that the target BA is cell-centered:
    BoxArray ba_cc_target = convert(u_s.boxArray(), IntVect::TheCellVector());


    std::unique_ptr<iMultiFab> et = std::unique_ptr<iMultiFab>(
            new iMultiFab(ba_cc_target, u_s.DistributionMap(), 1, u_s.nGrow())
        );

    et->setVal(0);


    // Perhaps this should not be necessary? If we use IBCore properly, then
    // regridding like this is not necessary.
    RemakeLevel (lev, 0., ba_cc_target, u_s.DistributionMap());


#ifdef _OPENMP
#pragma omp parallel
#endif
    for(MFIter mfi(* et, true); mfi.isValid(); ++ mfi) {
        // Iterate over cell-centered MultiFab (et) as reference for
        // face-centered data

        const Box & tile_box = mfi.tilebox();

        FArrayBox & u_d_tile = u_d[mfi];
        FArrayBox & v_d_tile = v_d[mfi];
        FArrayBox & w_d_tile = w_d[mfi];

        auto & et_tile = (* et)[mfi];

        const FArrayBox & u_s_tile = u_s[mfi];
        const FArrayBox & v_s_tile = v_s[mfi];
        const FArrayBox & w_s_tile = w_s[mfi];

        const FArrayBox & ls_tile    = (* ls)[mfi];
        const IArrayBox & iface_tile = (* tag_interface)[mfi];
        const FArrayBox & vel_tile   = (* ls_vel)[mfi];

        interpolate_ib_staggered (BL_TO_FORTRAN_BOX(tile_box),
                                  BL_TO_FORTRAN_3D(u_d_tile),
                                  BL_TO_FORTRAN_3D(v_d_tile),
                                  BL_TO_FORTRAN_3D(w_d_tile),
                                  BL_TO_FORTRAN_3D(u_s_tile),
                                  BL_TO_FORTRAN_3D(v_s_tile),
                                  BL_TO_FORTRAN_3D(w_s_tile),
                                  BL_TO_FORTRAN_3D(et_tile),
                                  BL_TO_FORTRAN_3D(ls_tile),
                                  BL_TO_FORTRAN_3D(iface_tile),
                                  BL_TO_FORTRAN_3D(vel_tile));
    }
    
    u_d.FillBoundary(Geom(lev).periodicity());
    v_d.FillBoundary(Geom(lev).periodicity());
    w_d.FillBoundary(Geom(lev).periodicity());

    et->FillBoundary(Geom(lev).periodicity());


#ifdef _OPENMP
#pragma omp parallel
#endif
    for(MFIter mfi(* et, true); mfi.isValid(); ++ mfi) {
        // Iterate over cell-centered MultiFab (et) as reference for
        // face-centered data

        const Box & tile_box = mfi.tilebox();

        FArrayBox & f_u_tile = f_u[mfi];
        FArrayBox & f_v_tile = f_v[mfi];
        FArrayBox & f_w_tile = f_w[mfi];

        const auto & et_tile = (* et)[mfi];

        fill_fgds_ib (BL_TO_FORTRAN_BOX(tile_box),
        		      BL_TO_FORTRAN_3D(f_u_tile),
		              BL_TO_FORTRAN_3D(f_v_tile),
		              BL_TO_FORTRAN_3D(f_w_tile),
		              BL_TO_FORTRAN_3D(et_tile));

    }

    f_u.FillBoundary(Geom(lev).periodicity());
    f_v.FillBoundary(Geom(lev).periodicity());
    f_w.FillBoundary(Geom(lev).periodicity());


    if (m_verbose > 1) {
        Real stoptime = ParallelDescriptor::second() - strttime;

        ParallelDescriptor::ReduceRealMax(stoptime,ParallelDescriptor::IOProcessorNumber());

        amrex::Print() << "MFIXParticleContainer::PICMultiDeposition time: " << stoptime << '\n';
    }
}


void IBCore::InterpolateForce ( const std::array<MultiFab, AMREX_SPACEDIM> & force,
                                int lev, std::pair<int,int> part_index,
                                std::array<Real, AMREX_SPACEDIM> & f_trans) const {


    std::array<MultiFab, AMREX_SPACEDIM> force_ibm;
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        force_ibm[d].define(
                convert(grids[lev], nodal_flag_dir[d]), dmap[lev],
                1, ib_pc->get_nghost()
                ); // TODO: No need for so many ghost cells after Andrew's done
        force_ibm[d].setVal(0.);
        f_trans[d] = 0.;
    }


    bool has_part = false;
    int index_ibm = 0;

    auto part_it = std::find(part_index_loc.begin(), part_index_loc.end(), part_index);
    if ( part_it < part_index_loc.end() ) {
        has_part = true;
        index_ibm = std::distance(part_index_loc.begin(), part_it);

        std::cout << index_ibm << std::endl;
    }


    if ( has_part ) {

        MultiFab dummy(grids[lev], dmap[lev], 1, 1);

#ifdef _OPENMP
#pragma omp parallel
#endif
        for(MFIter mfi(dummy , true); mfi.isValid(); ++ mfi) {
            // Iterate over cell-centered MultiFab `dummy` as reference for
            // face-centered data

            const Box & tile_box = mfi.growntilebox();

            FArrayBox & f_u_tile = force_ibm[0][mfi];
#if (AMREX_SPACEDIM > 1)
            FArrayBox & f_v_tile = force_ibm[1][mfi];
#endif
#if (AMREX_SPACEDIM > 2)
            FArrayBox & f_w_tile = force_ibm[2][mfi];
#endif

            const auto & tag_tile = iface_tags_loc[index_ibm][mfi];

            fill_fgds_ib (BL_TO_FORTRAN_BOX(tile_box),
                          BL_TO_FORTRAN_ANYD(f_u_tile),
#if (AMREX_SPACEDIM > 1)
                          BL_TO_FORTRAN_ANYD(f_v_tile),
#endif
#if (AMREX_SPACEDIM > 2)
                          BL_TO_FORTRAN_ANYD(f_w_tile),
#endif
                          BL_TO_FORTRAN_ANYD(tag_tile));

        }

    }

    ParallelDescriptor::Barrier();

    // Using growntilebox => don't fill boundary:
    // for (int d=0; d<AMREX_SPACEDIM; ++d)
    //     force_ibm[d].FillBoundary(Geom(lev).periodicity());

    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        MultiFab::Multiply(force_ibm[d], force[d], 0, 0, 1, 0);

        // Not multiplying in the ghost region => do fill boundary
        force_ibm[d].FillBoundary(Geom(lev).periodicity());


        VisMF::Write(force_ibm[d], "force_ibm_" + std::to_string(d));

        f_trans[d] = force_ibm[d].sum();
    }
}



void IBCore::save_levelset_data() {
    // VisMF::Write only does MultiFabs, so we need to convert the ls_id
    // iMultiFab -> MultiFab

    std::unique_ptr<MultiFab> dbl_ls_id = IBMFUtil::convert<iMultiFab, MultiFab>(* ls_id);
    std::unique_ptr<MultiFab> dbl_tag_iface = IBMFUtil::convert<iMultiFab, MultiFab>(* tag_interface);


    MultiFab nodal_ls(amrex::convert(ls->boxArray(), IntVect{0, 0, 0}),
                      ls->DistributionMap(), ls->nComp(), ls->nGrow());
    amrex::average_node_to_cellcenter(nodal_ls, 0, * ls, 0, 1);

    MultiFab nodal_vel(amrex::convert(ls_vel->boxArray(), IntVect{0,0,0}),
                       ls_vel->DistributionMap(), ls_vel->nComp(), ls_vel->nGrow());
    amrex::average_node_to_cellcenter(nodal_vel, 0, * ls_vel, 0, 3);


    VisMF::Write(  nodal_ls,      "ib_data/ls");
    VisMF::Write(  nodal_vel,     "ib_data/ls_vel");
    VisMF::Write(* dbl_ls_id,     "ib_data/ls_id");
    VisMF::Write(* dbl_tag_iface, "ib_data/iface");

}
