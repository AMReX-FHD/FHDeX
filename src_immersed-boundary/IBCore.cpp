#include <AMReX.H>
#include <AMReX_NeighborParticles.H>


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
    grids[lev] = a_ba;
    ba_nd      = amrex::convert(grids[lev], IntVect{1,1,1});



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

    ls->define(grids[lev], dmap[lev], 1, n_pad);
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
    * Fill Level Set and tags                                                  *
    ****************************************************************************/

    RealVect dx = RealVect(
        AMREX_D_DECL( Geom(lev).CellSize(0), Geom(lev).CellSize(1), Geom(lev).CellSize(2) )
        );

    // Fill neighbors as IB particles can overlap from nearby grids
    ib_pc->fillNeighbors();

    // ParIter skips tiles without particles => Iterate over MultiFab instead
    // of ParticleIter for the following:

#ifdef _OPENMP
#pragma omp parallel
#endif
    for(MFIter mfi(* ls, ib_pc->tile_size); mfi.isValid(); ++mfi) {
        // NOTE: mfi's tile size must match the ParticleContainer tile size
        // MuliFabs are indexed using a pair: (BoxArray index, tile index):
        PairIndex index(mfi.index(), mfi.LocalTileIndex());

        Vector<IBP_info> info = ib_pc->get_IBParticle_info(lev, index);
        int np = info.size();
        test_interface(info.dataPtr(), & np);

        const Box & tile_box = mfi.tilebox();
        auto & phi_tile      = (* ls)[mfi];
        auto & tag_tile      = (* ls_id)[mfi];
        auto & vel_tile      = (* ls_vel)[mfi];

        fill_levelset_ib (tile_box.loVect(),  tile_box.hiVect(),
                          info.dataPtr(),     & np,
                          phi_tile.dataPtr(), phi_tile.loVect(), phi_tile.hiVect(),
                          tag_tile.dataPtr(), tag_tile.loVect(), tag_tile.hiVect(),
                          vel_tile.dataPtr(), vel_tile.loVect(), vel_tile.hiVect(),
                          dx.dataPtr()                                              );

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
    for(MFIter mfi(* tag_interface, true); mfi.isValid(); ++ mfi) {
        //const Box & tile_box = mfi.tilebox();

        const FArrayBox & phi_tile = (* ls)[mfi];
        const IArrayBox & tag_tile = (* ls_id)[mfi];
              IArrayBox & iface_tile = (* tag_interface)[mfi];

        tag_interface_ib (iface_tile.dataPtr(), iface_tile.loVect(), iface_tile.hiVect(),
                          phi_tile.dataPtr(),   phi_tile.loVect(),   phi_tile.hiVect(),
                          tag_tile.dataPtr(),   tag_tile.loVect(),   tag_tile.hiVect()     );
    }

    tag_interface->FillBoundary(Geom(lev).periodicity());

    save_levelset_data();
}


// TODO: Implement refinement. For now this does nothing.
void IBCore::MakeNewLevelFromCoarse (int lev, Real time,
                                     const BoxArray & ba, const DistributionMapping & dm ) {

}


void IBCore::RemakeLevel (int lev, Real time,
                          const BoxArray & a_ba, const DistributionMapping & a_dm) {

    Periodicity periodicity = Geom(lev).periodicity();

    //bool ba_changed = (ba_cc != a_ba);
    bool ba_changed = (grids[lev] != a_ba);
    //bool dm_changed = (dm    != a_dm);
    bool dm_changed = (dmap[lev]  != a_dm);

    if (ba_changed || dm_changed) {
        update_loc = true;

        //ba_cc = a_ba;
        grids[lev] = a_ba;
        //ba_nd = convert(a_ba, IntVect{1,1,1});
        ba_nd = convert(grids[lev], IntVect{1,1,1});


        //dm = a_dm;
        dmap[lev] = a_dm;
    }


    // Only overwrite level if BA/DM have changed
    if (update_loc) {

        std::unique_ptr<MultiFab> ls_new =
            //MFUtil::duplicate<MultiFab>(ba_nd, dm, * ls, periodicity);
            IBMFUtil::duplicate<MultiFab>(ba_nd, dmap[lev], * ls, periodicity);
        ls = std::move(ls_new);

        std::unique_ptr<iMultiFab> ls_id_new =
            //MFUtil::duplicate<iMultiFab>(ba_nd, dm, * ls_id, periodicity);
            IBMFUtil::duplicate<iMultiFab>(ba_nd, dmap[lev], * ls_id, periodicity);
        ls_id = std::move(ls_id_new);

        std::unique_ptr<iMultiFab> tag_interface_new =
            //MFUtil::duplicate<iMultiFab>(ba_cc, dm, * tag_interface, periodicity);
            IBMFUtil::duplicate<iMultiFab>(grids[lev], dmap[lev], * tag_interface, periodicity);
        tag_interface = std::move(tag_interface_new);

        std::unique_ptr<MultiFab> ls_vel_new =
            //MFUtil::duplicate<MultiFab>(ba_nd, dm, * ls_vel, periodicity);
            IBMFUtil::duplicate<MultiFab>(ba_nd, dmap[lev], * ls_vel, periodicity);
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

        interpolate_ib_staggered (tile_box.loVect(), tile_box.hiVect(),
                                  BL_TO_FORTRAN_3D(u_d_tile),
                                  BL_TO_FORTRAN_3D(v_d_tile),
                                  BL_TO_FORTRAN_3D(w_d_tile),
                                  BL_TO_FORTRAN_3D(u_s_tile),
                                  BL_TO_FORTRAN_3D(v_s_tile),
                                  BL_TO_FORTRAN_3D(w_s_tile),
                                  BL_TO_FORTRAN_3D(et_tile),
                                  BL_TO_FORTRAN_3D(ls_tile),
                                  BL_TO_FORTRAN_3D(iface_tile),
                                  BL_TO_FORTRAN_3D(vel_tile)            );
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


        fill_force_ib_staggered (tile_box.loVect(), tile_box.hiVect(),
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
                                 & dt                                   );

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

        interpolate_ib_cc (tile_box.loVect(), tile_box.hiVect(),
                           BL_TO_FORTRAN_3D(vel_d_tile),
                           BL_TO_FORTRAN_3D(vel_s_tile),
                           BL_TO_FORTRAN_3D(et_tile),
                           BL_TO_FORTRAN_3D(ls_tile),
                           BL_TO_FORTRAN_3D(iface_tile),
                           BL_TO_FORTRAN_3D(vel_tile)            );
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


        fill_force_ib_cc (tile_box.loVect(), tile_box.hiVect(),
                          BL_TO_FORTRAN_3D(force_tile),
                          BL_TO_FORTRAN_3D(vel_s_tile),
                          BL_TO_FORTRAN_3D(vel_d_tile),
                          BL_TO_FORTRAN_3D(et_tile),
                          & dt                                   );

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


    const Real      strttime    = ParallelDescriptor::second();


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

        interpolate_ib_staggered (tile_box.loVect(), tile_box.hiVect(),
                                  BL_TO_FORTRAN_3D(u_d_tile),
                                  BL_TO_FORTRAN_3D(v_d_tile),
                                  BL_TO_FORTRAN_3D(w_d_tile),
                                  BL_TO_FORTRAN_3D(u_s_tile),
                                  BL_TO_FORTRAN_3D(v_s_tile),
                                  BL_TO_FORTRAN_3D(w_s_tile),
                                  BL_TO_FORTRAN_3D(et_tile),
                                  BL_TO_FORTRAN_3D(ls_tile),
                                  BL_TO_FORTRAN_3D(iface_tile),
                                  BL_TO_FORTRAN_3D(vel_tile)            );
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

        fill_fgds_ib (tile_box.loVect(), tile_box.hiVect(),
		      BL_TO_FORTRAN_3D(f_u_tile),
		      BL_TO_FORTRAN_3D(f_v_tile),
		      BL_TO_FORTRAN_3D(f_w_tile),
		      BL_TO_FORTRAN_3D(et_tile)                );

    }

    if (m_verbose > 1) {
        Real stoptime = ParallelDescriptor::second() - strttime;

        ParallelDescriptor::ReduceRealMax(stoptime,ParallelDescriptor::IOProcessorNumber());

        amrex::Print() << "MFIXParticleContainer::PICMultiDeposition time: " << stoptime << '\n';
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
