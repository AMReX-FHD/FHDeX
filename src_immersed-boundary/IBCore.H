#ifndef _IBCore_H_
#define _IBCore_H_

#include <AMReX.H>
#include <AMReX_AmrCore.H>

#include <MFUtil.H>
#include <IBParticleInfo.H>
#include <IBParticleContainer.H>


using namespace amrex;



class IBCore
    : public AmrCore {

public:
    IBCore();

    virtual ~IBCore ();

    void InitData ();


    // TODO: Implement refinement. For now this does nothing.

    //! Tag cells for refinement. TagBoxArray tags is built on level lev grids.
    virtual void ErrorEst (int lev, TagBoxArray & tags, Real time,
                           int ngrow
        ) override;


    //! Make a new level from scratch using provided BoxArray and
    //! DistributionMapping. Only used during initialization.
    virtual void MakeNewLevelFromScratch (int lev, Real time,
                                          const BoxArray & ba, const DistributionMapping & dm
        ) override;

    // TODO: Implement refinement. For now this does nothing.

    //! Make a new level using provided BoxArray and DistributionMapping and
    //! fill with interpolated coarse level data.
    virtual void MakeNewLevelFromCoarse (int lev, Real time,
                                         const BoxArray & ba, const DistributionMapping & dm
        ) override;


    //! Remake an existing level using provided BoxArray and DistributionMapping
    //! and fill with existing fine and coarse data.
    virtual void RemakeLevel (int lev, Real time,
                              const BoxArray & ba, const DistributionMapping & dm
        ) override;

    // TODO: Implement refinement. For now this does nothing.

    //! Delete level data
    virtual void ClearLevel (int lev)
        override;

    ///
    /// Do the actual PIC deposition for CalcDragOnFluid
    ///

    // Explicit deposition for staggered grids
    void IBForceDeposition (      MultiFab & f_u,       MultiFab & f_v,       MultiFab & f_w,
                                  MultiFab & u_d,       MultiFab & v_d,       MultiFab & w_d,
                            const MultiFab & u_s, const MultiFab & v_s, const MultiFab & w_s,
                            int lev, Real dt);

    // Explicit deposition of cell-centered grids
    void IBForceDeposition ( MultiFab & force, MultiFab & vel_d, const MultiFab & vel_s, int lev, Real dt);

    // Implicit deposition for staggered grids
    void ImplicitDeposition (      MultiFab & f_u,       MultiFab & f_v,       MultiFab & f_w,
                                   MultiFab & u_d,       MultiFab & v_d,       MultiFab & w_d,
                             const MultiFab & u_s, const MultiFab & v_s, const MultiFab & w_s,
                             int lev, Real dt);


    // Set local pointer to IBParticleContainer (used for updating the level-set)
    void set_IBParticleContainer (IBParticleContainer * a_ib_pc) {
        ib_pc = a_ib_pc;
    }


    //Needs to be the same indexing system as IBParticleContainer
    using PairIndex = std::pair<int, int>;


private:

    /****************************************************************************
     * Underlying BoxArray and DistributionMapping                              *
     ***************************************************************************/
    BoxArray ba_cc, ba_nd;
    DistributionMapping dm;


    /****************************************************************************
     * Utilities                                                                *
     ***************************************************************************/

    void save_levelset_data ();


    /****************************************************************************
     * General settings                                                         *
     ***************************************************************************/

    int m_verbose = 0;


    /****************************************************************************
     * Data used for LEVEL SET description of the immersed boundaries           *
     ***************************************************************************/

    // Flag that is set to `true` whenever the BoxArray and DistributionMapping
    // have changed and we need to copy (regrid) the level-set data.

    bool update_loc = true;

    Vector<PairIndex> part_index_loc;
    Vector< MultiFab> level_sets_loc;
    Vector<iMultiFab> iface_tags_loc;

    // Signed distance function to to neares particle surface
    std::unique_ptr<MultiFab>  ls;

    // Identifier of nearest particle (ID, INIT-CPU, INDEX)
    std::unique_ptr<iMultiFab> ls_id;

    // Tags all level-set points within 1 step of the IB interface
    std::unique_ptr<iMultiFab> tag_interface;

    // Stores the IB node velocity corresponding to the closest point on IB
    std::unique_ptr<MultiFab>  ls_vel;




    /****************************************************************************
     * IB Particle Information                                                  *
     ***************************************************************************/

    int  n_ibm_loc;
    IBParticleContainer * ib_pc;

};

#endif
