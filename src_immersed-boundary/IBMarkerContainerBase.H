#ifndef _IBMarkerContainerBase_H_
#define _IBMarkerContainerBase_H_


#include <AMReX.H>
#include <AMReX_Vector.H>
#include <AMReX_RealVect.H>
#include <AMReX_AmrParGDB.H>
#include <AMReX_Particles.H>
#include <AMReX_Periodicity.H>
#include <AMReX_NeighborParticles.H>
#include <IBParticleInfo.H>
#include <common_namespace.H>


using namespace amrex;



template <int NStructReal, int NStructInt>
class IBMarIterBase
    : public ParIter<NStructReal, NStructInt, 0, 0>
{

public:

    using ContainerType = ParticleContainer<NStructReal, NStructInt, 0, 0>;

    IBMarIterBase (ContainerType & pc, int level)
        : ParIter<NStructReal, NStructInt, 0, 0>(pc,level)
        {}

    IBMarIterBase (ContainerType & pc, int level, MFItInfo& info)
        : ParIter<NStructReal, NStructInt, 0, 0>(pc,level,info)
        {}
};



template <typename StructReal, typename StructInt>
class IBMarkerContainerBase
    : public NeighborParticleContainer<StructReal::count, StructInt::count>
{

public:

    using NeighborParticleContainer<StructReal::count, StructInt::count>
          ::NeighborParticleContainer;

    using MyConstIBMarIter = ParConstIter<StructReal::count, StructInt::count, 0, 0>;
    using MyIBMarIter      = IBMarIterBase<StructReal::count, StructInt::count>;

    using ParticleType = typename NeighborParticleContainer<StructReal::count, StructInt::count>::ParticleType;
    using PairIndex = typename NeighborParticleContainer<StructReal::count, StructInt::count>::PairIndex;
    using AoS = typename NeighborParticleContainer<StructReal::count, StructInt::count>::AoS;


    IBMarkerContainerBase(AmrCore * amr_core, int n_nbhd);
    IBMarkerContainerBase(const Geometry & geom, const DistributionMapping & dmap,
                          const BoxArray & ba, int n_nbhd);

    virtual ~IBMarkerContainerBase() {};


    // (ID, initial CPU) tuple: unique to each particle
    using MarkerIndex = std::pair<int, int>;


    // Get number of particles
    int NumberOfMarkers(IBMarIterBase<StructReal::count, StructInt::count> & pti){
        return pti.GetArrayOfStructs().size();
    };



    /****************************************************************************
     *                                                                          *
     * Marker Spreading and Interpolation                                       *
     *                                                                          *
     ***************************************************************************/


    //___________________________________________________________________________
    // Static spreading operators applied to raw data arrays
    static void SpreadMarkers(
                const Vector<RealVect> & f_in,
                const Vector<RealVect> & f_pos,
                      std::array<MultiFab, AMREX_SPACEDIM> & f_out,
                      std::array<MultiFab, AMREX_SPACEDIM> & f_weights,
                const std::array<MultiFab, AMREX_SPACEDIM> & coords,
                const Real* dx, int ghost
            );

    static void SpreadMarkers(
                const Vector<RealVect> & f_in,
                const Vector<RealVect> & f_pos,
                const Box & tile_box,
                      std::array<      FArrayBox *, AMREX_SPACEDIM> & f_out,
                      std::array<      FArrayBox *, AMREX_SPACEDIM> & f_weights,
                const std::array<const FArrayBox *, AMREX_SPACEDIM> & coords,
                const Real* dx, int ghost
            );
    //---------------------------------------------------------------------------



    //___________________________________________________________________________
    // Spreading operators applied to internal force fields (in amrex::Particles)
    void SpreadMarkers(int lev,
                       std::array<MultiFab, AMREX_SPACEDIM> & f_out) const;

    void SpreadMarkers(int lev,
                             std::array<MultiFab, AMREX_SPACEDIM> & f_out,
                       const std::array<MultiFab, AMREX_SPACEDIM> & coords,
                       const Real* dx, int ghost) const;

    void SpreadPredictor(int lev,
                         std::array<MultiFab, AMREX_SPACEDIM> & f_out) const;
    //---------------------------------------------------------------------------



    //___________________________________________________________________________
    // Static interpolation operators applied to raw data arrays
    static void InterpolateMarkers(
                      Vector<RealVect> & f_out,
                const Vector<RealVect> & f_pos,
                const std::array<MultiFab, AMREX_SPACEDIM> & f_in,
                const std::array<MultiFab, AMREX_SPACEDIM> & f_weights,
                const std::array<MultiFab, AMREX_SPACEDIM> & coords,
                const Real * dx, const int ghost
            );

    static void InterpolateMarkers(
                      Vector<RealVect> & f_out,
                const Vector<RealVect> & f_pos,
                const Box & tile_box,
                const std::array<const FArrayBox *, AMREX_SPACEDIM> & f_in,
                const std::array<const FArrayBox *, AMREX_SPACEDIM> & f_weights,
                const std::array<const FArrayBox *, AMREX_SPACEDIM> & coords,
                const Real * dx, const int ghost
            );
    //---------------------------------------------------------------------------



    //___________________________________________________________________________
    // Interpolation operators applied to internal force fields (in amrex::Particles)
    void InterpolateMarkers(int lev,
                            const std::array<MultiFab, AMREX_SPACEDIM> & f_in);

    void InterpolatePredictor(int lev,
                              const std::array<MultiFab, AMREX_SPACEDIM> & f_in);
    //---------------------------------------------------------------------------



    /****************************************************************************
     *                                                                          *
     * Physics                                                                  *
     *                                                                          *
     ***************************************************************************/


    void MoveMarkers(int lev, Real dt);
    void MovePredictor(int lev, Real dt);
    void ResetMarkers(int lev);
    void ResetPredictor(int lev);



    /****************************************************************************
     *                                                                          *
     * Stochastic drift                                                         *
     *                                                                          *
     ***************************************************************************/

    void RFD(int lev, const Real* dx,
                   std::array<MultiFab, AMREX_SPACEDIM> & f_out,
             const std::array<MultiFab, AMREX_SPACEDIM> & coords);



    /****************************************************************************
     *                                                                          *
     * Access marker (amrex particle) data                                      *
     *                                                                          *
     ***************************************************************************/

    void PullDown(int lev, Real * list, int element, int totalParticles);
    void PushUp(int lev, Real * list, int element, int totalParticles);

    int get_nghost() const {return nghost;};

    void PrintMarkerData(int lev) const;

    const Vector<std::array<MultiFab, AMREX_SPACEDIM>> & get_face_coords() const {
        return  face_coords;
    }
    const std::array<MultiFab, AMREX_SPACEDIM> & get_face_coords(int lev) const {
        return  face_coords[lev];
    }


protected:

    int nghost;

    void InitInternals(int ngrow);
    void ReadStaticParameters();

    AmrCore * m_amr_core;

private:

    // Positions on faces
    Vector<std::array<MultiFab, AMREX_SPACEDIM>> face_coords;
};



#include "IBMarkerContainerBaseI.H"


#endif
