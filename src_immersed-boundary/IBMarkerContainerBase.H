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
    : public amrex::ParIter<NStructReal, NStructInt, 0, 0>
{

public:

    using ParIter<NStructReal, NStructInt, 0, 0>::ParIter;
};



template <typename StructReal, typename StructInt>
class IBMarkerContainerBase
    : public amrex::NeighborParticleContainer<StructReal::count, StructInt::count>
{

public:

    using NeighborParticleContainer<StructReal::count, StructInt::count>
          ::NeighborParticleContainer;

    using MyConstParIter = ParConstIter<StructReal::count, StructInt::count, 0, 0>;

    using ParticleType = typename NeighborParticleContainer<StructReal::count, StructInt::count>::ParticleType;
    using PairIndex = typename NeighborParticleContainer<StructReal::count, StructInt::count>::PairIndex;
    using AoS = typename NeighborParticleContainer<StructReal::count, StructInt::count>::AoS;


    IBMarkerContainerBase(AmrCore * amr_core, int n_nbhd);
    IBMarkerContainerBase(const Geometry & geom, const DistributionMapping & dmap,
                          const BoxArray & ba, int n_nbhd);

    virtual ~IBMarkerContainerBase() {};


    // (ID, initial CPU) tuple: unique to each particle
    using MarkerIndex = std::pair<int, int>;


    // void InitList(int lev, const Vector<Real> & radius, const Vector<RealVect> & pos);
    // void InitSingle(int lev, Real radius, const RealVect & pos, int id, int cpu, int i);


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
    // Spreading operators applied to raw data arrays
    // void SpreadMarkers(int lev,
    //                    const Vector<RealVect> & f_in,
    //                    std::array<MultiFab, AMREX_SPACEDIM> & f_out) const;

    // void SpreadMarkers(int lev,
    //                    const Vector<RealVect> & f_in,
    //                    std::array<MultiFab, AMREX_SPACEDIM> & f_out,
    //                    std::array<MultiFab, AMREX_SPACEDIM> & f_weights) const;

    void SpreadMarkers(int lev,
                       const Vector<RealVect> & f_in,
                       const Vector<RealVect> & f_pos,
                       std::array<MultiFab, AMREX_SPACEDIM> & f_out,
                       std::array<MultiFab, AMREX_SPACEDIM> & f_weights) const;
    //---------------------------------------------------------------------------



    //___________________________________________________________________________
    // Spreading operators applied to internal force fields (in amrex::Particles)
    void SpreadMarkers(int lev,
                       std::array<MultiFab, AMREX_SPACEDIM> & f_out) const;

    void SpreadPredictor(int lev,
                         std::array<MultiFab, AMREX_SPACEDIM> & f_out) const;
    //---------------------------------------------------------------------------



    //___________________________________________________________________________
    // Interpolation operators applied to raw data arrays
    // void InterpolateMarkers(int lev,
    //                         Vector<RealVect> & f_out,
    //                         const std::array<MultiFab, AMREX_SPACEDIM> & f_in) const;

    // void InterpolateMarkers(int lev,
    //                         Vector<RealVect> & f_out,
    //                         const std::array<MultiFab, AMREX_SPACEDIM> & f_in,
    //                         const std::array<MultiFab, AMREX_SPACEDIM> & f_weights) const;

    void InterpolateMarkers(int lev,
                            Vector<RealVect> & f_out,
                            const Vector<RealVect> & f_pos,
                            const std::array<MultiFab, AMREX_SPACEDIM> & f_in,
                            const std::array<MultiFab, AMREX_SPACEDIM> & f_weights) const;
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
     * Access marker (amrex particle) data                                      *
     *                                                                          *
     ***************************************************************************/


    // const IntVector & get_neighbor_list(int lev, const PairIndex & index) const {
    //     return neighbor_list[lev].at(index);
    // };

    // ParticleType * get_neighbors(int lev, const PairIndex & index) const {
    //     return (ParticleType *) neighbors[lev].at(index).dataPtr();
    // };

    // int FindConnectedMarkers(      AoS & particles,
    //                          const ParticleType & part,
    //                                ParticleVector & nbhd_data,
    //                          const Vector<int> & nbhd,
    //                          int nbhd_index,
    //                          ParticleType *& prev_marker,
    //                          ParticleType *& next_marker) const;



    int get_nghost() const {return nghost;};

    void PrintMarkerData(int lev) const;

    // void LocalIBMarkerInfo(Vector<IBM_info> & marker_info, int lev, PairIndex index,
    //                        bool unique = false) const;
    // Vector<IBM_info> LocalIBMarkerInfo(int lev, PairIndex index,
    //                                    bool unique = false) const;
    // Vector<IBM_info> LocalIBMarkerInfo(int lev, bool unique = false) const;


    // void NeighborIBMarkerInfo(Vector<IBM_info> & npart_info, int lev, PairIndex index,
    //                           bool unique = false) const ;
    // Vector<IBM_info> NeighborIBMarkerInfo(int lev, PairIndex index,
    //                                       bool unique = false) const;
    // Vector<IBM_info> NeighborIBMarkerInfo(int lev, bool unique = false) const;


    // void IBMarkerInfo(Vector<IBM_info> & part_info, int lev, PairIndex index,
    //                   bool unique = false) const;
    // Vector<IBM_info> IBMarkerInfo(int lev, PairIndex index,
    //                               bool unique = false) const;
    // Vector<IBM_info> IBMarkerInfo(int lev, bool unique = false) const;




    static bool use_neighbor_list;
    static bool sort_neighbor_list;



    // Needs to be specified once a concrete class is created
    // struct {
    //     bool operator()(const ParticleType & p1, const ParticleType & p2) const {
    //         return AMREX_D_TERM(   (p1.pos(0) - p2.pos(0))*(p1.pos(0) - p2.pos(0)) ,
    //                              + (p1.pos(1) - p2.pos(1))*(p1.pos(1) - p2.pos(1)) ,
    //                              + (p1.pos(2) - p2.pos(2))*(p1.pos(2) - p2.pos(2)) )
    //             <= 4.0*p1.rdata(0)*p1.rdata(0);
    //     }
    // } CheckPair;

    const Vector<std::array<MultiFab, AMREX_SPACEDIM>> & get_face_coords() const {
        return  face_coords;
    }


private:

    // Positions on faces
    Vector<std::array<MultiFab, AMREX_SPACEDIM>> face_coords;

    void InitInternals(int ngrow);
    void ReadStaticParameters();

    AmrCore * m_amr_core;
    int nghost;
};



#include "IBMarkerContainerBaseI.H"


#endif
