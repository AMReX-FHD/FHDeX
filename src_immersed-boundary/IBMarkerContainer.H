#ifndef _IBMarkerContainer_H_
#define _IBMarkerContainer_H_


#include <AMReX.H>
#include <AMReX_Vector.H>
#include <AMReX_RealVect.H>
#include <AMReX_AmrParGDB.H>
#include <AMReX_Particles.H>
#include <AMReX_Periodicity.H>
#include <AMReX_NeighborParticles.H>
#include <IBParticleInfo.H>
#include <common_namespace.H>
#include <IBMarkerContainerBase.H>


using namespace amrex;


// IBM => Immmersed Boundary Marker
struct IBM_realData {
    //Analogous to particle realData (p.m_data)
    enum {
        radius = 0,
        velx,
        vely,
        velz,
        forcex,
        forcey,
        forcez,
        pred_posx,
        pred_posy,
        pred_posz,
        pred_velx,
        pred_vely,
        pred_velz,
        pred_forcex,
        pred_forcey,
        pred_forcez,
        count    // Awesome little trick! (only works if first field is 0)
    };

    static Vector<std::string> names() {
        return Vector<std::string> {
            "radius",
            "velx",
            "vely",
            "velz",
            "forcex",
            "forcey",
            "forcez",
            "pred_posx",
            "pred_posy",
            "pred_posz",
            "pred_velx",
            "pred_vely",
            "pred_velz",
            "pred_forcex",
            "pred_forcey",
            "pred_forcez"
        };
    };
};



struct IBM_intData {
    //Analogous to particle intData (p.m_data)
    enum {
        id_0 = 0,
        cpu_0,
        id_1,
        cpu_1,
        count
    };

    static Vector<std::string> names() {
        return Vector<std::string> {
            "id_0",
            "cpu_0",
            "id_1",
            "cpu_1",
        };
    };
};



// class IBMarIter
//     : public IBMarIterBase<IBM_realData::count, IBM_intData::count>
// {
//
// public:
//
//     using IBMarIterBase<IBM_realData::count, IBM_intData::count>::ParIter;
// };

class IBMarIter
    : public ParIter<IBM_realData::count, IBM_intData::count, 0, 0>
{

public:

    using ParIter<IBM_realData::count, IBM_intData::count, 0, 0>::ParIter;
};




class IBMarkerContainer
    : public IBMarkerContainerBase<IBM_realData, IBM_intData>
{

public:

    using IBMarkerContainerBase<IBM_realData, IBM_intData>
        ::IBMarkerContainerBase;

    using MyConstIBMarIter = IBMarkerContainerBase<IBM_realData, IBM_intData>
        ::MyConstIBMarIter;



    IBMarkerContainer(AmrCore * amr_core, int n_nbhd);
    IBMarkerContainer(const Geometry & geom, const DistributionMapping & dmap,
                      const BoxArray & ba, int n_nbhd);

    virtual ~IBMarkerContainer() {};


    // (ID, initial CPU) tuple: unique to each particle
    using MarkerIndex = typename IBMarkerContainerBase<IBM_realData, IBM_intData>
        ::PairIndex;


    void InitList(int lev, const Vector<Real> & radius, const Vector<RealVect> & pos);
    void InitSingle(int lev, Real radius, const RealVect & pos, int id, int cpu, int i);


    // Get number of particles
    int NumberOfMarkers(IBMarIter & pti){ return pti.GetArrayOfStructs().size();};



    /****************************************************************************
     *                                                                          *
     * Marker Spreading and Interpolation                                       *
     *                                                                          *
     ***************************************************************************/


    //___________________________________________________________________________
    // Spreading operators applied to external forces
    void SpreadMarkers(int lev,
                       const Vector<RealVect> & f_in,
                       std::array<MultiFab, AMREX_SPACEDIM> & f_out) const;

    void SpreadMarkers(int lev,
                       const Vector<RealVect> & f_in,
                       std::array<MultiFab, AMREX_SPACEDIM> & f_out,
                       std::array<MultiFab, AMREX_SPACEDIM> & f_weights) const;
    //---------------------------------------------------------------------------



    //___________________________________________________________________________
    // Spreading operators applied to internal force fields (in amrex::Particles)
    using IBMarkerContainerBase<IBM_realData, IBM_intData>::SpreadMarkers;

    using IBMarkerContainerBase<IBM_realData, IBM_intData>::SpreadPredictor;
    //---------------------------------------------------------------------------



    //___________________________________________________________________________
    // Interpolation operators applied to external forces
    void InterpolateMarkers(int lev,
                            Vector<RealVect> & f_out,
                            const std::array<MultiFab, AMREX_SPACEDIM> & f_in) const;

    void InterpolateMarkers(int lev,
                            Vector<RealVect> & f_out,
                            const std::array<MultiFab, AMREX_SPACEDIM> & f_in,
                            const std::array<MultiFab, AMREX_SPACEDIM> & f_weights) const;
    //---------------------------------------------------------------------------



    //___________________________________________________________________________
    // Interpolation operators applied to internal force fields (in amrex::Particles)
    using IBMarkerContainerBase<IBM_realData, IBM_intData>::InterpolateMarkers;

    using IBMarkerContainerBase<IBM_realData, IBM_intData>::InterpolatePredictor;
    //---------------------------------------------------------------------------



    /****************************************************************************
     *                                                                          *
     * Access marker (amrex particle) data                                      *
     *                                                                          *
     ***************************************************************************/


    const IntVector & get_neighbor_list(int lev, const PairIndex & index) const {
        return neighbor_list[lev].at(index);
    };

    ParticleType * get_neighbors(int lev, const PairIndex & index) const {
        return (ParticleType *) neighbors[lev].at(index).dataPtr();
    };

    int FindConnectedMarkers(      AoS & particles,
                             const ParticleType & part,
                                   ParticleVector & nbhd_data,
                             const Vector<int> & nbhd,
                             int nbhd_index,
                             ParticleType *& prev_marker,
                             ParticleType *& next_marker) const;


    void LocalIBMarkerInfo(Vector<IBM_info> & marker_info, int lev, PairIndex index,
                           bool unique = false) const;
    Vector<IBM_info> LocalIBMarkerInfo(int lev, PairIndex index,
                                       bool unique = false) const;
    Vector<IBM_info> LocalIBMarkerInfo(int lev, bool unique = false) const;


    void NeighborIBMarkerInfo(Vector<IBM_info> & npart_info, int lev, PairIndex index,
                              bool unique = false) const ;
    Vector<IBM_info> NeighborIBMarkerInfo(int lev, PairIndex index,
                                          bool unique = false) const;
    Vector<IBM_info> NeighborIBMarkerInfo(int lev, bool unique = false) const;


    void IBMarkerInfo(Vector<IBM_info> & part_info, int lev, PairIndex index,
                        bool unique = false) const;
    Vector<IBM_info> IBMarkerInfo(int lev, PairIndex index,
                                  bool unique = false) const;
    Vector<IBM_info> IBMarkerInfo(int lev, bool unique = false) const;




    static bool use_neighbor_list;
    static bool sort_neighbor_list;



    struct {
        bool operator()(const ParticleType & p1, const ParticleType & p2) const {
            return AMREX_D_TERM(   (p1.pos(0) - p2.pos(0))*(p1.pos(0) - p2.pos(0)) ,
                                 + (p1.pos(1) - p2.pos(1))*(p1.pos(1) - p2.pos(1)) ,
                                 + (p1.pos(2) - p2.pos(2))*(p1.pos(2) - p2.pos(2)) )
                <= 4.0*p1.rdata(0)*p1.rdata(0);
        }
    } CheckPair;
};


#endif
