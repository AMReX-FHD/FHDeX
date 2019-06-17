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


using namespace amrex;


// IBM => Immmersed Boundary Marker
struct IBM_realData {
    //Analogous to particle realData (p.m_data)
    enum {
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
        count    // Awesome little trick!
    };

    static Vector<std::string> names() {
        return Vector<std::string> {
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



class IBMarIter
    : public amrex::ParIter<IBM_realData::count, IBM_intData::count, 0, 0>
{

public:

    using ParIter<IBM_realData::count, IBM_intData::count, 0, 0>::ParIter;
};



class IBMarkerContainer
    : public amrex::NeighborParticleContainer<IBM_realData::count, IBM_intData::count>
{

public:

    using NeighborParticleContainer<IBM_realData::count, IBM_intData::count>
        ::NeighborParticleContainer;

    IBMarkerContainer(AmrCore * amr_core, int n_nbhd);
    IBMarkerContainer(const Geometry & geom, const DistributionMapping & dmap,
                      const BoxArray & ba, int n_nbhd);

    virtual ~IBMarkerContainer() {};


    // (ID, initial CPU) tuple: unique to each particle
    using MarkerIndex = std::pair<int, int>;


    void InitList(int lev, const Vector<RealVect> & pos);





    // Get number of particles
    int NumberOfMarkers(IBMarIter & pti){ return pti.GetArrayOfStructs().size();};


    void SpreadMarkers(int lev, const Vector<RealVect> & f_in,
                       std::array<MultiFab, AMREX_SPACEDIM> & f_out) const;

    void SpreadMarkers(int lev, const Vector<RealVect> & f_in,
                       std::array<MultiFab, AMREX_SPACEDIM> & f_out,
                       std::array<MultiFab, AMREX_SPACEDIM> & f_weights) const;


    void SpreadMarkers(int lev, std::array<MultiFab, AMREX_SPACEDIM> & f_out) const;

    void SpreadPredictor(int lev, std::array<MultiFab, AMREX_SPACEDIM> & f_out) const;



    void InterpolateMarkers(int lev, Vector<RealVect> & f_out,
                            const std::array<MultiFab, AMREX_SPACEDIM> & f_in) const;

    void InterpolateMarkers(int lev, Vector<RealVect> & f_out,
                            const std::array<MultiFab, AMREX_SPACEDIM> & f_in,
                            const std::array<MultiFab, AMREX_SPACEDIM> & f_weights) const;


    void InterpolateMarkers(int lev, const std::array<MultiFab, AMREX_SPACEDIM> & f_in) const;

    void InterpolatePredictor(int lev, const std::array<MultiFab, AMREX_SPACEDIM> & f_in) const;




    /****************************************************************************
     * Physics                                                                   *
     ****************************************************************************/

    void MoveMarkers(int lev, Real dt);
    void MovePredictor(int lev, Real dt);



    int get_nghost() const {return nghost;};

    void PrintMarkerData(int lev) const;

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




    void InitInternals(int ngrow);

    void ReadStaticParameters();


    static bool use_neighbor_list;
    static bool sort_neighbor_list;


    // Posotions on faces
    Vector<std::array<MultiFab, AMREX_SPACEDIM>> face_coords;



    AmrCore * m_amr_core;
    int nghost;
};


#endif
