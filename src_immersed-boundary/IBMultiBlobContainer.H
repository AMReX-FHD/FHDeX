#ifndef _IBMultiBlobContainer_H_
#define _IBMultiBlobContainer_H_

#include <AMReX.H>
#include <AMReX_AmrCore.H>
#include <AMReX_NeighborParticles.H>

#include <IBMarkerContainer.H>


using namespace amrex;



struct IBBReal {
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
        ref_delx,
        ref_dely,
        ref_delz,
        count
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
            "pred_forcez",
            "ref_delx",
            "ref_dely",
            "ref_delz"
        };
    };
};


using IBBInt = IBMInt;


class BlobIter
    : public IBMarIterBase<IBBReal::count, IBBInt::count>
{

// public:
//     using IBMarIterBase<IBBReal::count, IBBInt::count>::IBMarIterBase;

};


class BlobContainer
    : public IBMarkerContainerBase<IBBReal, IBBInt>
{

public:

    using IBMarkerContainerBase<IBBReal, IBBInt>
        ::IBMarkerContainerBase;


    BlobContainer(AmrCore * amr_core, int n_nbhd);
    BlobContainer(const Geometry & geom, const DistributionMapping & dmap,
                  const BoxArray & ba, int n_nbhd);

    virtual ~BlobContainer() {};


    void InitSingle(int lev, Real radius, const RealVect & pos, int id, int cpu, int i);
};



//TODO: Check members lables with "Audit"


// IBP => Immmersed Boundary Particle
struct IBMB_realData {
    //Analogous to particle realData (p.m_data)
    // TODO: orientation?
    enum {
        radius = 0,
        velx,
        vely,
        velz,
        omegax,
        omegay,
        omegaz,
        dragx,    // ..................................... TODO: Audit
        dragy,    // ..................................... TODO: Audit
        dragz,    // ..................................... TODO: Audit
        count     // Awesome little trick!
    };

    static Vector<std::string> names() {
        return Vector<std::string> {
            "radius",
            "velx",
            "vely",
            "velz",
            "omegax",
            "omegay",
            "omegaz",
            "dragx",
            "dragy",
            "dragz"
        };
    };
};



struct IBMB_intData {
    //Analogous to particle intData (p.m_data)
    enum {
        phase = 0,  // ................................... TODO: Audit
        state,      // ................................... TODO: Audit
        count
    };

    static Vector<std::string> names() {
        return Vector<std::string> {
            "phase",
            "state"
        };
    };
};



class IBMBIter
    : public amrex::ParIter<IBMB_realData::count, IBMB_intData::count, 0, 0>
{

public:

    using ParIter<IBMB_realData::count, IBMB_intData::count, 0, 0>::ParIter;
};



class IBMultiBlobContainer
    : public amrex::NeighborParticleContainer<IBMB_realData::count, IBMB_intData::count>
{

public:

    IBMultiBlobContainer(AmrCore * amr_core, int n_nbhd);
    IBMultiBlobContainer(const Geometry & geom, const DistributionMapping & dmap,
                         const BoxArray & ba, int n_nbhd);

    // TODO: Audit
    virtual ~IBMultiBlobContainer() {};


    // (ID, initial CPU) tuple: unique to each particle
    using ParticleIndex = std::pair<int, int>;
    using IBMarkerMap   = std::map<ParticleIndex, Vector<RealVect>>;


    // Get number of particles
    int NumberOfParticles(IBMBIter & pti){
        return pti.GetArrayOfStructs().size();
    };


    void InitSingle(int lev, const RealVect & pos, Real r, Real rho);


    void FillMarkerPositions(int lev, int n_marker);


    /****************************************************************************
     *                                                                          *
     * Marker Spreading and Interpolation                                       *
     *                                                                          *
     ***************************************************************************/


    //___________________________________________________________________________
    // Spreading operators applied to marker force fields
    void SpreadMarkers(int lev,
                       std::array<MultiFab, AMREX_SPACEDIM> & f_out) const;

    void SpreadPredictor(int lev,
                         std::array<MultiFab, AMREX_SPACEDIM> & f_out) const;
    //---------------------------------------------------------------------------

    //___________________________________________________________________________
    // Interpolation operators applied to marker force fields
    void InterpolateMarkers(int lev,
                            const std::array<MultiFab, AMREX_SPACEDIM> & f_in);

    void InterpolatePredictor(int lev,
                              const std::array<MultiFab, AMREX_SPACEDIM> & f_in);
    //---------------------------------------------------------------------------



    static bool use_neighbor_list;
    static bool sort_neighbor_list;



    struct {
        bool operator() (const ParticleType & p1, const ParticleType & p2) const {
            return AMREX_D_TERM(   (p1.pos(0) - p2.pos(0))*(p1.pos(0) - p2.pos(0)) ,
                                 + (p1.pos(1) - p2.pos(1))*(p1.pos(1) - p2.pos(1)) ,
                                 + (p1.pos(2) - p2.pos(2))*(p1.pos(2) - p2.pos(2)) )
                <= 4.0*p1.rdata(0)*p1.rdata(0);

        }
    } CheckPair;


private:


    // IB marker data
    BlobContainer markers;

    void InitInternals(int ngrow);
    void ReadStaticParameters();


    // TODO: this might not be used anymore:
    AmrCore * m_amr_core;
    int nghost;


};


#endif
