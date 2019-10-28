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
    using IBMarIterBase<IBBReal::count, IBBInt::count>::IBMarIterBase;
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

    // Add marker to database **without** redistributing (this is safe to call
    // within a `ParIter`/`MFIter loop), but the user will have to call
    // `BlobContainer::Redistribute` outside the loop.
    void AddSingle(int lev, Real radius, const RealVect & pos, int id, int cpu, int i);
    // Add marker to database **and** redistribute this particle (this is
    // **unsafe** to call **within** a `ParIter`/`MFIter loop). The user will
    // not have to call `BlobContainer::Redistribute` outside the loop.
    void InitSingle(int lev, Real radius, const RealVect & pos, int id, int cpu, int i);

    // Calls parent class MovePredictor, initilizing pred_pos[x,y,z] from
    // ref_del[x,y,z] fields (rather than 0);
    void MovePredictor(int lev, Real dt);
    // Calls parent class MoveMarkers, but also updates the ref_del[x,y,z] fields
    void MoveMarkers(int lev, Real dt);

    // Computes flexible for model
    void PredictorForces(int lev, Real k);
    void MarkerForces(int lev, Real k);
};



//TODO: Check members lables with "Audit"


// IBP => Immmersed Boundary Particle
struct IBMBReal {
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



struct IBMBInt {
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
    : public amrex::ParIter<IBMBReal::count, IBMBInt::count, 0, 0>
{

public:

    using ParIter<IBMBReal::count, IBMBInt::count, 0, 0>::ParIter;
};




class IBMultiBlobContainer
    : public amrex::NeighborParticleContainer<IBMBReal::count, IBMBInt::count>
{

public:

    IBMultiBlobContainer(AmrCore * amr_core, int n_nbhd);
    IBMultiBlobContainer(const Geometry & geom, const DistributionMapping & dmap,
                         const BoxArray & ba, int n_nbhd);

    // TODO: Audit
    virtual ~IBMultiBlobContainer() {};


    // (ID, initial CPU) tuple: unique to each particle
    using MarkerIndex = typename IBMarkerContainerBase<IBMReal, IBMInt>::PairIndex;
    // indexing particles and corresponding spot in the neighbor list
    using MarkerListIndex = std::pair<int, int>;
    // indexing tiles (box index, local tile index)
    using TileIndex = std::pair<int, int>;

    // Get number of particles
    int NumberOfParticles(IBMBIter & pti){
        return pti.GetArrayOfStructs().size();
    };



    /****************************************************************************
     *                                                                          *
     * Initialization                                                           *
     *                                                                          *
     ***************************************************************************/


    void InitSingle(int lev, const RealVect & pos, Real r, Real rho);
    void FillMarkerPositions(int lev, int n_marker);



    /****************************************************************************
     *                                                                          *
     * Resetting functions for internal markers                                 *
     *                                                                          *
     ***************************************************************************/


    void ResetPredictor(int lev);
    void ResetMarkers(int lev);


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



    /****************************************************************************
     *                                                                          *
     * Physics                                                                  *
     *                                                                          *
     ***************************************************************************/


    void MoveMarkers(int lev, Real dt);
    void MovePredictor(int lev, Real dt);

    void PredictorForces(int lev, Real k);
    void MarkerForces(int lev, Real k);



    /****************************************************************************
     *                                                                          *
     * IO                                                                       *
     *                                                                          *
     ***************************************************************************/

    void WritePlotFile(const std::string & dir, const std::string & name,
                       const std::string & blob_name) const;


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
