#ifndef _IBMultiBlobContainer_H_
#define _IBMultiBlobContainer_H_

#include <AMReX.H>
#include <AMReX_AmrCore.H>
#include <AMReX_NeighborParticles.H>


using namespace amrex;



//TODO: Check members lables with "Audit"


// IBP => Immmersed Boundary Particle
struct IBMB_realData {
    //Analogous to particle realData (p.m_data)
    // TODO: orientation?
    enum {
        velx = 0,
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


    // Get number of particles
    int NumberOfParticles(IBMBIter & pti){
        return pti.GetArrayOfStructs().size();
    };


    void InitList(int lev,
                  const Vector<RealVect> & pos,
                  const Vector<Real> & r,
                  const Vector<Real> & rho);


private:

    // Positions on faces
    Vector<std::array<MultiFab, AMREX_SPACEDIM>> face_coords;


    static bool use_neighbor_list;
    static bool sort_neighbor_list;

    void InitInternals(int ngrow);
    void ReadStaticParameters();


    // TODO: this might bto be used anymore:
    AmrCore * m_amr_core;
    int nghost;


};


#endif
