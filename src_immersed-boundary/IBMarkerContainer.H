#ifndef _IBMarkerContainer_H_
#define _IBMarkerContainer_H_


#include <AMReX.H>
#include <AMReX_Vector.H>
#include <AMReX_RealVect.H>
#include <AMReX_AmrParGDB.H>
#include <AMReX_Particles.H>
#include <AMReX_Periodicity.H>
#include <AMReX_NeighborParticles.H>


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
        count     // Awesome little trick!
    };

    static Vector<std::string> names() {
        return Vector<std::string> {
            "velx",
            "vely",
            "velz",
            "forcex",
            "forcey",
            "forcez"
        };
    };
};



struct IBP_intData {
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



class IBMarIter final 
    : public amrex::ParIter<IBP_realData::count, IBP_intData::count, 0, 0>
{

public:

    using ParIter<IBP_realData::count, IBP_intData::count, 0, 0>::ParIter;
};



class IBMarkerContainer final
    : public amrex::NeighborParticleContainer<IBP_realData::count, IBP_intData::count>
{

public:

    using NeighborParticleContainer<IBP_realData::count, IBP_intData::count>
        ::NeighborParticleContainer;

    IBMarkerContainer(AmrCore * amr_core, int n_nbhd);
    IBMarkerContainer(const Geometry & geom, const DistributionMapping & dmap,
                      const BoxArray & ba, int n_nbhd);

    virtual ~IBMarkerContainer() {};


    // (ID, initial CPU) tuple: unique to each particle
    using MarkerIndex = std::pair<int, int>;
 




    // Get number of particles
    int NumberOfMarkers(IBParIter & pti){ return pti.GetArrayOfStructs().size();};



}





#endif
