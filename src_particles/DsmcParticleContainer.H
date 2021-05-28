#ifndef _DSMCParticleContainer_H_
#define _DSMCParticleContainer_H_

//#include "IBParticleInfo.H"
//#include "IBMarkerContainerBase.H"
//#include "species.H"

#include "paramPlane.H"
#include <AMReX_NeighborParticles.H>
//#include "paramplane_functions_F.H"

//#include "particle_functions_F.H"
//#include "particle_functions.H"

using namespace std;

// IBM => Immmersed Boundary Marker
struct Dsmc_realData {
    //Analogous to particle realData (p.m_data)
    enum {
        radius = 0,
        velx,
        vely,
        velz,
        boostx,
        boosty,
        boostz,
        mass,
        count    // Awesome little trick! (only works if first field is 0)
    };

    static Vector<std::string> names() {
        return Vector<std::string> {
            "radius",
            "velx",
            "vely",
            "velz",
            "boostx",
            "boosty",
            "boostz",
            "mass"
        };
    };
};



struct Dsmc_intData {
    //Analogous to particle intData (p.m_data)
    enum {
        sorted,
        i,
        j,
        k,
        species,
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

class DsmcParticleContainer
    : public amrex::NeighborParticleContainer<Dsmc_realData::count, Dsmc_intData::count>
{
public:

    using DsmcParIter = ParIter<Dsmc_realData::count, Dsmc_intData::count>;

    DsmcParticleContainer(const Geometry            & geom, 
                                  const DistributionMapping & dmap,
                                  const BoxArray            & ba,
                                  int ncells);

protected:
    // used to store vectors of particle indices on a cell-by-cell basis
//    std::map<int, amrex::BaseFab<std::vector<int> > > m_cell_vectors;
    std::map<int, std::vector<std::vector<int> > > m_cell_vectors;


};


#endif
