#ifndef _DSMCParticleContainer_H_
#define _DSMCParticleContainer_H_

#include "paramPlane.H"
#include <AMReX_NeighborParticles.H>

#include <common_namespace.H>

using namespace std;
using namespace common;

// IBM => Immmersed Boundary Marker
struct FHD_realData {
    //Analogous to particle realData (p.m_data)
    enum {
        velx = 0,
        vely,
        velz,
        boostx,
        boosty,
        boostz,
        R,
        timeFrac,
        omega,
        count    // Awesome little trick! (only works if first field is 0)
    };

    static Vector<std::string> names() {
        return Vector<std::string> {
            "velx",
            "vely",
            "velz",
            "boostx",
            "boosty",
            "boostz",
            "R",
            "timeFrac",
            "omega",
        };
    };
};



struct FHD_intData {
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
            "sorted",
            "i",
            "j",
            "k",
            "species"
        };
    };
};

typedef struct {

  int type;
  int total;
  double mass;
  double radius;
  double R;
  double T;
  double q;
  double n0;
  double mfp;
  double P;
  double Neff;
  double cp;
  double sigma;

} dsmcSpecies;

class FhdParticleContainer
    : public amrex::NeighborParticleContainer<FHD_realData::count, FHD_intData::count>
{
public:

    using FhdParIter = ParIter<FHD_realData::count, FHD_intData::count>;

    FhdParticleContainer(const Geometry            & geom, 
                                  const DistributionMapping & dmap,
                                  const BoxArray            & ba,
                                  int ncells);

    void InitParticles();
    void SortParticles();

    void Source(const Real dt, const paramPlane* paramPlaneList, const int paramPlaneCount);

    void MoveParticlesCPP(const Real dt, const paramPlane* paramPlaneList, const int paramPlaneCount);

    void EvaluateStats(MultiFab& particleInstant, MultiFab& particleMeans,
                                         MultiFab& particleVars, const Real delt, int steps);
                                         
    void PrintCellList(int i, int j, int k);
    void PrintCellListInternal(int i, int j, int k, int index, long imap);

    dsmcSpecies properties[MAX_SPECIES];
    long realParticles;
    long simParticles;
    long totalCollisionCells;
    Real domainVol;
    Real domSize[3];
    Real dx[3];


protected:
    // used to store vectors of particle indices on a cell-by-cell basis

    std::map<int, std::vector<std::vector<int> > > m_cell_vectors[MAX_SPECIES];

};


#endif
