#ifndef _DSMCParticleContainer_H_
#define _DSMCParticleContainer_H_

#include "paramPlane.H"
#include <AMReX_NeighborParticles.H>
#include <AMReX_MultiFab.H>
#include <common_namespace.H>

using namespace std;
using namespace common;

// IBM => Immmersed Boundary Marker
struct FHD_realData {
    //Analogous to particle realData (p.m_data)
    // Maybe add Nselect and cross section here 
    // But needs to have a variable size dependent on number of defined species
    // NSelect is 1-1,1-2,1-3,...,2-1,2-2,...
    
    // We will eventually want to track particle positions as well to select collision
    // .. partners for smaller collision cells (dense case)
    // Can be implemented now for dilute since this will be rare
    enum {
        velx = 0, // this is to indicate start index?
        vely,
        velz,
        boostx,
        boosty,
        boostz,
        R,
        count
    };

    static Vector<std::string> names() {
        return Vector<std::string> {
            "velx",
            "vely",
            "velz",
            "boostx",
            "boosty",
            "boostz",
            "R"
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
  double Neff;
  double partVol;
  double alpha_wall; // resitution with wall
  double beta_wall; // friction with wall (for later)
  
} dsmcSpecies;

// Properties between species
// Includes same species
typedef struct {

	double csx;
	double alpha_pp;
	double alpha_pw;
	double beta_pp;
	double beta_pw;
	
} dsmcInterSpecies;

class FhdParticleContainer
    : public amrex::NeighborParticleContainer<FHD_realData::count, FHD_intData::count>
{
public:

    using FhdParIter = ParIter<FHD_realData::count, FHD_intData::count>;

    FhdParticleContainer(const Geometry            & geom, 
                                  const DistributionMapping & dmap,
                                  const BoxArray            & ba,
                                  int ncells);

    void InitParticles(MultiFab & vrmax);

	 void InitCollisionCells();

    void SortParticles();

    void MoveParticlesCPP(const Real dt, const paramPlane* paramPlaneList, const int paramPlaneCount);
    
    void CollideParticles(MultiFab & mfselect, MultiFab & mfphi, MultiFab & mfvrmax);

    dsmcSpecies properties[MAX_SPECIES];
    dsmcInterSpecies propBetweenSpecies[(long)std::ceil(MAX_SPECIES*(MAX_SPECIES+1)*0.5)];

    long realParticles;
    long simParticles;
    long totalCollisionCells;
    Real domainVol;
    Real domSize[3];
    Real dx[3];

	 int getSpeciesIndex(int species1, int species2);
	 
	 // Evaluating radial distribution functino
	 Real g0_Ma_Ahmadi(int species1, int species2);
	 const Real chi_max = 1000.0; // value is arbitrary (change as needed or maybe as input?)
	
protected:
	std::map<int, std::vector<std::vector<int> > > m_cell_vectors; //[MAX_SPECIES];
	MultiFab max_relative_velocity_each_cell;


};


#endif
