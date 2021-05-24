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
  double T;
  double q;
  double n0;
  double mfp;
  double P;
  double Neff;
  double cp;
  double sigma;
  // Simulated particle volume for evaluating radial distribution function
  double partVol
  // Wall interactions
  // Add to common (to be read as inputs)
  // Perhaps these already accoutned for
  double alpha_wall // resitution with wall
  double beta_wall // friction with wall (for later)
  
} dsmcSpecies;

// Properties between species
// Includes same species
typedef struct {

	double alpha_pp
	double csx
	
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

    void InitParticles();
	 // Initialize Collision cells here?
	 void InitCollisionCells();

	 // Some clarity here for myself would be good
	 // maybe do not want to combine with realData because variable with runs
    dsmcSpecies properties[MAX_SPECIES];
    dsmcInterSpecies propBetweenSpecies[MAX_SPECIES*(MAX_SPECIES+1)*0.5]
    // ma
    // Why are these not all ints?
    // use long in case there's a lot of particles
    long realParticles; // likely unneeded
    long simParticles;
    long totalCollisionCells;
    Real domainVol;

	 int getSpeciesIndex(int species1, int species2);
	 
	 // Evaluating radial distribution functino
	 Real g0_Ma_Ahmadi(Real solidFrac);
	 const Real chi_max = 1000.0; // value is arbitrary (change as needed or maybe as input?)
	 
	 

protected:
    // used to store vectors of particle indices on a cell-by-cell basis
//    std::map<int, amrex::BaseFab<std::vector<int> > > m_cell_vectors;
    std::map<int, std::vector<std::vector<int> > > m_cell_vectors;
    // map< ???, 2D array of particles in each cell>
    //


};


#endif
