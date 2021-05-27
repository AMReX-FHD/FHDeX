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
  double part2cellVol;
  // double alpha_wall; // resitution with wall
  // double beta_wall; // friction with wall (for later)
  
} dsmcSpecies;

typedef struct {

  double inel;
  double csx;
  
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
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Functions
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Initialize particles and cells
	void InitParticles();
	void InitCollisionCells(); // not implemented yet

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Main DSMC routines
	void SortParticles();

	void CalcCollisionCells(Real dt);
	void CollideParticles();

	void MoveParticlesCPP(const Real dt, const paramPlane* paramPlaneList, const int paramPlaneCount);
    
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Outputs
	void EvaluateStats(MultiFab& particleInstant, MultiFab& particleMeans,
                                         MultiFab& particleVars, const Real delt, int steps);

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Scalars
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Particle geometry and material
	// Note dissipation and friction for species defined in common (alpha/restitution, beta/friction)
	dsmcSpecies properties[MAX_SPECIES];
	dsmcInterSpecies interproperties[MAX_SPECIES*MAX_SPECIES];

	// Replace these by MultiFab mfnspec which seperates by species
	long realParticles;
	long simParticles;
	
	long totalCollisionCells; // might only be helpful for statistics
	Real collisionCellVol; // assumed constant for all collision cells
	Real ocollisionCellVol;
	Real domainVol;
	Real domSize[3];
	Real dx[3];

	int getSpeciesIndex(int species1, int species2);
	
	// Evaluating radial distribution function
	Real g0_Ma_Ahmadi(int species1, int species2, Real phi1, Real phi2);
	const Real chi_max = 1000.0; // value is arbitrary (change as needed or maybe as input?)
	const Real phi_max = 0.643;

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// MFs for DSMC
	// Move within particle container to avoid continuously passing as input
	MultiFab mfselect; // needed to determine how many collisions
   MultiFab mfvrmax; // needed for calculating selections
   // currently a single constant but updates over time
   MultiFab mfnspec; // needed for calculating selections
   MultiFab mfphi; // needed for evaluating radial distribution function
	
protected:
	std::map<int, std::vector<std::vector<int> > > m_cell_vectors; //[MAX_SPECIES];

};


#endif
