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
        timeFrac,
        omega,
        radius,
        mass,
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
        species_change, // for concentration wall
        count
    };

    static Vector<std::string> names() {
        return Vector<std::string> {
            "sorted",
            "i",
            "j",
            "k",
            "species",
            "species_change"
        };
    };
};

typedef struct {

  int type;
  int total;
  double mass;
  double radius;
  double n0;
  double R;
  double Neff;
  double partVol;
  double part2cellVol;
  // double alpha_wall; // restitution with wall
  // double beta_wall; // friction with wall (for later)
  
} dsmcSpecies;

typedef struct {

  double alpha;
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
// DSMC Routines
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	void InitParticles(Real & dt);
	void ReInitParticles();
	void InitCollisionCells();

	void Source(const Real dt, const paramPlane* paramPlaneList, const int paramPlaneCount);

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Main DSMC routines
	void SortParticles();

	Real maxDiam;
	void CalcSelections(Real dt);
	void CollideParticles(Real dt);

	void MoveParticlesCPP(const Real dt, const paramPlane* paramPlaneList, const int paramPlaneCount);
	void externalForce(const Real dt);
	void updateTimeStep(const Geometry& geom, Real& dt);
	
	Real leftMean;
	int leftcnt;
	
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// MFs for DSMC and directly related quantities
	// May condense MFs in the future
	MultiFab mfselect; // store number of selections
	MultiFab mfvrmax; // store max relative speed found
	// currently a single constant but updates over time
	MultiFab mfphi; // needed for evaluating radial distribution function
	//MultiFab mfgrantemp;

	MultiFab mfCollisions; // keep track to ensure correct collision frequency
	int expectedCollisions[MAX_SPECIES], countedCollisions[MAX_SPECIES];

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Outputs
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	void EvaluateStats(MultiFab& cuInst,
		MultiFab& cuMeans,
		MultiFab& cuVars,
		MultiFab& primInst,
		MultiFab& primMeans,
		MultiFab& primVars,
		MultiFab& mfcvlInst,
		MultiFab& mfcvlMeans,
		MultiFab& mfQMeans,
		MultiFab& coVars,
		MultiFab& spatialCross1D,
		const int steps,
		Real time);

	void OutputParticles();
	void writePlotFile(const MultiFab& cuInst,
		const MultiFab& cuMeans,
		const MultiFab& cuVars,
		const MultiFab& primInst,
		const MultiFab& primMeans,
		const MultiFab& primVars,
		const MultiFab& coVars,
		const MultiFab& mfspatialCorr1d,
		const Geometry& geom,
		Real time,
		const int ncross,
		int steps);
   
   Real Tg0[MAX_SPECIES];
   std::vector<Real> virial;
   
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
	int NSel_spec[MAX_SPECIES*MAX_SPECIES];
	int totalNSel;

	int getSpeciesIndex(int species1, int species2);
	
	// Evaluating radial distribution function
	Real g0_Ma_Ahmadi(int species1, int species2, Real phi1, Real phi2);
	const Real chi_max = 500.0; // value is arbitrary (change as needed or maybe as input?)
	const Real phi_max = 0.643;

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Constants	
	Real pi_usr = 4.0*atan(1.0);
                      
protected:
    // used to store vectors of particle indices on a cell-by-cell basis
    std::map<int, std::vector<std::vector<int> > > m_cell_vectors[MAX_SPECIES];
};


#endif
