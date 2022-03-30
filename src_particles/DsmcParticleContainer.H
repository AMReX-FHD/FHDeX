#ifndef _DSMCParticleContainer_H_
#define _DSMCParticleContainer_H_

#include "paramPlane.H"
#include <AMReX_NeighborParticles.H>
#include <AMReX_MultiFab.H>
#include <common_namespace.H>

using namespace std;
using namespace common;

struct FHD_realData {
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
		lambda,
		radius,
		mass,
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
			"R",
			"timeFrac",
			"omega",
			"lambda",
			"radius",
			"mass",
		};
	};
};



struct FHD_intData {
	enum {
		sorted,
		i,
		j,
		k,
		species,
		newSpecies,
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
	double n0;
	double R;
	double Neff;
	double partVol;
	double part2cellVol;
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

	FhdParticleContainer(const Geometry & geom, 
		const DistributionMapping & dmap,
		const BoxArray & ba,
		int ncells);
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// DSMC Routines
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	void InitParticles(Real & dt);
	void ReInitParticles();
	void InitCollisionCells();

	void Source(const Real dt, const paramPlane* paramPlaneList, const int paramPlaneCount);
    void SourcePhonons(const Real dt, const paramPlane* paramPlaneList, const int paramPlaneCount);

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Main DSMC routines
	void SortParticles();

	void CalcSelections(Real dt);
	void CollideParticles(Real dt);
	void CollideParticles2(Real dt);

	void MoveParticlesCPP(const Real dt, const paramPlane* paramPlaneList, const int paramPlaneCount);
    void MovePhononsCPP(const Real dt, const paramPlane* paramPlaneList, const int paramPlaneCount, const int step);
	void externalForce(const Real dt);
	void updateTimeStep(const Geometry& geom, Real& dt);
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	MultiFab mfselect;
	MultiFab mfvrmax;
	MultiFab mfphi;

	MultiFab mfCollisions;
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
		
   void EvaluateStatsPhonon(MultiFab& cuInst,
		MultiFab& cuMeans,
		MultiFab& cuVars,
		const int steps,
		Real time);

	void OutputParticles();
		

  //Real Tg0[MAX_SPECIES];
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Scalars
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	dsmcSpecies properties[MAX_SPECIES];
	dsmcInterSpecies interproperties[MAX_SPECIES*MAX_SPECIES];

	long realParticles;
	long simParticles;

	long totalCollisionCells; // might only be helpful for statistics
	Real collisionCellVol; // assumed constant for all collision cells
	Real ocollisionCellVol;
	Real domainVol;
	Real domSize[3];
	//Real dx[3];

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
    std::map<int, std::vector<Gpu::ManagedVector<int> > > m_cell_vectors[MAX_SPECIES];
    //std::vector<std::map<int, std::vector<Gpu::ManagedVector<int> > >> m_cell_vectors;
    //std::vector<std::vector<Gpu::ManagedVector<int>>>  m_cell_vectors[MAX_SPECIES];
};


#endif
