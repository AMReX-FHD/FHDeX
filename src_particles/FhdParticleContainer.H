#ifndef _FhdParticleContainer_H_
#define _FhdParticleContainer_H_

#include <AMReX_Particles.H>
#include "AMReX_NeighborParticles.H"

#include "gmres_functions.H"
#include "common_functions.H"
#include "species.H"
#include "surfaces.H"

#include "common_functions.H"
#include "common_functions_F.H"
#include "common_namespace.H"

#include "particle_functions_F.H"
#include "rng_functions_F.H"


//#include "common_namespace_declarations.H"

using namespace common;

namespace amrex {

struct RealData
{
    enum {
        vx = 0,
        vy,
        vz,
        fx,
        fy,
        fz,
        ux,
        uy,
        uz,
        mass,
        R,
        radius,
        q,
        accelFactor,
        dragFactor,
        ox,
        oy,
        oz,
        ax,
        ay,
        az,
        travelTime,
        diffAv,
        stepCount,
        multi,
        dryDiff,
        wetDiff,
        totalDiff,
        sigma,
        eepsilon,
        potential,
	p3m_radius,
        ncomps
    };
};

struct IntData
{
    enum {
        sorted,
        i,
        j,
        k,
        species,
        ncomps
    };
};


class FhdParticleContainer
    : public amrex::NeighborParticleContainer<RealData::ncomps, IntData::ncomps>
{
public:

    using FhdParIter = ParIter<RealData::ncomps, IntData::ncomps>;

    FhdParticleContainer(const Geometry            & geom, 
                                  const DistributionMapping & dmap,
                                  const BoxArray            & ba,
                                  int ncells);

    void InitParticles(species* particleInfo, const Real* dxp);
    void InitParticlesBrownian(species particleInfo);
    //void InitParticles(const IntVect& a_num_particles_per_cell);

    void updateParticles(const Real dt, const Real* dx, const std::array<MultiFab, AMREX_SPACEDIM>& umac,
                                           std::array<MultiFab, AMREX_SPACEDIM>& umacNodal,
                                           const std::array<MultiFab, AMREX_SPACEDIM>& RealFaceCoords,
                                           const MultiFab& betaCC, //Not necessary but may use later
                                           MultiFab& betaNodal, //Not necessary but may use later
                                           const MultiFab& rho, //Not necessary but may use later
                                           std::array<MultiFab, AMREX_SPACEDIM>& source,
                                           std::array<MultiFab, AMREX_SPACEDIM>& sourceTemp,
                                           iMultiFab& collisionCellMembers, iMultiFab& collisionCellLists, const Real* dxc, const int* hiVect, const int ppc);

    void computeForcesNL(const MultiFab& charge, const MultiFab& coords, const Real* dx);

    void InitCollisionCells(              MultiFab& collisionPairs,
                                          MultiFab& collisionFactor, 
                                          MultiFab& cellVols, const species particleInfo, const Real delt);

    void CollideParticles(              MultiFab& collisionPairs,
                                          MultiFab& collisionFactor, 
                                          MultiFab& cellVols, const species particleInfo, const Real delt);


    void InitializeFields(    MultiFab& particleInstant,
                              MultiFab& cellVols, const species particleInfo);

    void EvaluateStats(
                              MultiFab& particleInstant,
                              MultiFab& particleMeans,
                              MultiFab& particleVars,
                              MultiFab& cellVols, species particleInfo, const Real delt, int steps);

    void WriteParticlesAscii(std::string asciiName);

    void MoveParticlesDry(const Real dt, const Real* dxFluid, const std::array<MultiFab, AMREX_SPACEDIM>& umac, 
                                           const std::array<MultiFab, AMREX_SPACEDIM>& RealFaceCoords,
                                           std::array<MultiFab, AMREX_SPACEDIM>& source,
                                           std::array<MultiFab, AMREX_SPACEDIM>& sourceTemp,
                                           const surface* surfaceList, const int surfaceCount);

  void MoveParticlesDSMC(const Real dt, const surface* surfaceList, const int surfaceCount,Real time, int* flux);

    void MoveIons(const Real dt, const Real* dxFluid, const Real* dxE, const Geometry geomF, const std::array<MultiFab, AMREX_SPACEDIM>& umac, const std::array<MultiFab, AMREX_SPACEDIM>& efield,
                                           const std::array<MultiFab, AMREX_SPACEDIM>& RealFaceCoords,
                                           std::array<MultiFab, AMREX_SPACEDIM>& source,
                                           std::array<MultiFab, AMREX_SPACEDIM>& sourceTemp,
                                           const MultiFab& mobility,
                                           const surface* surfaceList, const int surfaceCount, int sw);

    void SpreadIons(const Real dt, const Real* dxFluid, const Real* dxE, const Geometry geomF, const std::array<MultiFab, AMREX_SPACEDIM>& umac, const std::array<MultiFab, AMREX_SPACEDIM>& efield,
                                           const MultiFab& charge,
                                           const std::array<MultiFab, AMREX_SPACEDIM>& RealFaceCoords,
                                           const MultiFab& cellCenters,
                                           std::array<MultiFab, AMREX_SPACEDIM>& source,
                                           std::array<MultiFab, AMREX_SPACEDIM>& sourceTemp,
                                           const surface* surfaceList, const int surfaceCount, int sw);

    void DoRFD(const Real dt, const Real* dxFluid, const Real* dxE, const Geometry geomF, const std::array<MultiFab, AMREX_SPACEDIM>& umac, const std::array<MultiFab, AMREX_SPACEDIM>& efield,
                                           const std::array<MultiFab, AMREX_SPACEDIM>& RealFaceCoords,
                                           const MultiFab& cellCenters,
                                           std::array<MultiFab, AMREX_SPACEDIM>& source,
                                           std::array<MultiFab, AMREX_SPACEDIM>& sourceTemp,
                                           const surface* surfaceList, const int surfaceCount, int sw);

    void collectFields(const Real dt, const Real* dxPotential, const MultiFab& RealCenterCoords, const Geometry geomF, MultiFab& charge, MultiFab& chargeTemp, MultiFab& mass, MultiFab& massTemp);


    void ReBin();

    //void MoveParticles();

    int numWrongCell();

    //Anthony Stuff
    void InitParticlesDSMC(species* particleInfo, int pL, int pR, Real tL, Real tR);
    void InitParticlesDSMCtest(species* particleInfo, int num_boxes, int pL, int pR, Real tL, Real tR);
    void getParticleMapping(std::map<int,int>& particleMap, int xCells, int num_boxes, int pL, int pR);
    void ApplyThermostat(species* particleInfo, MultiFab& cellVols, surface* surfaces, 
                          const int ns, const Real tL, const Real tR);
    void Resample(species* particleInfo, Real tL, Real tR);

protected:
    
    void UpdateFortranStructures();

    void UpdateCellVectors();

    // used to store vectors of particle indices on a cell-by-cell basis
    std::map<int, amrex::BaseFab<std::vector<int> > > m_cell_vectors;

    // primitive data for passing the above into Fortran
    std::map<int, amrex::BaseFab<int> > m_vector_size;
    std::map<int, amrex::BaseFab<int*> > m_vector_ptrs;

    // used to keep track of whether we need to rebuild the above or not
    bool m_vectors_initialized = false;
    amrex::BoxArray::RefID m_BARef;
    amrex::DistributionMapping::RefID m_DMRef;

private:

    void correctCellVectors(int old_index, int new_index, int grid, const ParticleType& p);

    struct {
        bool operator()(const ParticleType& p1, const ParticleType& p2) const
        {
//            return AMREX_D_TERM(   (p1.pos(0) - p2.pos(0))*(p1.pos(0) - p2.pos(0))   ,
//                                   + (p1.pos(1) - p2.pos(1))*(p1.pos(1) - p2.pos(1)) ,
//                                   + (p1.pos(2) - p2.pos(2))*(p1.pos(2) - p2.pos(2)) )
//                <= cut_off*cut_off;

              return true;
            }
    } CheckPair;
    
    //static constexpr Real cutoff = 1.122*sigma[0]; //NOTE that we're changing cutoff here now!
    //static constexpr Real cutoff = pow(2,1./6.)*1.8e-5; //NOTE that we're changing cutoff here now!
    //static constexpr Real min_r  = 1.e-12;

};
    
}

#endif
