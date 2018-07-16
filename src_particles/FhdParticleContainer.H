#ifndef _FhdParticleContainer_H_
#define _FhdParticleContainer_H_

#include <AMReX.H>
#include <AMReX_Particles.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ArrayLim.H>

#include "gmres_functions.H"
#include "common_functions.H"
#include "species.H"
#include "surfaces.H"

namespace amrex {

struct RealData
{
    enum {
        vx = 0,
        vy,
        vz,
        mass,
        R,
        radius,
        accelFactor,
        dragFactor,
        angularVel1,
        angularVel2,
        angularVel3,
        ncomps
    };
};

struct IntData
{
    enum {
        sorted = 0,
        species,
        ncomps
    };
};


class FhdParticleContainer
    : public amrex::ParticleContainer<RealData::ncomps, IntData::ncomps>
{
public:

    using FhdParIter = ParIter<RealData::ncomps, IntData::ncomps>;

    FhdParticleContainer(const Geometry            & geom, 
                                  const DistributionMapping & dmap,
                                  const BoxArray            & ba);

    void InitParticles(const int ppc, species particleInfo);
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

    void InitCollisionCells(              MultiFab& collisionPairs,
                                          MultiFab& collisionFactor, 
                                          MultiFab& cellVols, const species particleInfo, const Real delt);

    void CollideParticles(              MultiFab& collisionPairs,
                                          MultiFab& collisionFactor, 
                                          MultiFab& cellVols, const species particleInfo, const Real delt);

    void EvaluateFields(      MultiFab& particleMembers,
                              MultiFab& particleDensity,
                              std::array<MultiFab, 3>& particleVelocity,
                              MultiFab& particleTemperature,
                              std::array<MultiFab, 3>& particleMomentum,
                              MultiFab& particleEnergy,
                              MultiFab& cellVols, const Real Neff);

    void InitializeFields(      MultiFab& particleMembers,
                              MultiFab& particleDensity,
                              std::array<MultiFab, 3>& particleVelocity,
                              MultiFab& particleTemperature,
                              MultiFab& cellVols, const species particleInfo);

    void EvaluateStats(
                              MultiFab& particleMembers,
                              MultiFab& particleDensity,
                              std::array<MultiFab, 3>& particleVelocity,
                              MultiFab& particleTemperature,
                              std::array<MultiFab, 3>& particleMomentum,
                              MultiFab& particleEnergy,

                              MultiFab& particleMembersMean,
                              MultiFab& particleDensityMean,
                              std::array<MultiFab, 3>& particleVelocityMean,
                              MultiFab& particleTemperatureMean,
                              std::array<MultiFab, 3>& particleMomentumMean,
                              MultiFab& particleEnergyMean,

                              MultiFab& particleMembersVar,
                              MultiFab& particleDensityVar,
                              std::array<MultiFab, 3>& particleVelocityVar,
                              MultiFab& particleTemperatureVar,
                              std::array<MultiFab, 3>& particleMomentumVar,
                              MultiFab& particleEnergyVar,

                              MultiFab& particleGVar, 
                              MultiFab& particleKGCross,
                              MultiFab& particleKRhoCross, 
                              MultiFab& particleRhoGCross,
                              MultiFab& particleSpatialCross,

                              MultiFab& cellVols, species particleInfo, const Real delt, int steps);

    void WriteParticlesAscii(int n);

    void MoveParticles(const Real dt, const Real* dxFluid, const Real* ploFluid, const std::array<MultiFab, AMREX_SPACEDIM>& umac,
                                           std::array<MultiFab, AMREX_SPACEDIM>& umacNodal,
                                           const std::array<MultiFab, AMREX_SPACEDIM>& RealFaceCoords,
                                           const MultiFab& betaCC, //Not necessary but may use later
                                           MultiFab& betaNodal, //Not necessary but may use later
                                           const MultiFab& rho, //Not necessary but may use later
                                           std::array<MultiFab, AMREX_SPACEDIM>& source,
                                           std::array<MultiFab, AMREX_SPACEDIM>& sourceTemp,
                                           const surface* surfaceList, const int surfaceCount);

    void ReBin();

    void UpdateCellVectors();
    //void MoveParticles();

protected:
    
    void UpdateFortranStructures();

    // used to store vectors of particle indices on a cell-by-cell basis
    std::map<int, amrex::BaseFab<std::vector<int> > > m_cell_vectors;

    // primitive data for passing the above into Fortran
    std::map<int, amrex::BaseFab<int> > m_vector_size;
    std::map<int, amrex::BaseFab<int*> > m_vector_ptrs;

    // used to keep track of whether we need to rebuild the above or not
    bool m_vectors_initialized = false;
    amrex::BoxArray::RefID m_BARef;
    amrex::DistributionMapping::RefID m_DMRef;

};
    
}

#endif
