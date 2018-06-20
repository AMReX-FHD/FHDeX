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

namespace amrex {

class FhdParticleContainer
    : public ParticleContainer<10, 2+3*BL_SPACEDIM>
{
public:

    using FhdParIter = ParIter<10, 2+3*BL_SPACEDIM>;

    FhdParticleContainer(const Geometry            & geom, 
                                  const DistributionMapping & dmap,
                                  const BoxArray            & ba);

    void InitParticles(iMultiFab& collisionCellMembers, iMultiFab& collisionCellLists, const Real* dxc, const int ppc, const int ppb, species particleInfo);

    void updateParticles(const Real dt, const Real* dx, const std::array<MultiFab, AMREX_SPACEDIM>& umac,
                                           std::array<MultiFab, AMREX_SPACEDIM>& umacNodal,
                                           const std::array<MultiFab, AMREX_SPACEDIM>& RealFaceCoords,
                                           const MultiFab& betaCC, //Not necessary but may use later
                                           MultiFab& betaNodal, //Not necessary but may use later
                                           const MultiFab& rho, //Not necessary but may use later
                                           std::array<MultiFab, AMREX_SPACEDIM>& source,
                                           std::array<MultiFab, AMREX_SPACEDIM>& sourceTemp,
                                           iMultiFab& collisionCellMembers, iMultiFab& collisionCellLists, const Real* dxc, const int* hiVect, const int ppc);

    void InitCollisionCells(const iMultiFab& collisionCellMembers,
                                          const iMultiFab& collisionCellLists,
                                          MultiFab& collisionPairs,
                                          MultiFab& collisionFactor, 
                                          MultiFab& cellVols, const int ppc, const species particleInfo, const Real delt);

    void CollideParticles(const iMultiFab& collisionCellMembers,
                                          const iMultiFab& collisionCellLists,
                                          MultiFab& collisionPairs,
                                          MultiFab& collisionFactor, 
                                          MultiFab& cellVols, const int ppc, const Real Neff, const Real delt);

    void EvaluateFields(const iMultiFab& collisionCellMembers,
                              const iMultiFab& collisionCellLists,
                              MultiFab& particleDensity,
                              std::array<MultiFab, 3>& particleVelocity,
                              MultiFab& particleTemperature,
                              MultiFab& cellVols, const int ppc, const Real Neff, const Real delt);

    void WriteParticlesAscii(int n);

};
    
}

#endif
