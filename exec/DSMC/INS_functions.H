#ifndef _INS_functions_H_
#define _INS_functions_H_

#include "common_functions.H"
#include "gmres_functions.H"
#include "FhdParticleContainer.H"

///////////////////////////
// in WritePlotFile.cpp

void WritePlotFile(int step,
                   const amrex::Real time,
                   const amrex::Geometry geom,
                   const amrex::Geometry cgeom,
                   const MultiFab& rhotot,
                   const std::array< MultiFab, AMREX_SPACEDIM >& umac,
	           const MultiFab& div,
                   const MultiFab& particleMembers,
                   const MultiFab& particleDensity,
                   const std::array< MultiFab, 3 >& particleVelocity,
                   const MultiFab& particleTemperature,
                   const MultiFab& particlePressureMean,
                   const MultiFab& particleEnergyMean,
                   const MultiFab& spatialCorrelation,
                   const MultiFab& densityVar,
                   const MultiFab& density,
                   const MultiFab& particleMembraneFlux,
                   FhdParticleContainer& particles);

///////////////////////////


void eulerStep(const MultiFab& betaCC, const MultiFab& gammaCC,
               std::array<amrex::MultiFab, NUM_EDGE>& betaEdge, 
               std::array<MultiFab, AMREX_SPACEDIM>& umacIn, 
               std::array<MultiFab, AMREX_SPACEDIM>& umacOut,
               std::array<MultiFab, AMREX_SPACEDIM>& umacNew,
               std::array<MultiFab, AMREX_SPACEDIM>& alpha,
               const Geometry geom,
               Real* dt);

///////////////////////////


#endif
