#ifndef _particle_functions_H_
#define _particle_functions_H_

#include <AMReX_Particles.H>

#include "species.H"
#include "particle_functions_F.H"
#include "common_namespace.H"

using namespace amrex;
using namespace common;

void ComputeDryMobility(MultiFab & dryMobility, species* particleInfo, const Geometry & geom);

void MultiFabPhysBCCharge(MultiFab& charge, const Geometry& geom);

void MultiFabPhysBCDomainStress(MultiFab& stress, const amrex::Geometry& geom, int dim);

void MultiFabPhysBCMacStress(MultiFab& stress, const Geometry& geom, int dim);

///////////////////////////
// in FindCoords.cpp

void FindFaceCoords(std::array<MultiFab, AMREX_SPACEDIM> & RealFaceCoords,
                    const Geometry & geom);

void FindCenterCoords(MultiFab & RealCenterCoords, const Geometry & geom);

#endif
