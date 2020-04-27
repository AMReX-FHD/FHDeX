#ifndef _particle_functions_H_
#define _particle_functions_H_

#include <AMReX.H>
#include <AMReX_Particles.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include "AMReX_NeighborParticles.H"
#include <AMReX_ArrayLim.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Utility.H>
#include "species.H"

using namespace amrex;

//#include <AMReX_ParticleInit.H>

void ComputeDryMobility(amrex::MultiFab & dryMobility, 
                        species* particleInfo, 
                        const amrex::Geometry & geom);

#endif
