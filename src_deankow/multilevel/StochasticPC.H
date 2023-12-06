#ifndef STOCHASTIC_PC_H_
#define STOCHASTIC_PC_H_

#include <AMReX_Particles.H>

struct RealIdx
{
    enum {
        xold = 0,
        yold, zold,
        ncomps
    };
};

struct IntIdx
{
    enum {
        ncomps = 0
    };
};

class StochasticPC
    : public amrex::ParticleContainer<RealIdx::ncomps, IntIdx::ncomps, 0, 0,
                                      amrex::DefaultAllocator>
{

public:

    StochasticPC (amrex::ParGDBBase* gdb)
        : amrex::ParticleContainer<RealIdx::ncomps, IntIdx::ncomps, 0, 0,
                                   amrex::DefaultAllocator>(gdb)
        {}

    StochasticPC (const amrex::Geometry            & geom,
                  const amrex::DistributionMapping & dmap,
                  const amrex::BoxArray            & ba)
        : amrex::ParticleContainer<RealIdx::ncomps, IntIdx::ncomps, 0, 0,
                                   amrex::DefaultAllocator>(geom, dmap, ba)
        {}

    void InitParticles (amrex::MultiFab& phi_fine);
    void  AddParticles (amrex::MultiFab& phi_fine, amrex::BoxArray& ba_to_exclude);
    void AdvectWithRandomWalk (int lev, amrex::Real dt);
    void RemoveParticlesNotInBA (const amrex::BoxArray& ba_to_keep);
    void RefluxCrseToFine (const amrex::BoxArray& ba_to_keep, amrex::MultiFab& phi_for_reflux);
    void RefluxFineToCrse (const amrex::BoxArray& ba_to_keep, amrex::MultiFab& phi_for_reflux);
};

#endif
