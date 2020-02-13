#ifndef _StochMassFlux_H_
#define _StochMassFlux_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Vector.H>

#include <common_namespace.H>

using namespace amrex;
using namespace common;

class StochMassFlux {

    int n_rngs;

    Geometry geom;

    // Vectors containing Gaussian white noise field at each stage:
    // the vector is for "n_rngs", which is the number of random number stages
    // since these are face-centered, the inner vector is over faces
    Vector<std::array< MultiFab, AMREX_SPACEDIM > >  stoch_W_fc;

    // temporary storage of linear combinations of random number stages
    // multiplied by variance...
    std::array< MultiFab, AMREX_SPACEDIM >  stoch_W_fc_weighted;

public:

    // initialize n_rngs, geom
    // build MultiFabs to hold random numbers
    StochMassFlux(BoxArray,DistributionMapping,Geometry,int);
    
    // fill stoch_W_fc with random numbers
    void fillMassStochastic();
    
    void weightMassFlux(Vector< amrex::Real >);
    
    void StochMassFluxBC();
    
    void StochMassFluxDiv(std::array< amrex::MultiFab, AMREX_SPACEDIM >&,
                          const int&, const amrex::MultiFab&,
                          const std::array< amrex::MultiFab, NUM_EDGE >&, const amrex::MultiFab&,
                          const std::array< amrex::MultiFab, NUM_EDGE >&,
                          const amrex::Vector< amrex::Real >&, const amrex::Real&);

};

#endif
