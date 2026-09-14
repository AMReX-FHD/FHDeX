#ifndef _StochMassFlux_H_
#define _StochMassFlux_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Vector.H>

#include "common_functions.H"

using namespace amrex;

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

    // create weighted sum of stage RNGs and store in stoch_W_fc_weighted
    void weightMassFlux(Vector< amrex::Real >);

    // scale random numbers that lie on physical boundaries appropriately
    void StochMassFluxBC();

    void StochMassFluxDiv(const MultiFab& rho,
                          const MultiFab& rhotot,
                          const std::array<MultiFab, AMREX_SPACEDIM >& sqrtLonsager_fc,
                          MultiFab& stoch_mass_fluxdiv,
                          std::array<MultiFab, AMREX_SPACEDIM >& stoch_mass_flux,
                          const amrex::Real& dt,
                          const Vector< amrex::Real >& weights,
                          int increment=0);

};

#endif
